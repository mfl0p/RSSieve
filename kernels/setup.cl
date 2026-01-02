/*

	setup.cl - Bryan Little 9/2025, montgomery arithmetic by Yves Gallot

	setup for testing each sieve prime, then do the baby steps of the BSGS algorithm

*/

typedef struct {
	ulong hash;
	short idx;
} hash_entry;

// finalization mix step used in MurmurHash3 128-bit to avalanche the bits and produce a well-distributed hash output
ulong fmix64(ulong k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdUL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53UL;
    k ^= k >> 33;
    return k;
}

// Simple hash insert (linear probing, power-of-2 sized table)
void hash_insert(__global hash_entry *table, ulong hashed, int idx, uint offset) {
	uint pos = hashed & MASK;	          	// faster than val % size
	uint poff = pos+offset;
	while(table[poff].idx != -1) {		// used
		if(table[poff].hash == hashed) return;	// keep the first index
		pos = (pos + 1) & MASK;            	// wrap around with MASK
		poff = pos+offset;
	}
	table[poff].hash = hashed;
	table[poff].idx = idx;
}

// r0 + 2^64 * r1 = a * b
ulong2 mul_wide(const ulong a, const ulong b){
	ulong2 r;
#ifdef __NV_CL_C_VERSION
	const uint a0 = (uint)(a), a1 = (uint)(a >> 32);
	const uint b0 = (uint)(b), b1 = (uint)(b >> 32);
	uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c3) : "r" (a1), "r" (b1));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c3) : "r" (c3));
	r.s0 = upsample(c1, c0); r.s1 = upsample(c3, c2);
#else
	r.s0 = a * b; r.s1 = mul_hi(a, b);
#endif
	return r;
}

ulong m_mul(ulong a, ulong b, ulong p, ulong q){
	ulong2 ab = mul_wide(a,b);
	ulong m = ab.s0 * q;
	ulong mp = mul_hi(m,p);
	ulong r = ab.s1 - mp;
	return ( ab.s1 < mp ) ? r + p : r;
}

ulong add(ulong a, ulong b, ulong p){
	ulong r;
	ulong c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}

// left to right powmod montgomerizedbase^exp mod P, with 64 bit exponent
ulong powmodlg(ulong mbase, ulong exp, ulong p, ulong q) {
	ulong curBit = 0x8000000000000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = mbase;
	while( curBit )	{
		a = m_mul(a,a,p,q);
		if(exp & curBit){
			a = m_mul(a,mbase,p,q);
		}
		curBit >>= 1;
	}
	return a;
}

// left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
ulong powmodsm(ulong mbase, uint exp, ulong p, ulong q) {
	uint curBit = 0x80000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = mbase;
	while( curBit )	{
		a = m_mul(a,a,p,q);
		if(exp & curBit){
			a = m_mul(a,mbase,p,q);
		}
		curBit >>= 1;
	}
	return a;
}

// left to right powmod 2^exp mod P, with 64 bit exponent
ulong pow2modlg(ulong two, ulong exp, ulong p, ulong q) {
	ulong curBit = 0x8000000000000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = two;
	while( curBit ){
		a = m_mul(a, a, p, q);
		if(exp & curBit){
			a = add(a, a, p);	// base 2 we can add
		}
		curBit >>= 1;
	}
	return a;
}

// left to right powmod 2^exp mod P, with 32 bit exponent
ulong pow2modsm(ulong two, uint exp, ulong p, ulong q) {
	uint curBit = 0x80000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = two;
	while( curBit ){
		a = m_mul(a, a, p, q);
		if(exp & curBit){
			a = add(a, a, p);	// base 2 we can add
		}
		curBit >>= 1;
	}
	return a;
}


__constant int pres[12] = { 29, 23, 19, 17, 13, 11, 9, 8, 7, 5, 4, 3 };
// prefilter check for solvability
int prefilter(ulong hk_inv, ulong p, ulong q, ulong one, ulong pmo, ulong *power, int powcnt, int rg) {

	// power residue tests, decending
	// powmod continues from previous powmod
	ulong exp, ra, ra2, previous=0;
	for(int i=0; i<powcnt; ++i){
		exp = power[i] - previous;
		ra2 = powmodlg(hk_inv, exp, p, q);
		ra = (previous) ? m_mul(ra, ra2, p, q) : ra2;
		previous = power[i];
		if( ra!=one ) return 0;		// impossible
	}

	// quadratic
	exp = power[12] - previous;
	ra2 = powmodlg(hk_inv, exp, p, q);
	ra = (previous) ? m_mul(ra, ra2, p, q) : ra2;

						// there can be a factor when:
	if(rg == -1) {
		if(ra == one) return 2;		// restrict n to even
		else if(ra == pmo) return 3;	// restrict n to odd
		else return 0;			// impossible
	}
	if(rg == 1 && ra == pmo) {
		return 0;			// impossible
	}

	return 1;				// n is full range
}

__kernel void setup(	__global ulong4 * g_prime,
			__global uint * g_primecount,
			__global hash_entry * g_htable,
			__global hash_entry * g_htable_even,
			__global hash_entry * g_htable_odd,
			__global ulong * g_hadj,
			__global char * g_parity,
			__global ulong * g_sum ) {

	const uint gid = get_global_id(0);
	const uint pcnt = g_primecount[0]; 

	if(gid >= pcnt) return;

	if(gid == 0){
		// add primecount to total primecount
		g_sum[0] += pcnt;
		// store largest kernel prime count for array bounds check
		if( pcnt > g_primecount[1] ){
			g_primecount[1] = pcnt;
		}
	}

	// .s0=p, .s1=q, .s2=one, .s3=two
	ulong4 prime = g_prime[gid];
	const ulong pmo = prime.s0 - prime.s2;	// montgomerized p-1
	const ulong pm = prime.s0 - 1;
	const uint hashoffset = gid*HSIZE;
	const uint adjoffset = gid*KCOUNT;

	ulong four = add(prime.s3, prime.s3, prime.s0);

	// setup r2
	ulong r2 = m_mul(four, four, prime.s0, prime.s1);
	r2 = m_mul(r2, r2, prime.s0, prime.s1);
	r2 = m_mul(r2, r2, prime.s0, prime.s1);
	r2 = m_mul(r2, r2, prime.s0, prime.s1);
	r2 = m_mul(r2, r2, prime.s0, prime.s1);	// 4^{2^5} = 2^64

	// set .s3 to montgomery form of BASE. no setup is needed for base 2. 
#if BASE == 3
	prime.s3 = add(prime.s2, prime.s3, prime.s0); // base = montgomerized 1+2
#elif BASE == 5
	prime.s3 = add(four, prime.s2, prime.s0); // base = montgomerized 4+1
#elif BASE > 5
	prime.s3 = m_mul(BASE, r2, prime.s0, prime.s1); // base = montgomerized BASE
#endif

	// for batch inversion
	ulong mk[KCOUNT+1];
	ulong prefix[KCOUNT+1];

#if BASE == 2
	ulong gQ = pow2modsm(prime.s3, Q, prime.s0, prime.s1);
#else
	ulong gQ = powmodsm(prime.s3, Q, prime.s0, prime.s1);
#endif

	mk[KCOUNT] = gQ;	// use last index for gQ_inverse

	// compute mk[i] = montgomerized k
	for (int i = 0; i < KCOUNT; ++i) {
		ulong aki = (klist[i] < 0) ? -klist[i] : klist[i];
		mk[i] = m_mul(aki, r2, prime.s0, prime.s1);
	}

	// prefix products: prefix[i] = mk[0]*...*mk[i]
	prefix[0] = mk[0];
	for (int i = 1; i < KCOUNT+1; ++i) {
		prefix[i] = m_mul(prefix[i-1], mk[i], prime.s0, prime.s1);
	}

	// invert total product once: inv_total = prefix[KCOUNT-1]^(p-2)
	ulong inv_total = powmodlg(prefix[KCOUNT], prime.s0-2, prime.s0, prime.s1);

	// test if inverse exists
	ulong one = m_mul(prefix[KCOUNT], inv_total, prime.s0, prime.s1);
	if(one != prime.s2){
		// inverse doesn't exist, we can skip this 2-prp
		atomic_inc(&g_primecount[3]);
		g_prime[gid].s0 = 0;
		return;
	}

	ulong prev = prefix[KCOUNT-1];
	ulong gQ_inv = m_mul(inv_total, prev, prime.s0, prime.s1);
	inv_total = m_mul(inv_total, mk[KCOUNT], prime.s0, prime.s1);

	// parity type counters
	int count[4] = {0, 0, 0, 0};

	// setup power residue tests
	ulong expo[13];
	int r = 0;
	for(int i=0; i<12; ++i){
//		expo[r] = ( pm%pres[i] ) ? 0 : pm/pres[i];
		ulong quot = pm/pres[i];
		expo[r] = ( pm != (quot*pres[i]) ) ? 0 : quot;

		if(expo[r]){
#if BASE == 2
			ulong rg = pow2modlg(prime.s3, expo[r], prime.s0, prime.s1);
#else
			ulong rg = powmodlg(prime.s3, expo[r], prime.s0, prime.s1);
#endif
			if(rg == prime.s2){
				++r;			// test is valid for this P
			}
		}
	}

	// always generate for quadratic test
	expo[12] = pm>>1;
#if BASE == 2
	ulong rg = pow2modlg(prime.s3, expo[12], prime.s0, prime.s1);
#else
	ulong rg = powmodlg(prime.s3, expo[12], prime.s0, prime.s1);
#endif

	int resg = (rg == prime.s2) - (rg == pmo);

	// walk backward to get each k_inv
	for(int i = KCOUNT - 1; i >= 0; --i) {
		prev = (i == 0) ? prime.s2 : prefix[i - 1];      // prime.s2 is 'one' in montgomery
		// k_inv = inv_total * prev
		ulong k_inv = m_mul(inv_total, prev, prime.s0, prime.s1);
		// update inv_total = inv_total * mk[i] for next iteration
		inv_total = m_mul(inv_total, mk[i], prime.s0, prime.s1);

		// compute hk_inv
		ulong hk_inv = (klist[i] > 0) ? m_mul(pmo, k_inv, prime.s0, prime.s1) : k_inv;

		// prefilter skips unsolvable cases
		int parity = prefilter(hk_inv, prime.s0, prime.s1, prime.s2, pmo, expo, r, resg);
		g_parity[adjoffset+i] = parity;
		g_hadj[adjoffset+i] = hk_inv;
		count[parity]++;
	}

	// we can skip this p
	if(count[0] == KCOUNT){
		atomic_inc(&g_primecount[3]);
		g_prime[gid].s0 = 0;
		return;
	}

	// counters for testing
	if(count[0]){		// skipped
		atomic_add(&g_primecount[5],count[0]);
	}
	if(count[1]){		// reg
		atomic_add(&g_primecount[11],count[1]);
	}
	if(count[2]){		// even
		atomic_add(&g_primecount[7],count[2]);
	}
	if(count[3]){		// odd
		atomic_add(&g_primecount[8],count[3]);
	}

	// build baby steps, offset by even NMIN
	// BASE^NMIN mod P
#if BASE == 2
	ulong gj = pow2modsm(prime.s3, NMIN, prime.s0, prime.s1);
#else
	ulong gj = powmodsm(prime.s3, NMIN, prime.s0, prime.s1);
#endif

	for(int j = 0; j < QQ; j++) {
		if(j>=Q && !count[2] && !count[3]) break;

		ulong hashed = fmix64(gj);

		if( j < Q && count[1] ){
			hash_insert(g_htable, hashed, j, hashoffset);
		}
		if( (j&1) ){
			if( count[3] ){
				hash_insert(g_htable_odd, hashed, j-1, hashoffset);
			}
		}
		else{
			if( count[2] ){
				hash_insert(g_htable_even, hashed, j, hashoffset);
			}
		}

		// gj * BASE mod P
#if BASE == 2
		gj = add(gj, gj, prime.s0);
#else
		gj = m_mul(gj, prime.s3, prime.s0, prime.s1);
#endif
	}

	// .s0=p, .s1=q, .s2=one, .s3=gQ_inv
	g_prime[gid].s3 = gQ_inv;

}


