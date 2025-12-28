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

void do_powmods( ulong power, ulong *previous, ulong *ra, ulong *rg, ulong hk_inv, ulong p, ulong q, ulong base ){
	ulong exp = power - *previous;
	ulong ra2 = powmodlg(hk_inv, exp, p, q);
#if BASE == 2
	ulong rg2 = pow2modlg(base, exp, p, q);
#else
	ulong rg2 = powmodlg(base, exp, p, q);
#endif
	*ra = *previous ? m_mul(*ra, ra2, p, q) : ra2;
	*rg = *previous ? m_mul(*rg, rg2, p, q) : rg2;
	*previous = power;
} 

// prefilter check for solvability
int prefilter(ulong hk_inv, ulong p, ulong q, ulong one, ulong base, ulong pmo, ulong pm) {

	ulong previous = 0;
	ulong ra, rg;

	// div should get converted to mul by inverse during compile?
	ulong expo[8], quot;
	// expo[0] = (pm%13) ? 0 : pm/13;   ...
	quot = pm/13;
	expo[0] = ( pm != (quot*13) ) ? 0 : quot;
	quot = pm/11;
	expo[1] = ( pm != (quot*11) ) ? 0 : quot;
	quot = pm/9;
	expo[2] = ( pm != (quot*9) ) ? 0 : quot;
	quot = pm>>3;
	expo[3] = ( pm != (quot<<3) ) ? 0 : quot;
	quot = pm/7;
	expo[4] = ( pm != (quot*7) ) ? 0 : quot;
	quot = pm/5;
	expo[5] = ( pm != (quot*5) ) ? 0 : quot;
	quot = pm>>2;
	expo[6] = ( pm != (quot<<2) ) ? 0 : quot;
	quot = pm/3;
	expo[7] = ( pm != (quot*3) ) ? 0 : quot;

	// power residue tests, decending tridecic to cubic
	// powmod continues from previous powmod
	for(int i=0; i<8; ++i){
		if( expo[i] ){
			do_powmods( expo[i], &previous, &ra, &rg, hk_inv, p, q, base );    
			if( rg==one && ra!=one ) return 0;	// impossible: base yields only tridecic...cubic but a is not one
		}
	}

	// quadratic
	do_powmods( pm>>1, &previous, &ra, &rg, hk_inv, p, q, base );    

/*
	if(ra == one) ls_a = 1;
	else if(ra == pmo) ls_a = -1;
	if(rg == one) ls_g = 1;
	else if(rg == pmo) ls_g = -1;
*/
	int ls_a = (ra == one) - (ra == pmo);
	int ls_g = (rg == one) - (rg == pmo);

	// Quadratic logic ...
	if (ls_g == -1) {
		if (ls_a == 1){
//			printf("even");
			return 2;   // restrict n to even
		}
		else if (ls_a == -1){
//			printf("odd");			
			return 3; // restrict n to odd
		}
		else return 0;
	}
	if (ls_g == 1 && ls_a == -1) {
		return 0; // impossible
	}

	return 1;
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

	// walk backward to get each k_inv
	for (int i = KCOUNT - 1; i >= 0; --i) {
		prev = (i == 0) ? prime.s2 : prefix[i - 1];      // prime.s2 is 'one' in montgomery
		// k_inv = inv_total * prev
		ulong k_inv = m_mul(inv_total, prev, prime.s0, prime.s1);
		// update inv_total = inv_total * mk[i] for next iteration
		inv_total = m_mul(inv_total, mk[i], prime.s0, prime.s1);

		// compute hk_inv
		ulong hk_inv = (klist[i] > 0) ? m_mul(pmo, k_inv, prime.s0, prime.s1) : k_inv;

		// prefilter skips unsolvable cases
		int parity = prefilter(hk_inv, prime.s0, prime.s1, prime.s2, prime.s3, pmo, pm);
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


