/*

	setup.cl - Bryan Little 2/2026, montgomery arithmetic by Yves Gallot

	setup for testing each sieve prime

*/

typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

// note: removed Nvidia asm, not needed
// The scheduler recognizes the dependency of lo and hi multiply.
// It executes one multiply and routes results to two registers.
ulong m_mul(ulong a, ulong b, ulong p, ulong q){
	ulong lo = a*b;
	ulong hi = mul_hi(a,b);
	ulong m = lo * q;
	ulong mp = mul_hi(m,p);
	ulong r = hi - mp;
	return ( hi < mp ) ? r + p : r;
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

// left to right powmod montgomerizedbase^exp mod P, with 64 bit exponent
ulong basepowmodlg(ulong mbase, ulong exp, ulong p, ulong q) {
	ulong curBit = 0x8000000000000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = mbase;
	while( curBit )	{
		a = m_mul(a,a,p,q);
		if(exp & curBit){
#if BASE == 2
			a = add(a, a, p);	// a * 2
#elif BASE == 3
			ulong b = add(a, a, p);
			a = add(a, b, p);	// a * 3
#elif BASE == 5
			ulong b = add(a, a, p);
			b = add(b, b, p);
			a = add(a, b, p);	// a * 5
#elif BASE > 5
			a = m_mul(a,mbase,p,q);	// a * BASE
#endif
		}
		curBit >>= 1;
	}
	return a;
}

// left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
ulong basepowmodsm(ulong mbase, uint exp, ulong p, ulong q) {
	uint curBit = 0x80000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = mbase;
	while( curBit )	{
		a = m_mul(a,a,p,q);
		if(exp & curBit){
#if BASE == 2
			a = add(a, a, p);	// a * 2
#elif BASE == 3
			ulong b = add(a, a, p);
			a = add(a, b, p);	// a * 3
#elif BASE == 5
			ulong b = add(a, a, p);
			b = add(b, b, p);
			a = add(a, b, p);	// a * 5
#elif BASE > 5
			a = m_mul(a,mbase,p,q);	// a * BASE
#endif
		}
		curBit >>= 1;
	}
	return a;
}

int good_pr(ulong8 prime, ulong exponent){
	ulong rg = basepowmodlg(prime.s3, exponent, prime.s0, prime.s1);
	if(rg == prime.s2){
		return 1;	// test is valid for this P
	}
	return 0;
}

// setup for power residue tests
// literal division will be converted to something faster by the compiler
// probably a mul+shift
int setup_pr(ulong8 prime, ulong pmo, ulong *expo, int *resg){
	const ulong pm = prime.s0 - 1;
	int r = 0;
	expo[r] = pm/3;
	if( (pm == expo[r]*3) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm>>2; // pm/4
	if( ((pm&3) == 0) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/5;
	if( (pm == expo[r]*5) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/7;
	if( (pm == expo[r]*7) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm>>3; // pm/8
	if( ((pm&7) == 0) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/9;
	if( (pm == expo[r]*9) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/11;
	if( (pm == expo[r]*11) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/13;
	if( (pm == expo[r]*13) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/17;
	if( (pm == expo[r]*17) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/19;
	if( (pm == expo[r]*19) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/23;
	if( (pm == expo[r]*23) && good_pr(prime, expo[r]) ) r++;
	expo[r] = pm/29;
	if( (pm == expo[r]*29) && good_pr(prime, expo[r]) ) r++;

	// always generate for quadratic test
	expo[12] = pm>>1;
	ulong rg = basepowmodlg(prime.s3, expo[12], prime.s0, prime.s1);
	*resg = (rg == prime.s2) - (rg == pmo);

	return r;
}

// prefilter check for solvability
int prefilter(ulong hk_inv, ulong p, ulong q, ulong one, ulong pmo, ulong *power, int powcnt, int rg) {

	// power residue tests
	ulong ra;
	for(int i=0; i<powcnt; ++i){
		ra = powmodlg(hk_inv, power[i], p, q);
		if( ra!=one ) return 0;		// impossible
	}

	// quadratic
	ra = powmodlg(hk_inv, power[12], p, q);

	// If g is quadratic non-residue, g^n alternates residues/nonresidues with n parity
	if(rg == -1) {
		// if a is QR, then n must be even
		if(ra == one) return 2;
		// if a is QNR, then n must be odd
		else if(ra == pmo) return 3;
		// impossible
		else return 0;
	}
	// If g is QR and a is QNR there is no solution
	if(rg == 1 && ra == pmo) {
		// impossible
		return 0;
	}
	// n is full range
	return 1;
}

__kernel void setup(	__global ulong8 * g_prime,
			__global uint * g_primecount,
			__global kdata * g_k,
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

	// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=pmo, .s5=gQ_inv
	ulong8 prime = g_prime[gid];
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
	prime.s3 = m_mul(BASE, r2, prime.s0, prime.s1); // base = BASE * 2^64
#endif
#if BASE != 2
	g_prime[gid].s3 = prime.s3; // store montgomerized base to global
#endif

	// for batch inversion
	ulong mk[KCOUNT+1];
	ulong prefix[KCOUNT+1];

	ulong gQ = basepowmodsm(prime.s3, Q, prime.s0, prime.s1);

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
	g_prime[gid].s5 = m_mul(inv_total, prev, prime.s0, prime.s1);	// store gQ_inv to global
	inv_total = m_mul(inv_total, mk[KCOUNT], prime.s0, prime.s1);

	// setup for power residue testing
	int resg;
	ulong expo[13];
	int tests = setup_pr(prime, prime.s4, expo, &resg);

	// parity type counters
	int count[4] = {0, 0, 0, 0};

	// k data array position, pack k's that pass power residue/parity testing
	uint pos = adjoffset;

	// walk backward to get each k_inv
	for(int i = KCOUNT - 1; i >= 0; --i) {
		prev = (i == 0) ? prime.s2 : prefix[i - 1];      // prime.s2 is 'one' in montgomery
		// k_inv = inv_total * prev
		ulong k_inv = m_mul(inv_total, prev, prime.s0, prime.s1);
		// update inv_total = inv_total * mk[i] for next iteration
		inv_total = m_mul(inv_total, mk[i], prime.s0, prime.s1);

		// compute hk_inv
		ulong hk_inv = (klist[i] > 0) ? m_mul(prime.s4, k_inv, prime.s0, prime.s1) : k_inv;

		// prefilter skips unsolvable cases
		int parity = prefilter(hk_inv, prime.s0, prime.s1, prime.s2, prime.s4, expo, tests, resg);
		count[parity]++;
		if(parity){
			kdata thek = { hk_inv, parity, i };
			g_k[pos++] = thek; 
		}
	}

	// array not full?  add zero marker
	if(pos-adjoffset<KCOUNT){
		kdata thek = { 0, 0, 0 };
		g_k[pos] = thek;		
	}

	// we can skip this p if all k have no solutions
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

}


