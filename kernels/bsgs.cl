/*

	giant.cl - Bryan Little 9/2025, montgomery arithmetic by Yves Gallot

*/

typedef struct {
	ulong p;
	int nc;
	int k;
} factor;

typedef struct {
	ulong val;
	int idx;
} hash_entry;

// Hash lookup
int hash_lookup(hash_entry *table, ulong val, int *idx_out, uint offset) {
	uint pos = val & MASK;
	uint start = pos;
	while(table[pos+offset].idx != -1) {		// used
		if(table[pos+offset].val == val) {
			*idx_out = table[pos+offset].idx;
			return 1;
		}
		pos = (pos + 1) & MASK;            	// wrap around with MASK
		if(pos == start) break;           	// full cycle
	}
	return 0;
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

// left to right powmod 2^exp mod P, with 32 bit exponent
ulong pow2modsm(uint exp, ulong p, ulong q, ulong two) {
	uint curBit = 0x80000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a = two;
	while( curBit ){
		a = m_mul(a, a, p, q);
		if(exp & curBit){
			a = add(a, a, p);		// base 2 we can add
		}
		curBit >>= 1;
	}
	return a;
}

__kernel void giant(	__global const ulong8 * g_prime,
			__global uint * g_primecount,
			__global factor * g_factor,
			__global const int * g_klist,
			__global hash_entry * g_htable,
			const int kcount,
			const int L,
			const int m,
			const int nmin,
			const int nmax ) {

	const uint gid = get_global_id(0);

	if(gid >= g_primecount[0]) return;

	const uint hoffset = gid*HSIZE;

	// .s0=p, .s1=q, .s2=r2, .s3=one, .s4=two, .s5=nmo
	ulong8 prime = g_prime[gid];

	// gm^(p-2) mod P
	ulong g_m_inv = 

	// giant steps
	// gm_step = montgomerized 1
	ulong gm_step = prime.s3;

//	for(int step = 0; step * m <= L; step++) {
	for(int step = 0; step <= L; step+=m) {
		for(int l=0; l<kcount; ++l){
			ulong gamma = m_mul(h_adj[l], gm_step, prime.s0, prime.s1);
			int j;
			if(hash_lookup(g_htable, gamma, &j, hoffset)) {
//				int n = nmin + step * m + j;
				int n = nmin + step + j;
				if(n <= nmax) {
					uint f = atomic_inc(&g_primecount[2]);
					int k = g_klist[l];
					if(k<0){
						k = -k;
						n = -n;
					}
					factor fac = {prime.s0, n, k};
					g_factor[f] = fac;
				}
			}
		}
		gm_step = m_mul(gm_step, g_m_inv, prime.s0, prime.s1);
	}
	
}


