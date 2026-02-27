/*

	giant.cl - Bryan Little 2/2026, montgomery arithmetic by Yves Gallot

	the baby-step-giant-step algorithm in local memory

	testing shows optimal hash table size is 8192 with 4096 elements inserted.
	nvidia gpus usually expose 49152 bytes of lmem to the kernels.
	since we cannot partition the L1/lmem in OpenCL like CUDA, we have to get creative making the hash table fit.
	we store the full 64 bit key and index to global memory.
	a 32 bit fingerprint is created by mixing the hi and lo bits of the key.
	the fingerprint is stored in lmem for fast atomic setup and lookups.
	we check the full 64 bit key and index in global memory only if the fingerprint matches.
	testing on rtx 5090 shows performance is about 97% of storing the entire 64 bit keys in partitioned lmem.

*/

typedef struct {
	ulong p;
	int n;
	int k;
} factor;

typedef struct {
	ulong hadj;
	int kidx;
} kparity;

void hash_insert(__global ulong *htable, __global int *hidx, ulong val, int idx, const uint offset, __local uint *ltable){
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint mixed = lo ^ (hi * 0x9E3779B1u);
	uint pos = mixed & MASK;			// mixed % HSIZE
	while(true){
		// Slot looks free, try to claim it
		if(atomic_cmpxchg( (volatile __local uint *)&ltable[pos], 0, mixed) == 0){
			// We successfully claimed the slot and stored the 32 bit fingerprint in local mem
			uint poff = pos+offset;
			htable[poff] = val;		// store the 64 bit key in global mem
			hidx[poff] = idx;		// store the index in global mem
			return;
		}
		// Lost the race. keep probing
		pos = (pos + 1) & MASK;			// wrap around
	}
}

int hash_lookup(__global const ulong *htable, __global const int *hidx, ulong val, int *idx_out, const uint offset, __local const uint *ltable) {
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint mixed = lo ^ (hi * 0x9E3779B1u);
	uint pos = mixed & MASK;			// mixed % HSIZE
	uint start = pos;
	while(ltable[pos] != 0){			// slot is used
		if(ltable[pos] == mixed){		// matched the 32 bit fingerprint
			uint poff = pos+offset;
			if(htable[poff] == val) {	// check if the 64 bit key also matches
				*idx_out = hidx[poff];	// use the index
				return 1;
			}
		}
		pos = (pos + 1) & MASK;			// wrap around
		if(pos == start) break;			// full cycle
	}
	return 0;
}

ulong add(ulong a, ulong b, ulong p){
	ulong r;
	ulong c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}

ulong m_mul(ulong a, ulong b, ulong p, ulong q){
	ulong lo = a*b;
	ulong hi = mul_hi(a,b);
	ulong m = lo * q;
	ulong mp = mul_hi(m,p);
	ulong r = hi - mp;
	return ( hi < mp ) ? r + p : r;
}

// left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
ulong powmodsm(ulong mbase, uint exp, ulong p, ulong q, ulong one) {
	if(!exp)return one;
	if(exp==1)return mbase;
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

// left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
ulong basepowmodsm(ulong mbase, uint exp, ulong p, ulong q, ulong one) {
	if(!exp)return one;
	if(exp==1)return mbase;
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

// dual powmod, base^lid and gQ_inv^lid
// left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
void dualbasepowmodsm(ulong mbase1, ulong mbase2, uint exp, ulong p, ulong q, ulong one, ulong *res1, ulong *res2) {
	if(!exp){
		*res1 = one;
		*res2 = one;
		return;
	}
	if(exp==1){
		*res1 = mbase1;
		*res2 = mbase2;
		return;
	}
	uint curBit = 0x80000000;
	curBit >>= ( clz(exp) + 1 );
	ulong a1 = mbase1;
	ulong a2 = mbase2;
	while( curBit )	{
		a1 = m_mul(a1,a1,p,q);
		a2 = m_mul(a2,a2,p,q);
		if(exp & curBit){
#if BASE == 2
			a1 = add(a1, a1, p);		// a*2
#elif BASE == 3
			ulong b1 = add(a1, a1, p);
			a1 = add(a1, b1, p);		// a*3
#elif BASE == 5
			ulong b1 = add(a1, a1, p);
			b1 = add(b1, b1, p);
			a1 = add(a1, b1, p);		// a*5
#elif BASE > 5
			a1 = m_mul(a1,mbase1,p,q);	// a*BASE
#endif
			a2 = m_mul(a2,mbase2,p,q);
		}
		curBit >>= 1;
	}
	*res1 = a1;
	*res2 = a2;
}

__kernel __attribute__((work_group_size_hint(1024, 1, 1))) void giantfull(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const kparity * g_k,
				__global const int * g_kcount,
				__global ulong * g_htable,
				__global int * g_hidx,
				const uint start ) {

	const int group = get_group_id(0);
	const int primepos = start + group;
	const int lid = get_local_id(0);
	const int ls = get_local_size(0);
	// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=gj_start
	const ulong8 prime = g_prime[primepos];
	const uint koffset = primepos*KCOUNT;
	const uint hashoffset = group*HSIZE;
	__local uint l_htable[HSIZE];
	__local ulong l_k_hadj[KCOUNT];
	const int numk = g_kcount[primepos];

	// zero hash table
	for(int j=lid; j<HSIZE; j+=ls){
		l_htable[j] = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	ulong gj, thread_gm_step;
	dualbasepowmodsm(prime.s3, prime.s5, lid, prime.s0, prime.s1, prime.s2, &gj, &thread_gm_step);
	gj = m_mul(prime.s7, gj, prime.s0, prime.s1); // base^NMIN * base^lid 

	// do the baby steps, with local mem atomic inserts
 	for(int j=lid; j<Q; j+=ls) {
		hash_insert(g_htable, g_hidx, gj, j, hashoffset, l_htable);
		// gj * gj_inc mod P
		gj = m_mul(gj, prime.s4, prime.s0, prime.s1);
	}

	// copy k data to local cache
	if(lid < numk){
		l_k_hadj[lid] = g_k[lid + koffset].hadj;
	}
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	// giant steps
	for(int q=lid; q<M; q+=ls){
		for(int i=0; i<numk; ++i){
			ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
			int r;
			if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)) {
				int n = NMIN + q*Q + r;
				if(n <= NMAX) {
					uint f = atomic_inc(&g_primecount[2]);
					factor fac = {prime.s0, n, klist[g_k[i + koffset].kidx]};
					g_factor[f] = fac;
				}
			}
		}
		thread_gm_step = m_mul(thread_gm_step, prime.s6, prime.s0, prime.s1);
	}
}

__kernel __attribute__((work_group_size_hint(1024, 1, 1))) void giantparity(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const kparity * g_k,
				__global const int * g_kcount,
				__global ulong * g_htable,
				__global int * g_hidx,
				const int parity,
				const int start ) {

	const int group = get_group_id(0);
	const int primepos = start + group;
	const int lid = get_local_id(0);
	const int ls = get_local_size(0);
	// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=gj_start
	const ulong8 prime = g_prime[primepos];
	const uint koffset = primepos*KCOUNT;
	const uint hashoffset = group*HSIZE;
	__local uint l_htable[HSIZE];
	__local ulong l_k_hadj[KCOUNT];
	const int numk = g_kcount[primepos];
	const int nstart = (parity==3) ? NMIN+1 : NMIN;
	const int threadInc = ls<<1;
	int threadj = lid<<1;

	// zero hash table
	for(int j=lid; j<HSIZE; j+=ls){
		l_htable[j] = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// do the baby steps, with local mem atomic inserts
	int delta = (parity==3) ? 1 : 0; // start at NMIN+1 for odd
	ulong gj = basepowmodsm(prime.s3, delta+threadj, prime.s0, prime.s1, prime.s2);
	gj = m_mul(prime.s7, gj, prime.s0, prime.s1); // base^NMIN * base^lid 

 	for(; threadj < QQ; threadj+=threadInc) {
		hash_insert(g_htable, g_hidx, gj, threadj, hashoffset, l_htable);
		// gj * gj_inc mod P
		gj = m_mul(gj, prime.s4, prime.s0, prime.s1);
	}

	// copy k data to local cache
	if(lid < numk){
		l_k_hadj[lid] = g_k[lid + koffset].hadj;
	}
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	ulong thread_gm_step = powmodsm(prime.s5, lid, prime.s0, prime.s1, prime.s2);

	// giant steps
	for(int q=lid; q<MM; q+=ls){
		for(int i=0; i<numk; ++i){
			ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
			int r;
			if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)) {
				int n = nstart + q*QQ + r;
				if(n <= NMAX) {
					uint f = atomic_inc(&g_primecount[2]);
					factor fac = {prime.s0, n, klist[g_k[i + koffset].kidx]};
					g_factor[f] = fac;
				}
			}
		}
		thread_gm_step = m_mul(thread_gm_step, prime.s6, prime.s0, prime.s1);
	}
}



