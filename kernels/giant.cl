/*

	giant.cl - Bryan Little 3/2026, montgomery arithmetic by Yves Gallot

	the baby-step-giant-step algorithm in local memory

	testing shows optimal hash table size is 8192 with 4096 elements inserted.
	nvidia gpus usually expose 49152 bytes of lmem to the kernels.
	since we cannot partition the L1/lmem in OpenCL like CUDA, we have to get creative making the hash table fit.
	we store the upper 32 bits of the value and the index to global memory.
	the lower 32 bits are stored to local memory for fast atomic setup and lookup.
	we check the hi 32 bits in global memory only if the lower 32 bits match.
	if the hi bits match, we have a solution and can use the index stored global memory.
	upper 32 bits are needed as we expect a large amount of collisions during the giant steps. 
	testing on rtx 5090 shows performance is about 97% of storing the entire 64 bit values in partitioned lmem.

*/

void hash_insert(__global uint *htable, __global short *hidx, ulong val, int idx, const uint offset, __local uint *ltable){
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint pos = lo & MASK; 	// lo % HSIZE
	while(true){
		// Slot looks free, try to claim it
		if(atomic_cmpxchg(&ltable[pos], 0, lo) == 0){
			// We successfully claimed the slot and stored lo to local mem
			uint poff = pos+offset;
			htable[poff] = hi;		// store the upper 32 bits in global mem
			hidx[poff] = idx;		// store the index in global mem
			return;
		}
		// Lost the race. keep probing
		pos = (pos + 1) & MASK;			// wrap around
	}
}

int hash_lookup(__global const uint *htable, __global const short *hidx, ulong val, int *idx_out, const uint offset, __local const uint *ltable) {
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint pos = lo & MASK; // lo % HSIZE
	uint start = pos;
	for(uint entry = ltable[pos]; entry; entry = ltable[pos]){	// slot is used
		if(entry == lo){ 
			uint poff = pos+offset;
			if(htable[poff] == hi) {			// check if the upper 32 bits also match
				*idx_out = hidx[poff];			// use the index
				return 1;
			}
		}
		pos = (pos + 1) & MASK;					// wrap around
		if(pos == start) break;					// full cycle
	}
	return 0;
}


__kernel __attribute__((work_group_size_hint(1024, 1, 1))) void giantfull(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const kparity * g_k,
				__global const int * g_kcount,
				__global uint * g_htable,
				__global short * g_hidx,
				const uint start ) {

	const int group = get_group_id(0);
	const int primepos = start + group;
	if(!start){	// first iteration of kernel could overflow
		if(primepos >= g_primecount[3]) return;
	}
	const int lid = get_local_id(0);
	const int ls = get_local_size(0);
	// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=gj_start
	const ulong8 prime = g_prime[primepos];
	const uint koffset = primepos*KCOUNT;
	const uint hashoffset = group*HSIZE;
	__local uint l_htable[HSIZE];
	const int numk = g_kcount[primepos];
#if CACHEK
	__local ulong l_k_hadj[KCOUNT];
#endif
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

#if CACHEK
	// copy k data to local cache
	if(lid < numk) l_k_hadj[lid] = g_k[lid + koffset].hadj;
#endif
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	// giant steps
	for(int q=lid; q<M; q+=ls){
		for(int i=0; i<numk; ++i){
#if CACHEK
			ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
#else
			ulong gamma = m_mul(g_k[koffset+i].hadj, thread_gm_step, prime.s0, prime.s1);
#endif
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
				__global uint * g_htable,
				__global short * g_hidx,
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
	const int numk = g_kcount[primepos];
	const int nstart = (parity==3) ? NMIN+1 : NMIN;
	const int threadInc = ls<<1;
	int threadj = lid<<1;
#if CACHEK
	__local ulong l_k_hadj[KCOUNT];
#endif
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

#if CACHEK
	// copy k data to local cache
	if(lid < numk) l_k_hadj[lid] = g_k[lid + koffset].hadj;
#endif
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	ulong thread_gm_step = powmodsm(prime.s5, lid, prime.s0, prime.s1, prime.s2);

	// giant steps
	for(int q=lid; q<MM; q+=ls){
		for(int i=0; i<numk; ++i){
#if CACHEK
			ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
#else
			ulong gamma = m_mul(g_k[koffset+i].hadj, thread_gm_step, prime.s0, prime.s1);
#endif
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



