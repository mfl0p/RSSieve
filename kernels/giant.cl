/*

	giant.cl - Bryan Little 4/2026, montgomery arithmetic by Yves Gallot

	the baby-step-giant-step algorithm in local memory

	testing shows optimal hash table size is 8192 with 4096 elements inserted.
	nvidia gpus usually expose 49152 bytes of lmem to the kernels, AMD can be 32768 bytes.
	since we cannot partition the L1/lmem in OpenCL like CUDA, we have to get creative making the hash table fit.
	we store the upper 32 bits of the value and the index to global memory.
	the lower 32 bits are stored to local memory for fast atomic setup and lookup.
	we check the hi 32 bits in global memory only if the lower 32 bits match.
	if the hi bits match, we have a solution and can use the index stored in global memory.

	if the order of b mod p is less than the number of hash table entries (Q or QQ) then no giant steps are needed.

	testing on rtx 5090 shows performance is about 97% of storing the entire 64 bit values in partitioned lmem.

*/

int claim_factor_bitmap( __global uint *g_bitmap, const uint seq_idx, const uint n ){

	if(n < NMIN || n > NMAX)
		return 0;

	uint n_off = n - NMIN;
	uint word_idx = seq_idx * WORDS_PER_SEQ + (n_off >> 5);
	uint bit_mask  = 1u << (n_off & 31u);

	/*
		Atomic test-and-clear.
		old has the bitmap word value before clearing.
		If old had bit_mask set, this work-item owns the factor.
		If not, the factor was already used, not present, or duplicated.
	*/
	uint old = atomic_and(&g_bitmap[word_idx], ~bit_mask);

	return (old & bit_mask) != 0;
}

void hash_insert(__global uint *htable, __global short *hidx, ulong val, int idx, const uint offset, __local uint *ltable){
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint pos = lo & MASK; 						// lo % HSIZE
	while(true){
		if( atomic_cmpxchg(&ltable[pos], 0, lo) == 0 ){		// Slot looks free, try to claim it
			// We successfully claimed the slot and stored lo to local mem
			uint poff = pos+offset;
			htable[poff] = hi;				// store the upper 32 bits in global mem
			hidx[poff] = idx;				// store the index in global mem
			return;
		}
		// Lost the race. keep probing
		pos = (pos + 1) & MASK;					// wrap around
	}
}

int hash_lookup(__global const uint *htable, __global const short *hidx, ulong val, int *out_idx, const uint offset, __local const uint *ltable){
	uint lo = (uint)val;
	uint hi = (uint)(val >> 32);
	uint pos = lo & MASK; 						// lo % HSIZE
	uint start = pos;
	for(uint entry = ltable[pos]; entry; entry = ltable[pos]){	// slot is used
		if(entry == lo){
			uint poff = pos + offset;
			if(htable[poff] == hi){				// check if the upper 32 bits also match
				*out_idx = (int)hidx[poff];
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
				__global uint *g_bitmap,
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
	__local ulong l_idxzero;
	__local uint l_found_order;
	__local uint l_duplicates;
	const int numk = g_kcount[primepos];

	// zero hash table
	for(int j=lid; j<HSIZE; j+=ls){
		l_htable[j] = 0;
	}

	ulong gj, thread_gm_step;
	dualbasepowmodsm(prime.s3, prime.s5, lid, prime.s0, prime.s1, prime.s2, &gj, &thread_gm_step);
	gj = m_mul(prime.s7, gj, prime.s0, prime.s1); // base^NMIN * base^lid 

	if(!lid){
		l_idxzero = gj;
		l_found_order = 0xffffffffu;
		l_duplicates = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// do the baby steps, with local mem atomic inserts
	// 'active' guarantees all threads hit the barrier
	for(uint l = 0; l < Q; l += ls) {
		uint j = l + lid;
		int active = (j < Q);
		if(active && j && gj == l_idxzero) {
			atomic_min(&l_found_order, j);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if(active && j < l_found_order) {
			hash_insert(g_htable, g_hidx, gj, j, hashoffset, l_htable);
			gj = m_mul(gj, prime.s4, prime.s0, prime.s1);
		} 
	}

#if CACHEK
	// copy k data to local cache
	__local ulong l_k_hadj[KCOUNT];
	if(lid < numk) l_k_hadj[lid] = g_k[lid + koffset].hadj;
#endif
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	int order = l_found_order == 0xffffffffu ? 0 : l_found_order;

	// if order of b mod p < Q, no giant steps are needed
	if(order){
		if(lid < numk){
#if CACHEK
			ulong gamma = l_k_hadj[lid];
#else
			ulong gamma = g_k[koffset+lid].hadj;
#endif
			int r;
			if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)){
				int kidx = g_k[lid + koffset].kidx;
				int n = NMIN + r;
				while(n <= NMAX){
					if( claim_factor_bitmap(g_bitmap, kidx, n) ){ 
						uint f = atomic_inc(&g_primecount[2]);
						factor fac = {prime.s0, n, klist[kidx]};
						g_factor[f] = fac;
					}
					else atomic_inc(&l_duplicates);
					n += order;
				}
			}
		}
	}
	// giant steps
	else{
		for(int q=lid; q<M; q+=ls){
			for(int i=0; i<numk; ++i){
#if CACHEK
				ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
#else
				ulong gamma = m_mul(g_k[koffset+i].hadj, thread_gm_step, prime.s0, prime.s1);
#endif
				int r;
				if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)){
					int kidx = g_k[i + koffset].kidx;
					int n = NMIN + q*Q + r;
					if( claim_factor_bitmap(g_bitmap, kidx, n) ){
						uint f = atomic_inc(&g_primecount[2]);
						factor fac = {prime.s0, n, klist[kidx]};
						g_factor[f] = fac;
					}
					else atomic_inc(&l_duplicates);
				}
			}
			thread_gm_step = m_mul(thread_gm_step, prime.s6, prime.s0, prime.s1);
		}
	}

	// copy duplicate counter to global mem
	barrier(CLK_LOCAL_MEM_FENCE);
	if(!lid && l_duplicates) atomic_add(&g_primecount[6], l_duplicates);
}

__kernel __attribute__((work_group_size_hint(1024, 1, 1))) void giantparity(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const kparity * g_k,
				__global const int * g_kcount,
				__global uint * g_htable,
				__global short * g_hidx,
				__global uint *g_bitmap,
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
	__local ulong l_idxzero;
	__local uint l_found_order;
	__local uint l_duplicates;
	const int numk = g_kcount[primepos];
	const int nstart = (parity==3) ? NMIN+1 : NMIN;
	const int threadInc = ls<<1;
	int threadj = lid<<1;

	// zero hash table
	for(int j=lid; j<HSIZE; j+=ls){
		l_htable[j] = 0;
	}

	// do the baby steps, with local mem atomic inserts
	int delta = (parity==3) ? 1 : 0; // start at NMIN+1 for odd
	ulong gj = basepowmodsm(prime.s3, delta+threadj, prime.s0, prime.s1, prime.s2);
	gj = m_mul(prime.s7, gj, prime.s0, prime.s1); // base^NMIN * base^lid 

	if(!lid){
		l_idxzero = gj;
		l_found_order = 0xffffffffu;
		l_duplicates = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// do the baby steps, with local mem atomic inserts
	// 'active' guarantees all threads hit the barrier
	for(uint l = 0; l < QQ; l += threadInc) {
		uint j = l + threadj;
		int active = (j < QQ);
		if(active && j && gj == l_idxzero) {
			atomic_min(&l_found_order, j);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if(active && j < l_found_order) {
			hash_insert(g_htable, g_hidx, gj, j, hashoffset, l_htable);
			gj = m_mul(gj, prime.s4, prime.s0, prime.s1);
		} 
	}

#if CACHEK
	// copy k data to local cache
	__local ulong l_k_hadj[KCOUNT];
	if(lid < numk) l_k_hadj[lid] = g_k[lid + koffset].hadj;
#endif
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	int order = l_found_order == 0xffffffffu ? 0 : l_found_order;

	// if order of b mod p < QQ, no giant steps are needed
	if(order){
		if(lid < numk){
#if CACHEK
			ulong gamma = l_k_hadj[lid];
#else
			ulong gamma = g_k[koffset+lid].hadj;
#endif
			int r;
			if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)){
				int kidx = g_k[lid + koffset].kidx;
				int n = nstart + r;
				while(n <= NMAX){
					if( claim_factor_bitmap(g_bitmap, kidx, n) ){ 
						uint f = atomic_inc(&g_primecount[2]);
						factor fac = {prime.s0, n, klist[kidx]};
						g_factor[f] = fac;
					}
					else atomic_inc(&l_duplicates);
					n += order;
				}
			}
		}
	}
	// giant steps
	else{
		ulong thread_gm_step = powmodsm(prime.s5, lid, prime.s0, prime.s1, prime.s2);

		for(int q=lid; q<MM; q+=ls){
			for(int i=0; i<numk; ++i){
#if CACHEK
				ulong gamma = m_mul(l_k_hadj[i], thread_gm_step, prime.s0, prime.s1);
#else
				ulong gamma = m_mul(g_k[koffset+i].hadj, thread_gm_step, prime.s0, prime.s1);
#endif
				int r;
				if(hash_lookup(g_htable, g_hidx, gamma, &r, hashoffset, l_htable)){
					int kidx = g_k[i + koffset].kidx;
					int n = nstart + q*QQ + r;
					if( claim_factor_bitmap(g_bitmap, kidx, n) ){
						uint f = atomic_inc(&g_primecount[2]);
						factor fac = {prime.s0, n, klist[kidx]};
						g_factor[f] = fac;
					}
					else atomic_inc(&l_duplicates);
				}
			}
			thread_gm_step = m_mul(thread_gm_step, prime.s6, prime.s0, prime.s1);
		}
	}

	// copy duplicate counter to global mem
	barrier(CLK_LOCAL_MEM_FENCE);
	if(!lid && l_duplicates) atomic_add(&g_primecount[6], l_duplicates);
}



