/*

	giant.cl - Bryan Little 2/2026, montgomery arithmetic by Yves Gallot

	the baby-step-giant-step algorithm in local memory

*/

typedef struct {
	ulong p;
	int n;
	int k;
} factor;

typedef struct {
	ulong val;
	int idx;
} hash_entry;

typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

// Simple hash insert (linear probing, power-of-2 sized table)
void hash_insert(__local hash_entry *table, ulong val, int idx){
	uint pos = val & MASK;
	while(true){
		// Slot looks free, try to claim it
		if(atomic_cmpxchg( (volatile __local int *)&table[pos].idx, -1, idx) == -1){
			// We successfully claimed the slot
			table[pos].val = val;
			return;
		}
		// Lost the race. keep probing
		pos = (pos + 1) & MASK;
	}
}

// Hash lookup
int hash_lookup(__local const hash_entry *table, ulong val, int *idx_out) {
	uint pos = val & MASK;
	uint start = pos;
	while(table[pos].idx != -1) {		// used
		if(table[pos].val == val) {
			*idx_out = table[pos].idx;
			return 1;
		}
		pos = (pos + 1) & MASK;            	// wrap around with MASK
		if(pos == start) break;           	// full cycle
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

__kernel __attribute__ ((reqd_work_group_size(1024, 1, 1))) void giantparity(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const kdata * g_k,
				const int parity ) {

	const uint gid = get_global_id(0);
	const uint pcnt = (parity==1) ? g_primecount[20] : (parity==2) ? g_primecount[21] : g_primecount[22];
	const int primepos = get_group_id(0);
	const int lid = get_local_id(0);
	const int ls = get_local_size(0);
	if(primepos >= pcnt) return;
	// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=kcount_parity
	const ulong8 prime = g_prime[primepos];
	const uint koffset = primepos*KCOUNT;
	__local hash_entry l_htable[HSIZE];
	__local kdata l_k[KCOUNT];
	const int numk = prime.s7;
	const int nstart = (parity==3) ? NMIN+1 : NMIN;
	const int theQ = (parity==1) ? Q : QQ;
	const int theInc = (parity==1) ? 1 : 2;
	const int threadInc = ls*theInc;
	int threadj = lid*theInc;

	// zero hash table
	for(int j=lid; j<HSIZE; j+=ls){
		l_htable[j].idx = -1;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// do the baby steps, with local mem atomic inserts
	ulong gj = basepowmodsm(prime.s3, nstart+threadj, prime.s0, prime.s1, prime.s2);

 	for(; threadj < theQ; threadj+=threadInc) {
		hash_insert(l_htable, gj, threadj);
		// gj * gj_inc mod P
		gj = m_mul(gj, prime.s4, prime.s0, prime.s1);
	}

	// copy k list to local cache
	if(lid<numk){
		l_k[lid] = g_k[koffset+lid];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	ulong thread_gm_step = powmodsm(prime.s5, lid, prime.s0, prime.s1, prime.s2);

	// giant steps
	int end = (parity==1) ? M : MM;
	for(int q=lid; q<end; q+=ls){
		for(int i=0; i<numk; ++i){
			ulong gamma = m_mul(l_k[i].hadj, thread_gm_step, prime.s0, prime.s1);
			int r;
			if(hash_lookup(l_htable, gamma, &r)) {
				int n = nstart + q*theQ + r;
				if(n <= NMAX) {
					uint f = atomic_inc(&g_primecount[2]);
					factor fac = {prime.s0, n, klist[l_k[i].kidx]};
					g_factor[f] = fac;
				}
			}
		}
		thread_gm_step = m_mul(thread_gm_step, prime.s6, prime.s0, prime.s1);
	}
}
