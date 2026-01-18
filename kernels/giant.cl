/*

	giant.cl - Bryan Little 9/2025, montgomery arithmetic by Yves Gallot

	the giant steps of the BSGS algorithm

*/

typedef struct {
	ulong p;
	int n;
	int k;
} factor;

typedef struct {
	ulong hash;
	short idx;
} hash_entry;

typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

// finalization mix step used in MurmurHash3 128-bit to avalanche the bits and produce a well-distributed hash output
ulong fmix64(ulong k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdUL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53UL;
    k ^= k >> 33;
    return k;
}

// Hash lookup
int hash_lookup(__local const hash_entry *table, ulong val, int *idx_out) {
	ulong hashed = fmix64(val);
	uint pos = hashed & MASK;
	uint start = pos;
	while(table[pos].idx != -1) {		// used
		if(table[pos].hash == hashed) {
			*idx_out = table[pos].idx;
			return 1;
		}
		pos = (pos + 1) & MASK;            	// wrap around with MASK
		if(pos == start) break;           	// full cycle
	}
	return 0;
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


__kernel __attribute__ ((reqd_work_group_size(1024, 1, 1))) void giantparity(
				__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong8 * g_prime,
				__global const hash_entry * g_htable,
				__global const kdata * g_k,
				const int parity ) {

	const uint gid = get_global_id(0);
	const uint pcnt = (parity==1) ? g_primecount[20] : (parity==2) ? g_primecount[21] : g_primecount[22];
	const int primepos = get_group_id(0);
	const int lid = get_local_id(0);
	const int ls = get_local_size(0);
	if(primepos >= pcnt) return;
	// .s0=p, .s1=q, .s2=one, .s3=gQ_inv, .s4=gQ_step_inc, .s5=hashoffset, .s6=numk
	const ulong8 prime = g_prime[primepos];
	const uint koffset = primepos*KCOUNT;
	__local hash_entry l_htable[HSIZE];
	__local kdata l_k[KCOUNT];
	const uint hashoffset = prime.s5;
	const int numk = prime.s6;

	for(int i=lid; i<HSIZE; i+=ls){
		l_htable[i] = g_htable[i+hashoffset];
	}

	if(lid<numk){
		l_k[lid] = g_k[koffset+lid];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	ulong thread_gm_step = powmodsm(prime.s3, lid, prime.s0, prime.s1, prime.s2);

	int end = (parity==1) ? M : MM;
	for(int q=lid; q<end; q+=ls){
		for(int i=0; i<numk; ++i){
			ulong gamma = m_mul(l_k[i].hadj, thread_gm_step, prime.s0, prime.s1);
			int r;
			if(hash_lookup(l_htable, gamma, &r)) {
				int nstart = (parity==3) ? NMIN+1 : NMIN;
				int theQ = (parity==1) ? Q : QQ;
				int n = nstart + q*theQ + r;
				if(n <= NMAX) {
					uint f = atomic_inc(&g_primecount[2]);
					factor fac = {prime.s0, n, klist[l_k[i].kidx]};
					g_factor[f] = fac;
				}
			}
		}
		thread_gm_step = m_mul(thread_gm_step, prime.s4, prime.s0, prime.s1);
	}
}
