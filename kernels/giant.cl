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
int hash_lookup(__global const hash_entry *table, ulong val, int *idx_out, uint offset) {
	ulong hashed = fmix64(val);
	uint pos = hashed & MASK;
	uint start = pos;
	uint poff = pos+offset;
	while(table[poff].idx != -1) {		// used
		if(table[poff].hash == hashed) {
			*idx_out = table[poff].idx;
			return 1;
		}
		pos = (pos + 1) & MASK;            	// wrap around with MASK
		if(pos == start) break;           	// full cycle
		poff = pos+offset;
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


__kernel void giantparity(	__global uint * g_primecount,
				__global factor * g_factor,
				__global const ulong4 * g_prime,
				__global const hash_entry * g_htable,
				__global const hash_entry * g_htable_e,
				__global const hash_entry * g_htable_o,
				__global kdata * g_k ) {

	const uint gid = get_global_id(0);
	const uint primepos = gid/MR;
	const uint pcnt = g_primecount[0]; 
	if(primepos >= pcnt) return;
	// .s0=p, .s1=q, .s2=one, .s3=gQ_inv
	const ulong4 prime = g_prime[primepos];
	if(!prime.s0) return;
	int q = (gid%MR)*STEPS;
	const uint hashoffset = primepos*HSIZE;
	const uint adjoffset = primepos*KCOUNT;

	ulong gm_step[STEPS];
	ulong gm_step_parity[STEPS];

	gm_step[0] = powmodsm(prime.s3, q, prime.s0, prime.s1, prime.s2);
	gm_step_parity[0] = m_mul(gm_step[0], gm_step[0], prime.s0, prime.s1);
	for(int i=1; i<STEPS; i++){
		gm_step[i] = m_mul(gm_step[i-1], prime.s3, prime.s0, prime.s1);
		gm_step_parity[i] = m_mul(gm_step[i], gm_step[i], prime.s0, prime.s1);
	}

	int r;
	ulong gamma;
	uint offset = adjoffset;
	for(int l=0; l<KCOUNT; ++l){
		kdata thek = g_k[offset++];
		if(!thek.parity) return;
		if(thek.parity>1 && q>=MM) continue;

		__global const hash_entry *htable = (thek.parity==1) ? g_htable : ( (thek.parity==2) ? g_htable_e : g_htable_o );
		for(int i=0; i<STEPS; i++){
			gamma = (thek.parity==1) ? gm_step[i] : gm_step_parity[i];
			gamma = m_mul(thek.hadj, gamma, prime.s0, prime.s1);
			if(hash_lookup(htable, gamma, &r, hashoffset)) {
				int nstart = (thek.parity==3) ? NMIN+1 : NMIN;
				int theQ = (thek.parity==1) ? Q : QQ;
				int n = nstart + (q+i)*theQ + r;
				if(n <= NMAX) {
					uint f = atomic_inc(&g_primecount[2]);
					factor fac = {prime.s0, n, klist[thek.kidx]};
					g_factor[f] = fac;
				}
			}
		}
	}
}

