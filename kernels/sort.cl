typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

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

__kernel void sort(	__global uint * g_primecount,
			__global const ulong4 * g_prime,
			__global ulong8 * g_prime_full,
			__global ulong8 * g_prime_even,
			__global ulong8 * g_prime_odd,
			__global const kdata * g_k,
			__global kdata * g_k_full,
			__global kdata * g_k_even,
			__global kdata * g_k_odd ) {

	const uint gid = get_global_id(0);
	const uint pcnt = g_primecount[0]; 
	if(gid >= pcnt) return;
	// .s0=p, .s1=q, .s2=one, .s3=gQ_inv
	const ulong4 prime = g_prime[gid];
	if(!prime.s0) return;
	uint hashoffset = gid*HSIZE;
	uint primepos_full;
	uint primepos_even;
	uint primepos_odd;
	uint kpos_full;
	uint kpos_even;
	uint kpos_odd;
	int kf=0;
	int ke=0;
	int ko=0;

	uint kpos = gid*KCOUNT;
	for(int i=0; i<KCOUNT; ++i){
		kdata thek = g_k[kpos++];
		if(!thek.parity){
			break;	// done
		}
		else if(thek.parity==1){
			if(!kf){
				primepos_full = atomic_inc(&g_primecount[20]);
				kpos_full = primepos_full*KCOUNT;
			}
			g_k_full[kpos_full++] = thek;
			++kf;
		}
		else if(thek.parity==2){
			if(!ke){
				primepos_even = atomic_inc(&g_primecount[21]);
				kpos_even = primepos_even*KCOUNT;
			}
			g_k_even[kpos_even++] = thek;
			++ke;
		}
		else if(thek.parity==3){
			if(!ko){
				primepos_odd = atomic_inc(&g_primecount[22]);
				kpos_odd = primepos_odd*KCOUNT;
			}
			g_k_odd[kpos_odd++] = thek;
			++ko;
		}
	}

	ulong gQQ_inv, gQQ_step_inc;

	if(ke || ko){
		gQQ_inv = m_mul(prime.s3, prime.s3, prime.s0, prime.s1);
		gQQ_step_inc = powmodsm(gQQ_inv, 1024, prime.s0, prime.s1, prime.s2);
	}

	if(kf){
		ulong gQ_step_inc = powmodsm(prime.s3, 1024, prime.s0, prime.s1, prime.s2);
		g_prime_full[primepos_full] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gQ_step_inc, hashoffset, kf, 0);
	}

	if(ke){
		g_prime_even[primepos_even] = (ulong8)(prime.s0, prime.s1, prime.s2, gQQ_inv, gQQ_step_inc, hashoffset, ke, 0);
	}

	if(ko){
		g_prime_odd[primepos_odd] = (ulong8)(prime.s0, prime.s1, prime.s2, gQQ_inv, gQQ_step_inc, hashoffset, ko, 0);
	}
}
