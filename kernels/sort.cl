typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

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

__kernel void sort(	__global uint * g_primecount,
			__global const ulong8 * g_prime,
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
	// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=pmo, .s5=gQ_inv
	const ulong8 prime = g_prime[gid];
	if(!prime.s0) return;
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

	ulong gQQ_inv, gQQ_step_inc, gjj_inc;

	if(ke || ko){
		gQQ_inv = m_mul(prime.s5, prime.s5, prime.s0, prime.s1);
		gQQ_step_inc = powmodsm(gQQ_inv, 1024, prime.s0, prime.s1);
		gjj_inc = basepowmodsm(prime.s3, 2048, prime.s0, prime.s1);
	}

	if(kf){
		ulong gQ_step_inc = powmodsm(prime.s5, 1024, prime.s0, prime.s1);
		ulong gj_inc = basepowmodsm(prime.s3, 1024, prime.s0, prime.s1);
		// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=kcount_full
		g_prime_full[primepos_full] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gj_inc, prime.s5, gQ_step_inc, kf);
	}

	if(ke){
		// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=gjj_inc, .s5=gQQ_inv, .s6=gQQ_step_inc, .s7=kcount_even
		g_prime_even[primepos_even] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gjj_inc, gQQ_inv, gQQ_step_inc, ke);
	}

	if(ko){
		// .s0=p, .s1=q, .s2=one, .s3=two/montgomerized base, .s4=gjj_inc, .s5=gQQ_inv, .s6=gQQ_step_inc, .s7=kcount_odd
		g_prime_odd[primepos_odd] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gjj_inc, gQQ_inv, gQQ_step_inc, ko);
	}
}
