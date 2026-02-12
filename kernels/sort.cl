/*

	sort.cl - Bryan Little 2/2026, montgomery arithmetic by Yves Gallot

	sort parity/non parity k to be tested by the BSGS kernels

	setup constants for BSGS kernels

*/

typedef struct {
	ulong hadj;
	int parity;
	int kidx;
} kdata;

typedef struct {
	ulong hadj;
	int kidx;
} kparity;

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

__kernel void sort(	__global uint * g_primecount,
			__global const ulong8 * g_prime,
			__global ulong8 * g_prime_full,
			__global ulong8 * g_prime_even,
			__global ulong8 * g_prime_odd,
			__global const kdata * g_k,
			__global kparity * g_k_full,
			__global kparity * g_k_even,
			__global kparity * g_k_odd,
			__global int * g_kcount_full,
			__global int * g_kcount_even,
			__global int * g_kcount_odd ) {

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
			kparity kout = {thek.hadj, thek.kidx};
			g_k_full[kpos_full++] = kout;
			++kf;
		}
		else if(thek.parity==2){
			if(!ke){
				primepos_even = atomic_inc(&g_primecount[21]);
				kpos_even = primepos_even*KCOUNT;
			}
			kparity kout = {thek.hadj, thek.kidx};
			g_k_even[kpos_even++] = kout;
			++ke;
		}
		else if(thek.parity==3){
			if(!ko){
				primepos_odd = atomic_inc(&g_primecount[22]);
				kpos_odd = primepos_odd*KCOUNT;
			}
			kparity kout = {thek.hadj, thek.kidx};
			g_k_odd[kpos_odd++] = kout;
			++ko;
		}
	}

	ulong gj_start = basepowmodsm(prime.s3, NMIN, prime.s0, prime.s1, prime.s2);	// base^NMIN, NMIN is even
	ulong gj_inc = basepowmodsm(prime.s3, LS, prime.s0, prime.s1, prime.s2);	// base^(localsize of giant kernel)
	ulong gjj_inc = m_mul(gj_inc, gj_inc, prime.s0, prime.s1);			// base^(localsize of giant kernel * 2)
	ulong gQ_step_inc = powmodsm(prime.s5, LS, prime.s0, prime.s1);			// gQ_inv^(localsize of giant kernel)
	ulong gQQ_inv = m_mul(prime.s5, prime.s5, prime.s0, prime.s1);			// gQ_inv * gQ_inv
	ulong gQQ_step_inc = m_mul(gQ_step_inc, gQ_step_inc, prime.s0, prime.s1);	// gQQ_inv^(localsize of giant kernel)

	// full range
	if(kf){
		// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gj_inc, .s5=gQ_inv, .s6=gQ_step_inc, .s7=gj_start
		g_prime_full[primepos_full] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gj_inc, prime.s5, gQ_step_inc, gj_start);
		g_kcount_full[primepos_full] = kf;
	}
	// parity restricted even
	if(ke){
		// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gjj_inc, .s5=gQQ_inv, .s6=gQQ_step_inc, .s7=gj_start
		g_prime_even[primepos_even] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gjj_inc, gQQ_inv, gQQ_step_inc, gj_start);
		g_kcount_even[primepos_even] = ke;
	}
	// parity restricted odd
	if(ko){
		// .s0=p, .s1=q, .s2=one, .s3=montgomerized base, .s4=gjj_inc, .s5=gQQ_inv, .s6=gQQ_step_inc, .s7=gj_start
		g_prime_odd[primepos_odd] = (ulong8)(prime.s0, prime.s1, prime.s2, prime.s3, gjj_inc, gQQ_inv, gQQ_step_inc, gj_start);
		g_kcount_odd[primepos_odd] = ko;
	}
}
