/*

	Common kernel functions
	Bryan Little 3/2026

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

typedef struct {
	ulong p;
	int n;
	int k;
} factor;

ulong invert(ulong p){
	ulong p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

// note: removed Nvidia asm, not needed
// The scheduler recognizes the dependency of lo and hi multiply.
// It executes one multiply and routes results to two registers.
ulong m_mul(ulong a, ulong b, ulong p, ulong q){
	ulong lo = a*b;
	ulong hi = mul_hi(a,b);
	ulong m = lo * q;
	ulong mp = mul_hi(m,p);
	ulong r = hi - mp;
	return ( hi < mp ) ? r + p : r;
}

ulong add(ulong a, ulong b, ulong p){
	ulong r;
	ulong c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}

// left to right powmod montgomerizedbase^exp mod P, with 64 bit exponent
// assume exp > 1
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

// left to right powmod mBASE^exp mod P, with 64 bit exponent
// assume exp > 1
ulong basepowmodlg(ulong mbase, ulong exp, ulong p, ulong q) {
	ulong curBit = 0x8000000000000000;
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

// left to right powmod mBASE^exp mod P, with 32 bit exponent
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

// dual left to right powmod montgomerizedbase^exp mod P, with 32 bit exponent
// mbase1 must be BASE, mbase2 can be any montgomerized number
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


