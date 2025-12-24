
/* 

	verify_factor.c

	Bryan Little Dec 2025

	montgomery arithmetic by Yves Gallot

*/

#include "verify_factor.h"

uint64_t invert(uint64_t p)
{
	uint64_t p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}


uint64_t m_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t q)
{
	unsigned __int128 res;

	res  = (unsigned __int128)a * b;
	uint64_t ab0 = (uint64_t)res;
	uint64_t ab1 = res >> 64;

	uint64_t m = ab0 * q;

	res = (unsigned __int128)m * p;
	uint64_t mp = res >> 64;

	uint64_t r = ab1 - mp;

	return ( ab1 < mp ) ? r + p : r;
}


uint64_t add(uint64_t a, uint64_t b, uint64_t p)
{
	uint64_t r;

	uint64_t c = (a >= p - b) ? p : 0;

	r = a + b - c;

	return r;
}


/* Used in the prime validator
   Returns 0 only if p is composite.
   Otherwise p is a strong probable prime to base a.
 */
bool strong_prp(uint32_t base, uint64_t p, uint64_t q, uint64_t one, uint64_t pmo, uint64_t r2, int t, uint64_t exp, uint64_t curBit)
{
	/* If p is prime and p = d*2^t+1, where d is odd, then either
		1.  a^d = 1 (mod p), or
		2.  a^(d*2^s) = -1 (mod p) for some s in 0 <= s < t    */

	uint64_t a = m_mul(base,r2,p,q);  // convert base to montgomery form
	const uint64_t mbase = a;

  	/* r <-- a^d mod p, assuming d odd */
	while( curBit )
	{
		a = m_mul(a,a,p,q);

		if(exp & curBit){
			a = m_mul(a,mbase,p,q);
		}

		curBit >>= 1;
	}

	/* Clause 1. and s = 0 case for clause 2. */
	if (a == one || a == pmo){
		return true;
	}

	/* 0 < s < t cases for clause 2. */
	for (int s = 1; s < t; ++s){

		a = m_mul(a,a,p,q);

		if(a == pmo){
	    		return true;
		}
	}


	return false;
}


// prime if the number passes prp test to 7 bases.  good to 2^64
// this is very fast
bool isPrime(uint64_t p, uint64_t q, uint64_t one, uint64_t pmo, uint64_t two, uint64_t r2)
{
	const uint32_t bases[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

	if (p % 2==0)
		return false;

	int t = __builtin_ctzll( (p-1) );
	uint64_t exp = p >> t;
	uint64_t curBit = 0x8000000000000000;
	curBit >>= ( __builtin_clzll(exp) + 1 );

	for (int i = 0; i < 7; ++i){

		uint32_t base = bases[i];

		// needed for composite bases
		if (base >= p){
			base %= p;
			if (base == 0)
				continue;
		}

		if (!strong_prp(base, p, q, one, pmo, r2, t, exp, curBit))
			return false;
	}

	return true;
}


int verify_factor(uint64_t p, uint64_t k, uint32_t n, int32_t c, uint32_t base){


	uint64_t q = invert(p);
	uint64_t one = (-p) % p;
	uint64_t pmo = p - one;	
	uint64_t two = add(one, one, p);
	uint64_t t = add(two, two, p);
	for (int i = 0; i < 5; ++i)
		t = m_mul(t, t, p, q);	// 4^{2^5} = 2^64
	uint64_t r2 = t;

	uint64_t mbase = (base==2) ? two : m_mul(base,r2,p,q);

	if(!isPrime(p, q, one, pmo, two, r2)){
		return -1;
	}

	uint32_t exp = n;
	uint32_t curBit = 0x80000000;
	curBit >>= ( __builtin_clz(exp) + 1 );

	uint64_t a = mbase;

	// a = base^n mod P
	while( curBit )
	{
		a = m_mul(a,a,p,q);

		if(exp & curBit){
			a = (base==2) ? add(a,a,p) : m_mul(a,mbase,p,q);
		}

		curBit >>= 1;
	}

	// convert k to montgomery form
	uint64_t Km = m_mul(k,r2,p,q);
 
	// b = k*base^n mod P
	uint64_t b = m_mul(a,Km,p,q);

	if(b == one && c == -1){
//		printf("%" PRIu64 " is a factor of %" PRIu64 "*%u^%u-1\n",p,k,base,n);
		return 1;
	}
	else if(b == pmo && c == 1){
//		printf("%" PRIu64 " is a factor of %" PRIu64 "*%u^%u+1\n",p,k,base,n);
		return 1;
	}

	return 0;

}



