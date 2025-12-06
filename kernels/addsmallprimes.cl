/*

	addsmallprimes.cl - Bryan Little 4/2025

	generate primes <= 113
*/

ulong add(ulong a, ulong b, ulong p){
	ulong r;
	ulong c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}

ulong invert(ulong p){
	ulong p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

__kernel void addsmallprimes(ulong low, ulong high, __global ulong8 *g_prime, __global uint *g_primecount){

	const uint gid = get_global_id(0);

	const uint primes[30] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113};

	if(gid > 29) return;

	ulong p = primes[gid];

	if(p < low || p >= high) return;

	ulong q = invert(p);
	ulong one = (-p) % p;
	ulong nmo = p - one;
	ulong two = add(one, one, p);

	g_prime[ atomic_inc(&g_primecount[0]) ] = (ulong8)( p, q, 0, one, two, nmo, 0, 0 );

}


