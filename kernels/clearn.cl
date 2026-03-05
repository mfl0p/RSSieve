/*

	clearn.cl - Bryan Little 2/2026

	Clears prime counter.

*/

__kernel void clearn(__global uint *primecount){

	const uint i = get_global_id(0);

	if(i==0){
		primecount[0]=0;	// count from the prp generator
		primecount[3]=0;	// count of full range primes
		primecount[4]=0;	// count of even parity primes
		primecount[5]=0;	// count of odd parity primes
	}
}



