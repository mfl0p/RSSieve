/*

	clearn.cl - Bryan Little 2/2026

	Clears prime counter.

*/


__kernel void clearn(__global uint *primecount, __global uint *bsgs_count){

	const uint i = get_global_id(0);

	if(i==0){
		primecount[0]=0;
		primecount[6]=0;
		primecount[9]=0;
		primecount[10]=0;

		bsgs_count[0]=0;
		bsgs_count[1]=0;
		bsgs_count[2]=0;
	}


}



