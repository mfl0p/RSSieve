/*

	clearresult.cl - Bryan Little 6/2024
	
	clear prime counters and checksum

*/


__kernel void clearresult(__global uint *g_primecount, __global ulong *g_sum){

	const uint gid = get_global_id(0);

	if(gid == 0){
		g_primecount[1] = 0;	// keep track of largest kernel prime count
		g_primecount[2] = 0;	// # of factors found
		g_primecount[3] = 0;	// flag set for power table error
		g_primecount[4] = 0;	// flag set for getsegprimes local memory overflow
		g_primecount[5] = 0;	// 
		g_primecount[7] = 0;	// 
		g_primecount[8] = 0;	//
		g_primecount[11] = 0;	//

		g_sum[0] = 0;		// total primecount between checkpoints  
	}

}

