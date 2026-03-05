/*

	clearresult.cl - Bryan Little 6/2024
	
	clear prime counters and checksum

*/

__kernel void clearresult(__global uint *g_primecount, __global ulong *g_sum){

	const uint gid = get_global_id(0);

	if(gid == 0){
		g_primecount[1] = 0;	// largest kernel prime count, for overflow checking
		g_primecount[2] = 0;	// # of factors found

		g_sum[0] = 0;		// total primecount between checkpoints  
	}

}

