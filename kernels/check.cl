/*

	check.cl - Bryan Little 4/2025, montgomery arithmetic by Yves Gallot

	Generate a checksum for boinc quorum 

*/

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

__kernel __attribute__ ((reqd_work_group_size(256, 1, 1))) void factorial_compositorial_check(	__global ulong8 * g_prime,
												__global uint * g_primecount,
												__global ulong * g_sum,
												const uint nmax ) {

	const uint gid = get_global_id(0);
	const uint lid = get_local_id(0);
	const uint pcnt = g_primecount[0];
	__local ulong sum[256];

	if(gid < pcnt){
		// .s0=p, .s1=q, .s2=r2, .s3=one, .s4=two, .s5=nmo, .s6=residue of final factorial, .s7= montgomery form of last n
		const ulong8 prime = g_prime[gid];

		sum[lid] = prime.s6 + prime.s7;

		// convert last n out of montgomery form
		uint result = (uint)m_mul(prime.s7, 1, prime.s0, prime.s1);

		// adjust result for case where nmax > pmin
		if(prime.s0 <= nmax){
			if(nmax % prime.s0 == result){
				result = nmax;
			}
		}

		if(result != nmax){
			atomic_or(&g_primecount[5], 1);
		}
	}
	else{
		sum[lid] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for(uint s = 128; s > 0; s >>= 1){
		if(lid < s){
			sum[lid] += sum[lid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if(lid == 0){
		uint index = get_group_id(0) + 1;
		g_sum[index] += sum[0];
	}

	if(gid == 0){
		// add primecount to total primecount
		g_sum[0] += pcnt;
		// store largest kernel prime count for array bounds check
		if( pcnt > g_primecount[1] ){
			g_primecount[1] = pcnt;
		}
	}

}


__kernel __attribute__ ((reqd_work_group_size(256, 1, 1))) void primorial_check(	__global ulong8 * g_prime,
											__global uint * g_primecount,
											__global ulong * g_sum ) {

	const uint gid = get_global_id(0);
	const uint lid = get_local_id(0);
	const uint pcnt = g_primecount[0];
	__local ulong sum[256];

	if(gid < pcnt){
		sum[lid] = g_prime[gid].s6;
	}
	else{
		sum[lid] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for(uint s = 128; s > 0; s >>= 1){
		if(lid < s){
			sum[lid] += sum[lid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if(lid == 0){
		uint index = get_group_id(0) + 1;
		g_sum[index] += sum[0];
	}

	if(gid == 0){
		// add primecount to total primecount
		g_sum[0] += pcnt;
		// store largest kernel prime count for array bounds check
		if( pcnt > g_primecount[1] ){
			g_primecount[1] = pcnt;
		}
	}

}


__kernel __attribute__ ((reqd_work_group_size(256, 1, 1))) void combined_check(	__global ulong8 * g_prime,
										__global uint * g_primecount,
										__global ulong * g_sum,
										const uint nmax ) {

	const uint gid = get_global_id(0);
	const uint lid = get_local_id(0);
	const uint pcnt = g_primecount[0];
	__local ulong sum[256];

	if(gid < pcnt){
		// .s0=p, .s1=q, .s2=r2, .s3=one, .s4=residue of final compositorial, .s5=nmo, .s6=residue of final factorial, .s7= montgomery form of last n
		const ulong8 prime = g_prime[gid];

		sum[lid] = prime.s4 + prime.s6 + prime.s7;

		// convert last n out of montgomery form
		uint result = (uint)m_mul(prime.s7, 1, prime.s0, prime.s1);

		// adjust result for case where nmax > pmin
		if(prime.s0 <= nmax){
			if(nmax % prime.s0 == result){
				result = nmax;
			}
		}

		if(result != nmax){
			atomic_or(&g_primecount[5], 1);
		}
	}
	else{
		sum[lid] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for(uint s = 128; s > 0; s >>= 1){
		if(lid < s){
			sum[lid] += sum[lid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if(lid == 0){
		uint index = get_group_id(0) + 1;
		g_sum[index] += sum[0];
	}

	if(gid == 0){
		// add primecount to total primecount
		g_sum[0] += pcnt;
		// store largest kernel prime count for array bounds check
		if( pcnt > g_primecount[1] ){
			g_primecount[1] = pcnt;
		}
	}

}








