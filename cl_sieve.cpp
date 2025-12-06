/*
	THIS IS AN INCOMPLETE VERSION THAT ONLY SUPPORTS BASE 2

	RSSieve
	Bryan Little, Dec 2025
	
	with contributions by Geoffrey Reynolds and Yves Gallot

	Required minimum OpenCL version is 1.1
	CL_TARGET_OPENCL_VERSION 110 in simpleCL.h

	Search limits:  P up to 2^64 and N up to 2^31

*/

#include <unistd.h>
#include <cinttypes>
#include <math.h>
#include <ctime>

#include "boinc_api.h"
#include "boinc_opencl.h"
#include "simpleCL.h"

#include "clearn.h"
#include "clearresult.h"
#include "getsegprimes.h"
#include "addsmallprimes.h"
#include "setup.h"
#include "init.h"
#include "giant.h"

#include "putil.h"
#include "cl_sieve.h"
#include "verify_factor.h"

#define MAX_SEQUENCES 200

#define RESULTS_FILENAME "factors.txt"
#define STATE_FILENAME_A "stateA.ckp"
#define STATE_FILENAME_B "stateB.ckp"

void handle_trickle_up(workStatus & st){
	if(boinc_is_standalone()) return;
	uint64_t now = (uint64_t)time(NULL);
	if( (now-st.last_trickle) > 86400 ){	// Once per day
		st.last_trickle = now;
		double progress = boinc_get_fraction_done();
		double cpu;
		boinc_wu_cpu_time(cpu);
		APP_INIT_DATA init_data;
		boinc_get_init_data(init_data);
		double run = boinc_elapsed_time() + init_data.starting_elapsed_time;
		char msg[512];
		sprintf(msg, "<trickle_up>\n"
			    "   <progress>%lf</progress>\n"
			    "   <cputime>%lf</cputime>\n"
			    "   <runtime>%lf</runtime>\n"
			    "</trickle_up>\n",
			     progress, cpu, run  );
		char variety[64];
		sprintf(variety, "rssieve_progress");
		boinc_send_trickle_up(variety, msg);
	}
}


FILE *my_fopen(const char * filename, const char * mode){
	char resolved_name[512];
	boinc_resolve_filename(filename,resolved_name,sizeof(resolved_name));
	return boinc_fopen(resolved_name,mode);
}


void cleanup( progData & pd, searchData & sd, workStatus & st ){
	sclReleaseMemObject(pd.d_factor);
	sclReleaseMemObject(pd.d_sum);
	sclReleaseMemObject(pd.d_primes);
	sclReleaseMemObject(pd.d_primecount);
	sclReleaseMemObject(pd.d_htable);
	sclReleaseMemObject(pd.d_htable_even);
	sclReleaseMemObject(pd.d_htable_odd);
	sclReleaseMemObject(pd.d_hadj);
	sclReleaseMemObject(pd.d_primes_out);
	sclReleaseMemObject(pd.d_parity);
//	sclReleaseClSoft(pd.check);
	sclReleaseClSoft(pd.clearn);
	sclReleaseClSoft(pd.clearresult);
        sclReleaseClSoft(pd.getsegprimes);
        sclReleaseClSoft(pd.addsmallprimes);
	sclReleaseClSoft(pd.init);
	sclReleaseClSoft(pd.setup);
	sclReleaseClSoft(pd.giant);
}


// using fast binary checkpoint files with checksum calculation
void write_state( workStatus & st, searchData & sd ){

	FILE * out;

	st.state_sum = st.pmin+st.pmax+st.p+st.checksum+st.primecount+st.factorcount+st.last_trickle+st.nmin+st.nmax;

        if (sd.write_state_a_next){
		if ((out = my_fopen(STATE_FILENAME_A,"wb")) == NULL)
			fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_A);
	}
	else{
                if ((out = my_fopen(STATE_FILENAME_B,"wb")) == NULL)
                        fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_B);
        }

	if(out != NULL){

		if( fwrite(&st, sizeof(workStatus), 1, out) != 1 ){
			fprintf(stderr,"Cannot write checkpoint to file. Continuing...\n");
			// Attempt to close, even though we failed to write
			fclose(out);
		}
		else{
			// If state file is closed OK, write to the other state file
			// next time around
			if (fclose(out) == 0) 
				sd.write_state_a_next = !sd.write_state_a_next; 
		}
	}
}


int read_state( workStatus & st, searchData & sd ){

	FILE * in;
	bool good_state_a = true;
	bool good_state_b = true;
	workStatus stat_a, stat_b;

        // Attempt to read state file A
	if ((in = my_fopen(STATE_FILENAME_A,"rb")) == NULL){
		good_state_a = false;
        }
	else{
		if( fread(&stat_a, sizeof(workStatus), 1, in) != 1 ){
			fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_A);
			printf("Cannot parse %s !!!\n",STATE_FILENAME_A);
			good_state_a = false;
		}
		else if(stat_a.pmin != st.pmin || stat_a.pmax != st.pmax || stat_a.nmin != st.nmin || stat_a.nmax != st.nmax
			|| stat_a.factorial != st.factorial || stat_a.primorial != st.primorial || stat_a.compositorial != st.compositorial){
			fprintf(stderr,"Invalid checkpoint file %s !!!\n",STATE_FILENAME_A);
			printf("Invalid checkpoint file %s !!!\n",STATE_FILENAME_A);
			good_state_a = false;
		}
		else{
			uint64_t state_sum = stat_a.pmin+stat_a.pmax+stat_a.p+stat_a.checksum+stat_a.primecount+stat_a.factorcount
						+stat_a.last_trickle+stat_a.nmin+stat_a.nmax;
			if(state_sum != stat_a.state_sum){
				fprintf(stderr,"Checksum error in %s !!!\n",STATE_FILENAME_A);
				printf("Checksum error in %s !!!\n",STATE_FILENAME_A);
				good_state_a = false;
			}
		}
		fclose(in);
	}

        // Attempt to read state file B
	if ((in = my_fopen(STATE_FILENAME_B,"rb")) == NULL){
		good_state_b = false;
        }
	else{
		if( fread(&stat_b, sizeof(workStatus), 1, in) != 1 ){
			fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_B);
			printf("Cannot parse %s !!!\n",STATE_FILENAME_B);
			good_state_b = false;
		}
		else if(stat_b.pmin != st.pmin || stat_b.pmax != st.pmax || stat_b.nmin != st.nmin || stat_b.nmax != st.nmax
			|| stat_b.factorial != st.factorial || stat_b.primorial != st.primorial || stat_b.compositorial != st.compositorial){
			fprintf(stderr,"Invalid checkpoint file %s !!!\n",STATE_FILENAME_B);
			printf("Invalid checkpoint file %s !!!\n",STATE_FILENAME_B);
			good_state_b = false;
		}
		else{
			uint64_t state_sum = stat_b.pmin+stat_b.pmax+stat_b.p+stat_b.checksum+stat_b.primecount+stat_b.factorcount
						+stat_b.last_trickle+stat_b.nmin+stat_b.nmax;
			if(state_sum != stat_b.state_sum){
				fprintf(stderr,"Checksum error in %s !!!\n",STATE_FILENAME_B);
				printf("Checksum error in %s !!!\n",STATE_FILENAME_B);
				good_state_b = false;
			}
		}
		fclose(in);
	}

        // If both state files are OK, check which is the most recent
	if (good_state_a && good_state_b)
	{
		if (stat_a.p > stat_b.p)
			good_state_b = false;
		else
			good_state_a = false;
	}

        // Use data from the most recent state file
	if (good_state_a && !good_state_b)
	{
		memcpy(&st, &stat_a, sizeof(workStatus));
		sd.write_state_a_next = false;
		if(boinc_is_standalone()){
			printf("Resuming from checkpoint in %s\n",STATE_FILENAME_A);
		}
		return 1;
	}
        if (good_state_b && !good_state_a)
        {
		memcpy(&st, &stat_b, sizeof(workStatus));
		sd.write_state_a_next = true;
		if(boinc_is_standalone()){
			printf("Resuming from checkpoint in %s\n",STATE_FILENAME_B);
		}
		return 1;
        }

	// If we got here, neither state file was good
	return 0;
}


void checkpoint( workStatus & st, searchData & sd ){
	handle_trickle_up( st );
	write_state( st, sd );
	if(boinc_is_standalone()){
		printf("Checkpoint, current p: %" PRIu64 "\n", st.p);
	}
	boinc_checkpoint_completed();
}


// sleep CPU thread while waiting on the specified event to complete in the command queue
// using critical sections to prevent BOINC from shutting down the program while kernels are running on the GPU
void waitOnEvent(sclHard hardware, cl_event event){

	cl_int err;
	cl_int info;
#ifdef _WIN32
#else
	struct timespec sleep_time;
	sleep_time.tv_sec = 0;
	sleep_time.tv_nsec = 1000000;	// 1ms
#endif

	boinc_begin_critical_section();

	err = clFlush(hardware.queue);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clFlush\n" );
		fprintf(stderr, "ERROR: clFlush\n" );
		sclPrintErrorFlags( err );
       	}

	while(true){

#ifdef _WIN32
		Sleep(1);
#else
		nanosleep(&sleep_time,NULL);
#endif

		err = clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
		if ( err != CL_SUCCESS ) {
			printf( "ERROR: clGetEventInfo\n" );
			fprintf(stderr, "ERROR: clGetEventInfo\n" );
			sclPrintErrorFlags( err );
	       	}

		if(info == CL_COMPLETE){
			err = clReleaseEvent(event);
			if ( err != CL_SUCCESS ) {
				printf( "ERROR: clReleaseEvent\n" );
				fprintf(stderr, "ERROR: clReleaseEvent\n" );
				sclPrintErrorFlags( err );
		       	}

			boinc_end_critical_section();

			return;
		}
	}
}


// queue a marker and sleep CPU thread until marker has been reached in the command queue
void sleepCPU(sclHard hardware){

	cl_event kernelsDone;
	cl_int err;
	cl_int info;
#ifdef _WIN32
#else
	struct timespec sleep_time;
	sleep_time.tv_sec = 0;
	sleep_time.tv_nsec = 1000000;	// 1ms
#endif

	boinc_begin_critical_section();

	// OpenCL v2.0
/*
	err = clEnqueueMarkerWithWaitList( hardware.queue, 0, NULL, &kernelsDone);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clEnqueueMarkerWithWaitList\n");
		fprintf(stderr, "ERROR: clEnqueueMarkerWithWaitList\n");
		sclPrintErrorFlags(err); 
	}
*/
	err = clEnqueueMarker( hardware.queue, &kernelsDone);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clEnqueueMarker\n");
		fprintf(stderr, "ERROR: clEnqueueMarker\n");
		sclPrintErrorFlags(err); 
	}

	err = clFlush(hardware.queue);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clFlush\n" );
		fprintf(stderr, "ERROR: clFlush\n" );
		sclPrintErrorFlags( err );
       	}

	while(true){

#ifdef _WIN32
		Sleep(1);
#else
		nanosleep(&sleep_time,NULL);
#endif

		err = clGetEventInfo(kernelsDone, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
		if ( err != CL_SUCCESS ) {
			printf( "ERROR: clGetEventInfo\n" );
			fprintf(stderr, "ERROR: clGetEventInfo\n" );
			sclPrintErrorFlags( err );
	       	}

		if(info == CL_COMPLETE){
			err = clReleaseEvent(kernelsDone);
			if ( err != CL_SUCCESS ) {
				printf( "ERROR: clReleaseEvent\n" );
				fprintf(stderr, "ERROR: clReleaseEvent\n" );
				sclPrintErrorFlags( err );
		       	}

			boinc_end_critical_section();

			return;
		}
	}
}



// find mod 30 wheel index based on starting N
// this is used by gpu threads to iterate over the number line
void findWheelOffset(uint64_t & start, int32_t & index){

	int32_t wheel[8] = {4, 2, 4, 2, 4, 6, 2, 6};
	int32_t idx = -1;

	// find starting number using mod 6 wheel
	// N=(k*6)-1, N=(k*6)+1 ...
	// where k, k+1, k+2 ...
	uint64_t k = start / 6;
	int32_t i = 1;
	uint64_t N = (k * 6)-1;


	while( N < start || N % 5 == 0 ){
		if(i){
			i = 0;
			N += 2;
		}
		else{
			i = 1;
			N += 4;
		}
	}

	start = N;

	// find mod 30 wheel index by iterating with a mod 6 wheel until finding N divisible by 5
	// forward to find index
	while(idx < 0){

		if(i){
			N += 2;
			i = 0;
			if(N % 5 == 0){
				N -= 2;
				idx = 5;
			}

		}
		else{
			N += 4;
			i = 1;
			if(N % 5 == 0){
				N -= 4;
				idx = 7;
			}
		}
	}

	// reverse to find starting index
	while(N != start){
		--idx;
		if(idx < 0)idx=7;
		N -= wheel[idx];
	}


	index = idx;

}


int factorcompare(const void *a, const void *b) {
  	factor *factA = (factor *)a;
	factor *factB = (factor *)b;
	if(factB->p < factA->p){
		return 1;
	}
	else if(factB->p == factA->p){
		if(factB->n < factA->n){
			return 1;
		}
	}
	return -1;
}


void getResults( progData & pd, workStatus & st, searchData & sd, sclHard hardware, uint32_t * h_primecount, uint64_t * h_sum ){
	// copy checksum and total prime count to host memory, non-blocking
	sclReadNB(hardware, sd.numgroups*sizeof(uint64_t), pd.d_sum, h_sum);
	// copy prime count to host memory, blocking
	sclRead(hardware, 12*sizeof(uint32_t), pd.d_primecount, h_primecount);
	// index 0 is the gpu's total prime count
	st.primecount += h_sum[0];
	// sum blocks
//	for(uint32_t i=1; i<sd.numgroups; ++i){
//		st.checksum += h_sum[i];
//	}
	// largest kernel prime count.  used to check array bounds
	if(h_primecount[1] > sd.psize){
		fprintf(stderr,"error: gpu prime array overflow\n");
		printf("error: gpu prime array overflow\n");
		exit(EXIT_FAILURE);
	}
	// flag set if there is a gpu overflow error
	if(h_primecount[4] == 1){
		fprintf(stderr,"error: getsegprimes kernel local memory overflow\n");
		printf("error: getsegprimes kernel local memory overflow\n");
		exit(EXIT_FAILURE);
	}
	printf("skipped %u primes\n",h_primecount[3]);
	printf("skipped %u k\n",h_primecount[5]);
	printf("%u k full range n\n",h_primecount[11]);
	printf("%u k restricted to even n\n",h_primecount[7]);
	printf("%u k restricted to odd n\n",h_primecount[8]);

	uint32_t numfactors = h_primecount[2];
	if(numfactors > 0){
		if(boinc_is_standalone()){
			printf("processing %u factors on CPU\n", numfactors);
		}
		if(numfactors > sd.numresults){
			fprintf(stderr,"Error: number of results (%u) overflowed array.\n", numfactors);
			exit(EXIT_FAILURE);
		}
		factor * h_factor = (factor *)malloc(numfactors * sizeof(factor));
		if( h_factor == NULL ){
			fprintf(stderr,"malloc error: h_factor\n");
			exit(EXIT_FAILURE);
		}
		// copy factors to host memory, blocking
		sclRead(hardware, numfactors * sizeof(factor), pd.d_factor, h_factor);
		// sort results by prime size if needed
		if(numfactors > 1){
			if(boinc_is_standalone()){
				printf("sorting factors\n");
			}
			qsort(h_factor, numfactors, sizeof(factor), factorcompare);
		}
		// verify all factors on CPU
		if(boinc_is_standalone()){
			printf("Verifying factors on CPU...\n");
		}

		uint32_t prpcount = 0;

		for(uint32_t i=0; i<numfactors; ++i){
			uint64_t fp = h_factor[i].p;
			uint32_t fn = h_factor[i].n;
			uint32_t fk = (h_factor[i].k < 0) ? -h_factor[i].k : h_factor[i].k;
			int32_t fc = (h_factor[i].k < 0) ? -1 : 1;
			int32_t vres = verify_factor( fp, fk, fn, fc ); 
			if( !vres ){
				fprintf(stderr,"CPU factor verification failed!  %" PRIu64 " is not a factor of %u*2^%u%+d\n", fp, fk, fn, fc);
				printf("CPU factor verification failed!  %" PRIu64 " is not a factor of %u*2^%u%+d\n", fp, fk, fn, fc);
				exit(EXIT_FAILURE);
			}
			else if( vres == -1 ){		// Unlikely
				h_factor[i].p = 0;
				++prpcount;
			}
		}

		fprintf(stderr,"Verified %u factors. Discarded %u 2-PRP factors.\n", numfactors-prpcount, prpcount);
		if(boinc_is_standalone()){
			printf("Verified %u factors. Discarded %u 2-PRP factors.\n", numfactors-prpcount, prpcount);
		}
		// write factors to file
		FILE * resfile = my_fopen(RESULTS_FILENAME,"a");
		if( resfile == NULL ){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
		if(boinc_is_standalone()){
			printf("writing factors to %s\n", RESULTS_FILENAME);
		}
		for(uint32_t i=0; i<numfactors; ++i){
			uint64_t fp = h_factor[i].p;
			uint32_t fn = h_factor[i].n;
			uint32_t fk = (h_factor[i].k < 0) ? -h_factor[i].k : h_factor[i].k;
			int32_t fc = (h_factor[i].k < 0) ? -1 : 1;
			if(fp){
				++st.factorcount;
				if( fprintf( resfile, "%" PRIu64 " | %u*2^%u%+d\n", fp, fk, fn, fc) < 0 ){
					fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
					exit(EXIT_FAILURE);
				}
				// add the factor to checksum
				st.checksum += fk + fn + fc;
			}
		}
		fclose(resfile);
		free(h_factor);
	}
}


void setupSearch(workStatus & st, searchData & sd){

	st.p = st.pmin;

	if(st.pmin == 0 || st.pmax == 0){
		printf("\n-p and -P arguments are required\nuse -h for help\n");
		fprintf(stderr, "-p and -P arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if(st.nmin == 0 || st.nmax == 0){
		printf("\n-n and -N arguments are required\nuse -h for help\n");
		fprintf(stderr, "-n and -N arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if (st.nmin > st.nmax){
		printf("nmin <= nmax is required\nuse -h for help\n");
		fprintf(stderr, "nmin <= nmax is required\n");
		exit(EXIT_FAILURE);
	}

	if (st.pmin > st.pmax){
		printf("pmin <= pmax is required\nuse -h for help\n");
		fprintf(stderr, "pmin <= pmax is required\n");
		exit(EXIT_FAILURE);
	}

	if (st.nmin&1){
		printf("even nmin is required\nuse -h for help\n");
		fprintf(stderr, "even nmin is required\n");
		exit(EXIT_FAILURE);
	}

	// increase result buffer at low P range
	// it's still possible to overflow this with a fast GPU and large search range
	if(st.pmin < 0xFFFFFFFF){
		sd.numresults = 30000000;
		sd.numresults = 300000000;
	}

	fprintf(stderr, "Starting sieve at p: %" PRIu64 " n: %u\nStopping sieve at P: %" PRIu64 " N: %u\n", st.pmin, st.nmin, st.pmax, st.nmax);
	if(boinc_is_standalone()){
		printf("Starting sieve at p: %" PRIu64 " n: %u\nStopping sieve at P: %" PRIu64 " N: %u\n", st.pmin, st.nmin, st.pmax, st.nmax);
	}

}



void profileGPU(progData & pd, workStatus & st, searchData & sd, sclHard hardware){

	// calculate approximate chunk size based on gpu's compute units
	uint64_t start = st.p;

	uint64_t range_primes = sd.computeunits * 256;

	double C = range_primes + start / log(start);
	// x is the initial guess
	double x = C * (log(C) + log(log(C)));

	// converge on "stop" using Newton iterations
	for (int i = 0; i < 10; i++) {
		double fx = x/log(x) - C;
		double dfx = (log(x) - 1) / (log(x)*log(x));
		x = x - fx/dfx;
	}
	uint64_t stop = (uint64_t)x;
	uint64_t calc_range = stop - start;

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	// calculate prime array size based on result
	uint64_t mem_size = (uint64_t)( 1.5 * (double)range_primes );
	// make it a multiple of check kernel's local size
//	mem_size = (mem_size / pd.check.local_size[0]) * pd.check.local_size[0];	

	if(mem_size > UINT32_MAX){
		fprintf(stderr, "ERROR: mem_size too large.\n");
                printf( "ERROR: mem_size too large.\n" );
		exit(EXIT_FAILURE);
	}

	sd.range = calc_range;
	sd.psize = mem_size;

	printf("range: %llu numprimesinrange: %llu\n",sd.range,sd.psize);

}


void finalizeResults( workStatus & st ){

	char line[256];
	uint32_t lc = 0;
	FILE * resfile;

	if(st.factorcount){
		// check result file has the same number of lines as the factor count
		resfile = my_fopen(RESULTS_FILENAME,"r");

		if(resfile == NULL){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}

		while(fgets(line, sizeof(line), resfile) != NULL) {
			++lc;
		}

		fclose(resfile);

		if(lc < st.factorcount){
			fprintf(stderr,"ERROR: Missing factors in %s !!!\n",RESULTS_FILENAME);
			printf("ERROR: Missing factors in %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}

	// print checksum
	resfile = my_fopen(RESULTS_FILENAME,"a");

	if(resfile == NULL){
		fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
		exit(EXIT_FAILURE);
	}

	if(st.factorcount){
		if( fprintf( resfile, "%016" PRIX64 "\n", st.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}
	else{
		if( fprintf( resfile, "no factors\n%016" PRIX64 "\n", st.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}

	fclose(resfile);
}

// Generates a string like "__constant uint lookup[N] = {a,b,c,...};\n"
char* generate_constant_array_string(const int *arr, size_t n, const char *name) {
    // estimate needed buffer size: each number ~11 chars, plus commas/braces, plus header
    size_t buf_size = n * 12 + 128;
    char *buf = (char *)malloc(buf_size);
    if (!buf) return NULL;

    // start with __constant declaration
    int offset = snprintf(buf, buf_size, "__constant int %s[%zu] = {", name, n);

    for (size_t i = 0; i < n; i++) {
        offset += snprintf(buf + offset, buf_size - offset, "%d%s", arr[i], (i < n - 1) ? "," : "");
    }

    offset += snprintf(buf + offset, buf_size - offset, "};\n");

    return buf;
}

void read_kfile(int *klist, int *kcount){
    FILE *fp = fopen("seq.txt", "r");
    if (!fp) {
		fprintf(stderr, "ERROR: cannot open seq.txt!\n");
                printf( "ERROR: cannot open seq.txt!\n" );
		exit(EXIT_FAILURE);
	}

    size_t count = 0;

    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        int k;

        if (count >= MAX_SEQUENCES) {
            fprintf(stderr, "Reached maximum of %d sequences, stopping.\n", MAX_SEQUENCES);
            break;
        }

        if (sscanf(line, "%d", &k) == 1) {
            klist[count++] = k;
        } else {
            fprintf(stderr, "Failed to parse line: %s", line);
        }
    }

    fclose(fp);

	*kcount = count;
}

void cl_sieve( sclHard hardware, workStatus & st, searchData & sd ){

	progData pd = {};
	time_t boinc_last, ckpt_last, time_curr;
	cl_int err = 0;

	// setup kernel parameters
	setupSearch(st,sd);

	int h_klist[MAX_SEQUENCES];
	int kcount;
	read_kfile(h_klist, &kcount);
	sd.kcount = kcount;

	printf("Read in %d k:\n",kcount);
	for(int i=0; i<kcount; ++i){
		printf("%d ",h_klist[i]);
	}
	printf("\n");

	// device arrays
	pd.d_primecount = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, 12*sizeof(cl_uint), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_factor = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numresults*sizeof(factor), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: d_factor array.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

        pd.clearn = sclGetCLSoftware(clearn_cl,"clearn",hardware, NULL);
        pd.clearresult = sclGetCLSoftware(clearresult_cl,"clearresult",hardware, NULL);
        pd.addsmallprimes = sclGetCLSoftware(addsmallprimes_cl,"addsmallprimes",hardware, NULL);
	pd.init = sclGetCLSoftware(init_cl,"init",hardware, NULL);
	if(st.pmax < 0xFFFFFFFFFF000000){
	        pd.getsegprimes = sclGetCLSoftware(getsegprimes_cl,"getsegprimes",hardware, NULL);
	}
	else{
	       	pd.getsegprimes = sclGetCLSoftware(getsegprimes_cl,"getsegprimes",hardware, "-D CKOVERFLOW=1" );
	}

	uint32_t Q = 2048; 
	uint32_t L = st.nmax - st.nmin + 1;
	uint32_t m = (uint32_t) ceil((double)L / Q);

        size_t hsize = 1;
        while (hsize < 2*Q) hsize <<= 1;
	sd.hsize = (int32_t)hsize;

	// m reduced, 4 per iteration in giant kernel
	uint32_t mr = (uint32_t) ceil((double)m / 4);
	mr = (uint32_t) (ceil((double)mr / 32))*32;		// scale by warp size

	// for parity restricted P
	uint32_t QQ = Q<<1;
	uint32_t mm = (uint32_t) ceil((double)L / QQ);

	printf("-DHSIZE=%d -DMASK=%d -DQ=%u -DM=%u -DKCOUNT=%d -DNMIN=%u -DNMAX=%u -DQQ=%u -DMM=%u -DMR=%u\n",
		sd.hsize, sd.hsize-1, Q, m, sd.kcount, st.nmin, st.nmax, QQ, mm, mr);

	char cldef[256];
	snprintf(cldef, sizeof(cldef), "-DHSIZE=%d -DMASK=%d -DQ=%u -DM=%u -DKCOUNT=%d -DNMIN=%u -DNMAX=%u -DQQ=%u -DMM=%u -DMR=%u",
		sd.hsize, sd.hsize-1, Q, m, sd.kcount, st.nmin, st.nmax, QQ, mm, mr);

	// setup each k as part of a __constant kernel array
	char *const_str = generate_constant_array_string(h_klist, sd.kcount, "klist");
	
	int len_s = strlen(setup_cl);
	int len_g = strlen(giant_cl);
	int larg_len = (len_s > len_g) ? len_s : len_g;
	int total_len = larg_len + strlen(const_str);
	char src_str[total_len];

	snprintf(src_str, sizeof(src_str), "%s%s", const_str, setup_cl);
	pd.setup = sclGetCLSoftware(src_str,"setup",hardware, cldef);

	src_str[0] = '\0';
	snprintf(src_str, sizeof(src_str), "%s%s", const_str, giant_cl);
	pd.giantparity = sclGetCLSoftware(src_str,"giantparity",hardware, cldef);

//	pd.check = sclGetCLSoftware(check_cl,"combined_check",hardware, NULL);
	// kernel has __attribute__ ((reqd_work_group_size(256, 1, 1)))
	// it's still possible the CL complier picked a different size
	if(pd.getsegprimes.local_size[0] != 256){
		pd.getsegprimes.local_size[0] = 256;
		fprintf(stderr, "Set getsegprimes kernel local size to 256\n");
	}
/*	if(pd.check.local_size[0] != 256){
		pd.check.local_size[0] = 256;
		fprintf(stderr, "Set check kernel local size to 256\n");
	}
*/

	if( sd.test ){
		// clear result file
		FILE * temp_file = my_fopen(RESULTS_FILENAME,"w");
		if (temp_file == NULL){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
		fclose(temp_file);
	}
	else{
		// Resume from checkpoint if there is one
		if( read_state( st, sd ) ){
			if(boinc_is_standalone()){
				printf("Current p: %" PRIu64 "\n", st.p);
			}
			fprintf(stderr,"Resuming from checkpoint, current p: %" PRIu64 "\n", st.p);

			//trying to resume a finished workunit
			if( st.p == st.pmax ){
				if(boinc_is_standalone()){
					printf("Workunit complete.\n");
				}
				fprintf(stderr,"Workunit complete.\n");
				boinc_finish(EXIT_SUCCESS);
			}
		}
		// starting from beginning
		else{
			// clear result file
			FILE * temp_file = my_fopen(RESULTS_FILENAME,"w");
			if (temp_file == NULL){
				fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
				exit(EXIT_FAILURE);
			}
			fclose(temp_file);

			// setup boinc trickle up
			st.last_trickle = (uint64_t)time(NULL);
		}
	}

	// kernel used in profileGPU, setup arg
	sclSetKernelArg(pd.clearn, 0, sizeof(cl_mem), &pd.d_primecount);
	sclSetGlobalSize( pd.clearn, 64 );

	profileGPU(pd,st,sd,hardware);

	// number of gpu workgroups, used to size the sum array on gpu
//	sd.numgroups = (sd.psize / pd.check.local_size[0]) + 1;
	sd.numgroups = 2;

	// host arrays used for data transfer from gpu during checkpoints
	uint64_t * h_sum = (uint64_t *)malloc(sd.numgroups*sizeof(uint64_t));
	if( h_sum == NULL ){
		fprintf(stderr,"malloc error: h_sum\n");
		exit(EXIT_FAILURE);
	}
	uint32_t * h_primecount = (uint32_t *)malloc(12*sizeof(uint32_t));
	if( h_primecount == NULL ){
		fprintf(stderr,"malloc error: h_primecount\n");
		exit(EXIT_FAILURE);
	}

	uint64_t hs = (uint64_t)sd.psize * sd.hsize;
	sd.numhash = hs;
	hs *= sizeof(hash_entry);
//	printf("m:%d L:%d psize:%u\n",sd.m,sd.L,sd.psize);
	printf("gpu hash table has %" PRIu64 " elements\n",sd.numhash);
	printf("gpu hash table size is %d elements per thread\n",sd.hsize);
	printf("gpu hash table memory used is %" PRIu64 " megabytes\n", hs/1000000);
	printf("gpu max malloc is %" PRIu64 " megabytes\n", sd.maxmalloc/1000000);
	if(hs > sd.maxmalloc){
		fprintf(stderr, "ERROR: hash table too large!\n");
                printf( "ERROR: hash table too large!\n");
		exit(EXIT_FAILURE);
	}
	pd.d_htable = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, hs, NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: hash table\n");
                printf( "ERROR: clCreateBuffer failure: hash table\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_htable_even = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, hs, NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: hash table\n");
                printf( "ERROR: clCreateBuffer failure: hash table\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_htable_odd = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, hs, NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: hash table\n");
                printf( "ERROR: clCreateBuffer failure: hash table\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_hadj = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, sd.psize*sd.kcount*sizeof(cl_ulong), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: hash table\n");
                printf( "ERROR: clCreateBuffer failure: hash table\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_parity = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, sd.psize*sd.kcount*sizeof(cl_char), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: hash table\n");
                printf( "ERROR: clCreateBuffer failure: hash table\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_primes = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong4), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_primes_out = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong4), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_sum = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numgroups*sizeof(cl_ulong), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

	// setup global sizes
	sclSetGlobalSize( pd.getsegprimes, (sd.range/60)+1 );
	sclSetGlobalSize( pd.addsmallprimes, 64 );
	sclSetGlobalSize( pd.setup, sd.psize );
	sclSetGlobalSize( pd.init, sd.numhash );
	sclSetGlobalSize( pd.giantparity, sd.psize*mr );

//	sclSetGlobalSize( pd.check, sd.psize );
	sclSetGlobalSize( pd.clearresult, sd.numgroups );

	// set static kernel args
	sclSetKernelArg(pd.clearresult, 0, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.clearresult, 1, sizeof(cl_mem), &pd.d_sum);
	sclSetKernelArg(pd.clearresult, 2, sizeof(uint32_t), &sd.numgroups);

	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &pd.d_primecount);

	sclSetKernelArg(pd.addsmallprimes, 2, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.addsmallprimes, 3, sizeof(cl_mem), &pd.d_primecount);

	int ai = 0;
	sclSetKernelArg(pd.init, ai++, sizeof(cl_mem), &pd.d_htable);
	sclSetKernelArg(pd.init, ai++, sizeof(cl_mem), &pd.d_htable_even);
	sclSetKernelArg(pd.init, ai++, sizeof(cl_mem), &pd.d_htable_odd);
	sclSetKernelArg(pd.init, ai++, sizeof(uint64_t), &sd.numhash);
	ai = 0;
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_htable);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_htable_even);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_htable_odd);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_hadj);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_parity);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_sum);
	ai = 0;
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_factor);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_htable);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_htable_even);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_htable_odd);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_hadj);
	sclSetKernelArg(pd.giantparity, ai++, sizeof(cl_mem), &pd.d_parity);
	
/*
	sclSetKernelArg(pd.check, 0, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.check, 1, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.check, 2, sizeof(cl_mem), &pd.d_sum);
*/

	time(&boinc_last);
	time(&ckpt_last);
	time_t totals, totalf;
	if(boinc_is_standalone()){
		time(&totals);
	}

	float kernel_ms;
	int kernelq = 0;
	const int maxq = sd.compute ? 20 : 100;		// target kernel queue depth is 1 second
	cl_event launchEvent = NULL;
	const double irsize = 1.0 / (double)(st.pmax-st.pmin);

	sclEnqueueKernel(hardware, pd.clearresult);

	fprintf(stderr,"Starting Sieve...\n");
	if(boinc_is_standalone()){
		printf("Starting Sieve...\n");
	}

	// main search loop
	while(st.p < st.pmax){

		uint64_t stop = st.p + sd.range;
		if(stop > st.pmax || stop < st.p){
			// ck overflow
			stop = st.pmax;
		}

		// clear prime count
		sclEnqueueKernel(hardware, pd.clearn);
		
		time(&time_curr);
		if( ((int)time_curr - (int)boinc_last) > 1 ){
			// update BOINC fraction done every 2 sec
    			double fd = (double)(st.p-st.pmin)*irsize;
			boinc_fraction_done(fd);
			if(boinc_is_standalone()) printf("Sieve Progress: %.1f%%\n",fd*100.0);
			boinc_last = time_curr;
			if( ((int)time_curr - (int)ckpt_last) > 60 ){
				// 1 minute checkpoint
				if(kernelq > 0){
					waitOnEvent(hardware, launchEvent);
					kernelq = 0;
				}
				sleepCPU(hardware);
				boinc_begin_critical_section();
				getResults(pd, st, sd, hardware, h_primecount, h_sum);
				checkpoint(st, sd);
				boinc_end_critical_section();
				ckpt_last = time_curr;
				// clear result arrays
				sclEnqueueKernel(hardware, pd.clearresult);
			}
		}

		// add small primes that cannot be generated with getsegprimes kernel
		if(st.p < 114){
			uint64_t stop_sm = (stop > 114) ? 114 : stop;
			sclSetKernelArg(pd.addsmallprimes, 0, sizeof(uint64_t), &st.p);
			sclSetKernelArg(pd.addsmallprimes, 1, sizeof(uint64_t), &stop_sm);
			sclEnqueueKernel(hardware, pd.addsmallprimes);
			st.p = stop_sm;
		}

		// get a segment of primes (2-PRPs)
		int32_t wheelidx;
		uint64_t kernel_start = st.p;
		findWheelOffset(kernel_start, wheelidx);
		
		sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &kernel_start);
		sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &stop);
		sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &wheelidx);
		sclEnqueueKernel(hardware, pd.getsegprimes);

		// initialize hash tables
		if(kernelq == 0){
			launchEvent = sclEnqueueKernelEvent(hardware, pd.init);
		}
		else{
			sclEnqueueKernel(hardware, pd.init);
		}
		if(++kernelq == maxq){
			// limit cl queue depth and sleep cpu
			waitOnEvent(hardware, launchEvent);
			kernelq = 0;
		}

		// setup sieve prime data and baby steps
		sclEnqueueKernel(hardware, pd.setup);
//		kernel_ms = ProfilesclEnqueueKernel(hardware, pd.setup);
//		printf("setup kernel time %0.2fms\n",kernel_ms);

		// giant steps
		sclEnqueueKernel(hardware, pd.giantparity);
//		kernel_ms = ProfilesclEnqueueKernel(hardware, pd.giantparity);
//		printf("parity giant kernel time %0.2fms\n",kernel_ms);

		// checksum kernel
//		sclEnqueueKernel(hardware, pd.check);

		st.p = stop;

	}

	// final checkpoint
	if(kernelq > 0){
		waitOnEvent(hardware, launchEvent);
	}
	sleepCPU(hardware);

	boinc_begin_critical_section();
	st.p = st.pmax;
	boinc_fraction_done(1.0);
	if(boinc_is_standalone()) printf("Sieve Progress: %.1f%%\n",100.0);
	getResults(pd, st, sd, hardware, h_primecount, h_sum);
	checkpoint(st, sd);
	finalizeResults(st);
	boinc_end_critical_section();

	fprintf(stderr,"Sieve complete.\nfactors %" PRIu64 ", prime count %" PRIu64 "\n", st.factorcount, st.primecount);

	if(boinc_is_standalone()){
		time(&totalf);
		printf("Sieve finished in %d sec.\n", (int)totalf - (int)totals);
		printf("factors %" PRIu64 ", prime count %" PRIu64 ", checksum %016" PRIX64 "\n", st.factorcount, st.primecount, st.checksum);
	}

	free(h_sum);
	free(h_primecount);
	cleanup(pd, sd, st);
}


void reset_data(workStatus & st, searchData & sd){
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;
	sd.scount = 0;
	sd.powcount = 0;
	sd.prodcount = 0;
	st.factorial = false;
	st.primorial = false;
	st.compositorial = false;
}


void run_test( sclHard hardware, workStatus & st, searchData & sd ){

	fprintf(stderr,"self test not implemented yet!\n");
	printf("self test not implemented yet!\n");
	exit(EXIT_FAILURE);


	int goodtest = 0;

	printf("Beginning self test of 16 ranges.\n");

	time_t start, finish;
	time(&start);

	printf("Starting Factorial tests\n\n");
//	-p 100e6 -P 101e6 -n 1e6 -N 2e6 -!
	reset_data(st, sd);
	st.factorial = true;
	st.pmin = 100000000;
	st.pmax = 101000000;
	st.nmin = 1000000;
	st.nmax = 2000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 1071 && st.primecount == 54211 && st.checksum == 0x000004F844B5103C ){
		printf("test case 1 passed.\n\n");
		fprintf(stderr,"test case 1 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 1 failed.\n\n");
		fprintf(stderr,"test case 1 failed.\n");
	}

//	-p 1e12 -P 100001e7 -n 10000 -N 2e6 -!
	reset_data(st, sd);
	st.factorial = true;
	st.pmin = 1000000000000;
	st.pmax = 1000010000000;
	st.nmin = 10000;
	st.nmax = 2000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 3 && st.primecount == 361727 && st.checksum == 0x0505A1C238896511 ){
		printf("test case 2 passed.\n\n");
		fprintf(stderr,"test case 2 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 2 failed.\n\n");
		fprintf(stderr,"test case 2 failed.\n");
	}

//	-p 101 -P 100000 -n 101 -N 1e6 -!
	reset_data(st, sd);
	st.factorial = true;
	st.pmin = 101;
	st.pmax = 100000;
	st.nmin = 101;
	st.nmax = 1000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 42821 && st.primecount == 9571 && st.checksum == 0x0000000065DDB8A0 ){
		printf("test case 3 passed.\n\n");
		fprintf(stderr,"test case 3 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 3 failed.\n\n");
		fprintf(stderr,"test case 3 failed.\n");
	}

//	-p 1e12 -P 1000001e6 -n 10e7 -N 11e7 -!
	reset_data(st, sd);
	st.factorial = true;
	st.pmin = 1000000000000;
	st.pmax = 1000001000000;
	st.nmin = 100000000;
	st.nmax = 110000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 3 && st.primecount == 36249 && st.checksum == 0x00804FE7D7AA6C09 ){
		printf("test case 4 passed.\n\n");
		fprintf(stderr,"test case 4 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 4 failed.\n\n");
		fprintf(stderr,"test case 4 failed.\n");
	}

	printf("Starting Primorial tests\n\n");
//	-p 100e6 -P 101e6 -n 101 -N 25e6 -#
	reset_data(st, sd);
	st.primorial = true;
	st.pmin = 100000000;
	st.pmax = 101000000;
	st.nmin = 101;
	st.nmax = 25000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 1703 && st.primecount == 54211 && st.checksum == 0x0000027EFF497990 ){
		printf("test case 5 passed.\n\n");
		fprintf(stderr,"test case 5 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 5 failed.\n\n");
		fprintf(stderr,"test case 5 failed.\n");
	}

//	-p 101 -P 2e6 -n 101 -N 2e6 -#
	reset_data(st, sd);
	st.primorial = true;
	st.pmin = 101;
	st.pmax = 2000000;
	st.nmin = 101;
	st.nmax = 2000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 24503 && st.primecount == 148954 && st.checksum == 0x000000027BF5B8E0 ){
		printf("test case 6 passed.\n\n");
		fprintf(stderr,"test case 6 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 6 failed.\n\n");
		fprintf(stderr,"test case 6 failed.\n");
	}

//	-p 1e11 -P 100005e6 -n 9e6 -N 11e7 -#
	reset_data(st, sd);
	st.primorial = true;
	st.pmin = 100000000000;
	st.pmax = 100005000000;
	st.nmin = 9000000;
	st.nmax = 110000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 32 && st.primecount == 197222 && st.checksum == 0x0022FE7C09210B4B ){
		printf("test case 7 passed.\n\n");
		fprintf(stderr,"test case 7 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 7 failed.\n\n");
		fprintf(stderr,"test case 7 failed.\n");
	}

//	-n 600000 -N 30e6 -p 1730720716e6 -P 1730720720e6 -#
	reset_data(st, sd);
	st.primorial = true;
	st.pmin = 1730720716000000;
	st.pmax = 1730720720000000;
	st.nmin = 600000;
	st.nmax = 30000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 1 && st.primecount == 114208 && st.checksum == 0x5CDCB47F7E9532C2 ){
		printf("test case 8 passed.\n\n");
		fprintf(stderr,"test case 8 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 8 failed.\n\n");
		fprintf(stderr,"test case 8 failed.\n");
	}

	printf("Starting Compositorial tests\n\n");
//	-p 200e6 -P 20001e4 -n 101 -N 26e6 -c
	reset_data(st, sd);
	st.compositorial = true;
	st.pmin = 200000000;
	st.pmax = 200010000;
	st.nmin = 101;
	st.nmax = 26000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 127 && st.primecount == 529 && st.checksum == 0x0000001848D8AFBB ){
		printf("test case 9 passed.\n\n");
		fprintf(stderr,"test case 9 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 9 failed.\n\n");
		fprintf(stderr,"test case 9 failed.\n");
	}

//	-p 101 -P 1e5 -n 101 -N 1e6 -c
	reset_data(st, sd);
	st.compositorial = true;
	st.pmin = 101;
	st.pmax = 100000;
	st.nmin = 101;
	st.nmax = 1000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 34271 && st.primecount == 9571 && st.checksum == 0x000000006FF88EAE ){
		printf("test case 10 passed.\n\n");
		fprintf(stderr,"test case 10 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 10 failed.\n\n");
		fprintf(stderr,"test case 10 failed.\n");
	}

//	-p 2e11 -P 200005e6 -n 15e6 -N 2e7 -c
	reset_data(st, sd);
	st.compositorial = true;
	st.pmin = 200000000000;
	st.pmax = 200005000000;
	st.nmin = 15000000;
	st.nmax = 20000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 13 && st.primecount == 192386 && st.checksum == 0x0088B59C23CD3E2B ){
		printf("test case 11 passed.\n\n");
		fprintf(stderr,"test case 11 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 11 failed.\n\n");
		fprintf(stderr,"test case 11 failed.\n");
	}

//	-n 700000 -N 25e6 -p 1e12 -P 1000001e6 -c
	reset_data(st, sd);
	st.compositorial = true;
	st.pmin = 1000000000000;
	st.pmax = 1000001000000;
	st.nmin = 700000;
	st.nmax = 25000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 2 && st.primecount == 36249 && st.checksum == 0x0080997AF3BF42FE ){
		printf("test case 12 passed.\n\n");
		fprintf(stderr,"test case 12 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 12 failed.\n\n");
		fprintf(stderr,"test case 12 failed.\n");
	}

	printf("Starting Combined Factorial+Compositorial tests\n\n");
//	-p 1e11 -P 10001e7 -n 96000 -N 2e6 -! -c
	reset_data(st, sd);
	st.factorial = true;
	st.compositorial = true;
	st.pmin = 100000000000;
	st.pmax = 100010000000;
	st.nmin = 96000;
	st.nmax = 2000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 27 && st.primecount == 394403 && st.checksum == 0x00D214CC0EF0ECB4 ){
		printf("test case 13 passed.\n\n");
		fprintf(stderr,"test case 13 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 13 failed.\n\n");
		fprintf(stderr,"test case 13 failed.\n");
	}

//	-p 101 -P 11e4 -n 101 -N 1e6 -c -!
	reset_data(st, sd);
	st.factorial = true;
	st.compositorial = true;
	st.pmin = 101;
	st.pmax = 110000;
	st.nmin = 101;
	st.nmax = 1000000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 84077 && st.primecount == 10433 && st.checksum == 0x00000000EFB634E9 ){
		printf("test case 14 passed.\n\n");
		fprintf(stderr,"test case 14 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 14 failed.\n\n");
		fprintf(stderr,"test case 14 failed.\n");
	}

//	-p 101e8 -P 10101e6 -n 115e5 -N 125e5 -! -c
	reset_data(st, sd);
	st.factorial = true;
	st.compositorial = true;
	st.pmin = 10100000000;
	st.pmax = 10101000000;
	st.nmin = 11500000;
	st.nmax = 12500000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 19 && st.primecount == 43374 && st.checksum == 0x0002578EA9FD63C7 ){
		printf("test case 15 passed.\n\n");
		fprintf(stderr,"test case 15 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 15 failed.\n\n");
		fprintf(stderr,"test case 15 failed.\n");
	}

//	-p 2e12 -P 200002e7 -n 670000 -N 2460000 -! -c
	reset_data(st, sd);
	st.factorial = true;
	st.compositorial = true;
	st.pmin = 2000000000000;
	st.pmax = 2000020000000;
	st.nmin = 670000;
	st.nmax = 2460000;
	cl_sieve( hardware, st, sd );
	if( st.factorcount == 3 && st.primecount == 706162 && st.checksum == 0x1D63BBC574E8D50F ){
		printf("test case 16 passed.\n\n");
		fprintf(stderr,"test case 16 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 16 failed.\n\n");
		fprintf(stderr,"test case 16 failed.\n");
	}

//	done
	if(goodtest == 16){
		printf("All test cases completed successfully!\n");
		fprintf(stderr, "All test cases completed successfully!\n");
	}
	else{
		printf("Self test FAILED!\n");
		fprintf(stderr, "Self test FAILED!\n");
	}

	time(&finish);
	printf("Elapsed time: %d sec.\n", (int)finish - (int)start);

}


