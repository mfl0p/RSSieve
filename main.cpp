/*
	THIS IS AN INCOMPLETE VERSION

	RSSieve
	Bryan Little, Feb 2026

	with contributions by Geoffrey Reynolds and Yves Gallot

	Required minimum OpenCL version is 1.1
	CL_TARGET_OPENCL_VERSION 110 in simpleCL.h

	Search limits:  P up to 2^64 and N up to 2^31


*/

#include <unistd.h>
#include <getopt.h>
#include <cinttypes>

#include "boinc_api.h"
#include "boinc_opencl.h"
#include "simpleCL.h"
#include "cl_sieve.h"

void help()
{
	printf("Welcome to RSSieve, an OpenCL Riesel Sierpinski Sieve for multiple sequences k*b^n+-1 using the BSGS algorithm.\n");
//	printf("A list of up to 100 k, one per line, are specified in seq.txt where negative k is for Riesel and positive k is for Sierpinski.\n");
	printf("Program usage:\n");
	printf("-b #	base, 2 <= b < 2^31\n");
	printf("-n #	Start n\n");
	printf("-N #	End N\n");
	printf("		N range is 101 <= -n < -N < 2^31, [-n, -N) exclusive\n");
	printf("-p #	Starting prime factor p\n");
	printf("-P #	End prime factor P\n");
	printf("		P range is 3 <= -p < -P < 2^64, [-p, -P) exclusive\n");
	printf("-i inputfile.abcd\n");
	printf("-s 	Perform self test to verify proper operation of the program with the current GPU.\n");
	printf("-h	Print this help\n");
        boinc_finish(EXIT_FAILURE);
}


/* Parse uint64 with optional 'e' exponent, returns 0=ok, -1=invalid, -2=out of range */
int parse_uint64(uint64_t *result, const char *str,
                 uint64_t lo, uint64_t hi)
{
    if (!str || !*str) return -1;

    char *tail;
    errno = 0;
    uint64_t num = strtoull(str, &tail, 10);
    if (errno != 0) return -2;

    if (*tail == 'e' || *tail == 'E') {
        tail++;
        errno = 0;
        uint64_t exp = strtoul(tail, &tail, 10);
        if (errno != 0) return -2;
        if (*tail != '\0') return -1;

        while (exp-- > 0) {
            if (num > hi / 10) return -2;
            num *= 10;
        }
    }
    else if (*tail != '\0') {
        return -1; // invalid character
    }

    if (num < lo || num > hi) return -2;

    *result = num;
    return 0;
}

/* parse_uint simply wraps parse_uint64 for 32-bit values */
int parse_uint(uint32_t *result, const char *str,
               uint32_t lo, uint32_t hi)
{
    uint64_t tmp;
    int status = parse_uint64(&tmp, str, lo, hi);
    if (status == 0) *result = (uint32_t)tmp;
    return status;
}

/* Full command line parser from string */
void parse_cmdline_string(const char *cmdline, workStatus *st, searchData *sd)
{
    if (!cmdline || !st || !sd) return;

    char buf[8192];
    strncpy(buf, cmdline, sizeof(buf)-1);
    buf[sizeof(buf)-1] = '\0';

    char *token = strtok(buf, " \t");
    while (token) {
        if (strcmp(token, "-b") == 0) {
            token = strtok(NULL, " \t");
            if (token && parse_uint(&st->base, token, 2, 0x7FFFFFFF) != 0) {
                fprintf(stderr, "Invalid value for -b: %s\n", token);
            }
        }
        else if (strcmp(token, "-p") == 0) {
            token = strtok(NULL, " \t");
            if (token && parse_uint64(&st->pmin, token, 1, 0xFFFFFFFFFFFFFFFFULL) != 0) {
                fprintf(stderr, "Invalid value for -p: %s\n", token);
            }
        }
        else if (strcmp(token, "-P") == 0) {
            token = strtok(NULL, " \t");
            if (token && parse_uint64(&st->pmax, token, 1, 0xFFFFFFFFFFFFFFFFULL) != 0) {
                fprintf(stderr, "Invalid value for -P: %s\n", token);
            }
        }
        else if (strcmp(token, "-n") == 0) {
            token = strtok(NULL, " \t");
            if (token && parse_uint(&st->nmin, token, 1, 0x7FFFFFFF) != 0) {
                fprintf(stderr, "Invalid value for -n: %s\n", token);
            }
        }
        else if (strcmp(token, "-N") == 0) {
            token = strtok(NULL, " \t");
            if (token && parse_uint(&st->nmax, token, 1, 0x7FFFFFFF) != 0) {
                fprintf(stderr, "Invalid value for -N: %s\n", token);
            }
        }
	else if (strcmp(token, "-i") == 0) {
	    token = strtok(NULL, " \t");
	    if (token) {
		sd->input_file = strdup(token);  // allocate a copy
		if (!sd->input_file) {
		    fprintf(stderr, "Failed to allocate memory for input file\n");
		}
	    }
	}
        else if (strcmp(token, "-s") == 0) {
            sd->test = true;
        }
        else if (strcmp(token, "-h") == 0) {
            fprintf(stderr, "Help requested\n");
            // call help() if you have it
        }
        // unknown flag or extra token → ignore silently
        token = strtok(NULL, " \t");
    }
}

/* Join argc/argv into a single string with spaces */
char* join_argv(int argc, char *argv[])
{
    if (argc <= 0) return NULL;

    // Estimate required buffer size
    size_t total = 1;
    for (int i = 1; i < argc; i++)
        total += strlen(argv[i]) + 1; // +1 for space or null terminator

    char *cmdline = (char *)malloc(total);
    if (!cmdline) return NULL;

    cmdline[0] = '\0';

    for (int i = 1; i < argc; i++) {
        strcat(cmdline, argv[i]);
        if (i < argc - 1) strcat(cmdline, " ");
    }

    return cmdline; // must free later
}


int main(int argc, char *argv[])
{ 
	sclHard hardware = {};
	searchData sd = {};
	sd.numresults = 10000000;
	sd.write_state_a_next = true;
	workStatus st = {};

        // Initialize BOINC
        BOINC_OPTIONS options;
        boinc_options_defaults(options);
        options.normal_thread_priority = true;
        boinc_init_options(&options);

	fprintf(stderr, "\nNOTE THIS IS AN INCOMPLETE VERSION!\n");
	printf("\nNOTE THIS IS AN INCOMPLETE VERSION!\n");

	fprintf(stderr, "\nRSSieve v%s.%s by Bryan Little\nwith contributions by Geoffrey Reynolds and Yves Gallot\n",VERSION_MAJOR,VERSION_MINOR);
	fprintf(stderr, "Compiled " __DATE__ " with GCC " __VERSION__ "\n");
	if(boinc_is_standalone()){
		printf("\nRSSieve v%s.%s by Bryan Little\nwith contributions by Geoffrey Reynolds and Yves Gallot\n",VERSION_MAJOR,VERSION_MINOR);
		printf("Compiled " __DATE__ " with GCC " __VERSION__ "\n");
	}

        // Print out cmd line for diagnostics
        fprintf(stderr, "Command line: ");
        for (int i = 0; i < argc; i++)
        	fprintf(stderr, "%s ", argv[i]);
        fprintf(stderr, "\n");

	char *cmdline = join_argv(argc, argv);
	if (!cmdline) {
		fprintf(stderr, "Failed to build command line string\n");
		return 1;
	}

	// hack to work around invalid args when running under app_info
	printf("%s\n",cmdline);
	parse_cmdline_string(cmdline, &st, &sd);
	free(cmdline);


	cl_platform_id platform = 0;
	cl_device_id device = 0;
	cl_context ctx;
	cl_command_queue queue;
	cl_int err = 0;

	int retval = 0;
	retval = boinc_get_opencl_ids(argc, argv, 0, &device, &platform);
	if (retval) {
		if(boinc_is_standalone()){
			printf("init_data.xml not found, using device 0.\n");

			err = clGetPlatformIDs(1, &platform, NULL);
			if (err != CL_SUCCESS) {
				printf( "clGetPlatformIDs() failed with %d\n", err );
				fprintf(stderr, "Error: clGetPlatformIDs() failed with %d\n", err );
				exit(EXIT_FAILURE);
			}
			err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
			if (err != CL_SUCCESS) {
				printf( "clGetDeviceIDs() failed with %d\n", err );
				fprintf(stderr, "Error: clGetDeviceIDs() failed with %d\n", err );
				exit(EXIT_FAILURE);
			}
		}
		else{
			fprintf(stderr, "Error: boinc_get_opencl_ids() failed with error %d\n", retval );
			exit(EXIT_FAILURE);
		}
	}

	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

	ctx = clCreateContext(cps, 1, &device, NULL, NULL, &err);
	if (err != CL_SUCCESS) {
		fprintf(stderr, "Error: clCreateContext() returned %d\n", err);
        	exit(EXIT_FAILURE); 
   	}

	// OpenCL v2.0
	//cl_queue_properties qp[] = { CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0 };
	//queue = clCreateCommandQueueWithProperties(ctx, device, qp, &err);

	queue = clCreateCommandQueue(ctx, device, CL_QUEUE_PROFILING_ENABLE, &err);	
	if(err != CL_SUCCESS) { 
		fprintf(stderr, "Error: Creating Command Queue. (clCreateCommandQueueWithProperties) returned %d\n", err );
		exit(EXIT_FAILURE);
    	}

	hardware.platform = platform;
	hardware.device = device;
	hardware.queue = queue;
	hardware.context = ctx;

 	char device_name[1024];
 	char device_vend[1024];
 	char device_driver[1024];
	cl_uint CUs;
	cl_ulong LMS = 0;

	err = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), &device_name, NULL);
	if (err != CL_SUCCESS) {
		printf( "clGetDeviceInfo failed with %d\n", err );
		exit(EXIT_FAILURE);
	}
	err = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(device_vend), &device_vend, NULL);
	if (err != CL_SUCCESS) {
		printf( "clGetDeviceInfo failed with %d\n", err );
		exit(EXIT_FAILURE);
	}
	err = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(device_driver), &device_driver, NULL);
	if (err != CL_SUCCESS) {
		printf( "clGetDeviceInfo failed with %d\n", err );
		exit(EXIT_FAILURE);
	}
	err = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &CUs, NULL);
	if (err != CL_SUCCESS) {
		printf( "clGetDeviceInfo failed with %d\n", err );
		exit(EXIT_FAILURE);
	}
	err = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &LMS, NULL);
	if (err != CL_SUCCESS) {
		printf( "clGetDeviceInfo failed with %d\n", err );
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "GPU Info:\n  Name: \t\t%s\n  Vendor: \t\t%s\n  Driver: \t\t%s\n  Compute Units: \t%u\n  Local Mem Size: \t%u bytes\n",
		device_name, device_vend, device_driver, CUs, (uint32_t)LMS);
	if(boinc_is_standalone()){
		printf("GPU Info:\n  Name: \t\t%s\n  Vendor: \t\t%s\n  Driver: \t\t%s\n  Compute Units: \t%u\n  Local Mem Size: \t%u bytes\n",
			device_name, device_vend, device_driver, CUs, (uint32_t)LMS);
	}

	sd.computeunits = (uint32_t)CUs;
	sd.lmemsize = (uint32_t)LMS;

	char intel_s[] = "Intel";
	char arc_s[] = "Arc";
	char nvidia_s[] = "NVIDIA";	

	if(strstr((char*)device_vend, (char*)nvidia_s) != NULL){
		sd.nvidia = true;
	}
	// TODO: normalize intel compute units
	else if( strstr((char*)device_vend, (char*)intel_s) != NULL ){
		if( strstr((char*)device_name, (char*)arc_s) != NULL ){
			sd.computeunits /= 10;
		}
		else{
			sd.computeunits /= 20;
	                fprintf(stderr,"Detected Intel integrated graphics\n");	
		}
	}

	if(!sd.computeunits) sd.computeunits++;
	
	if(sd.test){
		run_test(hardware, st, sd);
	}
	else{
		cl_sieve(hardware, st, sd);
	}

        sclReleaseClHard(hardware);

	boinc_finish(EXIT_SUCCESS);

	return 0; 
} 

