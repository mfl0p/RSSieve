# RSSieve

RSSieve by Bryan Little

A BOINC enabled OpenCL Riesel Sierpinski Sieve for multiple sequences k*b^n+-1 using the BSGS algorithm.

With contributions by
* Geoffrey Reynolds
* Yves Gallot

## Requirements

* OpenCL v1.1
* 64 bit operating system

## How it works

1. Search parameters are given on the command line and input sr2sieve/sr5sieve style ABCD file.
2. A small group of sieve primes are generated on the GPU.
3. The group of primes are tested for factors using the BSGS algorithm.
4. Repeat #2-3 until checkpoint.  Gather and verify factors from GPU.
5. Report any factors in psp_sr2sieve.out

## Running the program
```
command line options
* -p #			Starting prime factor p
* -P #			End prime factor P
* 			P range is 2^32 <= -p < -P < 2^64, [-p, -P) exclusive
* -i inputfile		Use specified sr2sieve ABCD input file with a maximum of 100 sequences
* -h			Print help

Note that you will need to use sr2sieve/sr5sieve for P below 2^32.

Program gets the OpenCL GPU device index from BOINC.  To run stand-alone, the program will
default to GPU 0 unless an init_data.xml is in the same directory with the format:

<app_init_data>
<gpu_type>NVIDIA</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>

or

<app_init_data>
<gpu_type>ATI</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>
```

## Related Links
* [Yves Gallot on GitHub](https://github.com/galloty)

