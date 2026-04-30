# RSSieve

RSSieve by Bryan Little

A BOINC enabled OpenCL Riesel Sierpinski Sieve for multiple sequences k*b^n+-1 using the Baby-Step Giant-Step algorithm.

With contributions by
* Geoffrey Reynolds
* Yves Gallot

## Requirements

* OpenCL v1.1
* 64 bit operating system

## How it works

1. Search parameters are given on the command line and input sr2sieve/sr5sieve style ABCD file.
2. A small group of sieve primes are generated on the GPU.
3. The group of primes are tested for factors using the Baby-Step Giant-Step algorithm.
4. Repeat #2-3 until checkpoint.  Gather and verify factors from GPU.
5. Report any factors in factors file.

A note about comparing factor output to other programs at small P and high factor density.
Due to the parallel nature of the program and the use of atomic operations, only the first factor
found for a sequence will be reported in the factor file.  Multiple tests of the same range will
report different factors each time!  A better way to compare the factor files with another program
would be to compare the sequence k*b^n+-1 of the factor, not the prime factor itself.
Basically, you are comparing what sequence is filtered out, not the primes.  If both program's
factor files are filtering the same sequences, then output is effectively the same.
At higher P ranges with low factor density, the factor files will match.

## Running the program
```
command line options
* -p #			Starting prime factor p
* -P #			End prime factor P
* 			P range is 3 <= -p < -P < 2^64, [-p, -P) exclusive
* -i inputfile		Use specified sr2sieve ABCD input file with a maximum of 100 sequences
* -f factorfile		Override default factor file name (psp_sr2sieve.out) with specified file name.
* -h			Print help

Note that you will need to use sr2sieve/sr5sieve to generate the initial ABCD input file.

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

