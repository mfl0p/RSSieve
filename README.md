# RSSieve

RSSieve by Bryan Little

A BOINC enabled OpenCL standalone Riesel Sierpinski Sieve for multiple sequences k*b^n+-1 using the BSGS algorithm.

With contributions by
* Geoffrey Reynolds
* Yves Gallot

	THIS IS AN INCOMPLETE VERSION!

## Requirements

* OpenCL v1.1
* 64 bit operating system

## How it works

1. Search parameters are given on the command line and seq.txt file.
2. A small group of sieve primes are generated on the GPU.
3. The group of primes are tested for factors in the N range specified.
4. Repeat #2-3 until checkpoint.  Gather factors and checksum data from GPU.
5. Report any factors in factors.txt, along with a checksum at the end.
6. Checksum can be used to compare results in a BOINC quorum.

## Running the program
```
A list of up to 200 k, one per line, are specified in seq.txt where negative k is for Riesel and positive k is for Sierpinski.
command line options
* -b #  base
* -n #	Start n
* -N #	End N
* 		N range is 101 <= -n < -N < 2^31, [-n, -N) exclusive
* -p #	Starting prime factor p
* -P #	End prime factor P
* 		P range is 3 <= -p < -P < 2^64, [-p, -P) exclusive
* -s 	Perform self test to verify proper operation of the program with the current GPU.
* -h	Print help

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

