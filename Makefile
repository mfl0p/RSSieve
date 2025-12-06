CC = g++
LD = $(CC)

.SUFFIXES:
.SUFFIXES: .o .c .h .cl .cpp

VERSION_MAJOR := 1
VERSION_MINOR := 4
date := $(shell powershell.exe get-date -format FileDate)

APP = PFCSieve-win64-v$(VERSION_MAJOR).$(VERSION_MINOR)-$(date).exe

SRC = main.cpp cl_sieve.cpp cl_sieve.h simpleCL.c simpleCL.h kernels/check.cl kernels/clearn.cl kernels/clearresult.cl kernels/getsegprimes.cl kernels/addsmallprimes.cl kernels/iterate.cl kernels/setup.cl kernels/verifyslow.cl kernels/verify.cl kernels/verifyresult.cl putil.c putil.h verifyprime.cpp verifyprime.h
KERNEL_HEADERS = kernels/check.h kernels/clearn.h kernels/clearresult.h kernels/iterate.h kernels/setup.h kernels/getsegprimes.h kernels/addsmallprimes.h kernels/verifyslow.h kernels/verify.h kernels/verifyresult.h
OBJ = main.o cl_sieve.o simpleCL.o putil.o verifyprime.o

LIBS = OpenCL.dll libprimesievewin.a

BOINC_DIR = C:/mingwbuilds/boinc
BOINC_INC = -I$(BOINC_DIR)/lib -I$(BOINC_DIR)/api -I$(BOINC_DIR) -I$(BOINC_DIR)/win_build
BOINC_LIB = -L$(BOINC_DIR)/lib -L$(BOINC_DIR)/api -L$(BOINC_DIR) -lboinc_opencl -lboinc_api -lboinc

CFLAGS  = -I . -I kernels -O3 -m64 -Wall -DVERSION_MAJOR=\"$(VERSION_MAJOR)\" -DVERSION_MINOR=\"$(VERSION_MINOR)\"
LDFLAGS = $(CFLAGS) -lstdc++ -static -fopenmp

all : clean $(APP)

$(APP) : $(OBJ)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(BOINC_LIB) -o $@

main.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -fopenmp -c -o $@ main.cpp

cl_sieve.o : $(SRC) $(KERNEL_HEADERS)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -fopenmp -c -o $@ cl_sieve.cpp

simpleCL.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ simpleCL.c

putil.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ putil.c

verifyprime.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ verifyprime.cpp

.cl.h:
	perl cltoh.pl $< > $@

clean :
	del *.o
	del kernels\*.h
	del $(APP)

