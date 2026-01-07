
// cl_sieve.h

typedef struct {
	cl_ulong hash;
	cl_short idx;
} hash_entry;

typedef struct {
	cl_ulong p;
	cl_int n;
	cl_int k;
} factor;

typedef struct {
	cl_ulong hadj;
	cl_int parity;
	cl_int kidx;
} kdata;

typedef struct {
	uint64_t pmin, pmax, p, checksum, primecount, factorcount, last_trickle, state_sum;
	uint32_t nmin, nmax, base;
} workStatus;

typedef struct {
	uint64_t maxmalloc, numhash;
	uint32_t computeunits, nstep, sstep, powcount, prodcount, scount, numresults, range, psize, numgroups, nlimit;
	int32_t kcount, L, m, hsize;
	bool test, compute, write_state_a_next;
} searchData;

typedef struct {
	cl_mem d_factor;
	cl_mem d_sum;
	cl_mem d_primes;
	cl_mem d_primecount;
	cl_mem d_htable, d_htable_even, d_htable_odd;
	cl_mem d_k;
	sclSoft setup, giant, clearn, clearresult, getsegprimes, addsmallprimes, init, giantparity;
} progData;

void cl_sieve( sclHard hardware, workStatus & st, searchData & sd );

void run_test( sclHard hardware, workStatus & st, searchData & sd );
