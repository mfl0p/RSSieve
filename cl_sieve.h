
#define MAX_SEQUENCES 100

typedef struct {
    long long K;
    char sign;
    long long N0;
    long long lastN;
    size_t nbits;        // bitmap bits
    uint8_t *bitmap;     // 1 = can be used
} Sequence;

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
	cl_ulong hadj;
	cl_int kidx;
} kparity;

typedef struct {
	uint64_t pmin, pmax, p, checksum, primecount, factorcount, last_trickle, state_sum;
	uint32_t nmin, nmax, base;
} workStatus;

typedef struct {
	uint32_t computeunits, numresults, range, psize, numgroups, Q, m, QQ, mm, hsize;
	int32_t kcount;
	bool test, compute, write_state_a_next, nvidia;
	Sequence sequences[MAX_SEQUENCES];
	int klist[MAX_SEQUENCES];
	const char *input_file;
} searchData;

typedef struct {
	cl_mem d_factor;
	cl_mem d_sum;
	cl_mem d_primes, d_primes_full, d_primes_even, d_primes_odd, d_htable, d_hidx;
	cl_mem d_primecount;
	cl_mem d_k, d_k_full, d_k_even, d_k_odd;
	cl_mem d_kcount_full, d_kcount_even, d_kcount_odd;
	sclSoft setup, sort, clearn, clearresult, getsegprimes, addsmallprimes, giantparity, giantfull;
} progData;

void cl_sieve( sclHard hardware, workStatus & st, searchData & sd );

void run_test( sclHard hardware, workStatus & st, searchData & sd );

void read_input(workStatus & st, searchData & sd);

int factor_can_be_used(Sequence *sequences, size_t count, long long K, char sign, long long n);

void mark_factor_used(Sequence *sequences, size_t count, long long K, char sign, long long n);


