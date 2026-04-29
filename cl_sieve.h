
#define MAX_SEQUENCES 100

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
	uint64_t pmin, pmax, p, checksum, primecount, factorcount, last_trickle, state_sum, dupcount;
	uint32_t nmin, nmax, base;
	int32_t kcount;
	int klist[MAX_SEQUENCES];
} workStatus;

typedef struct {
	uint32_t computeunits, numresults, range, psize, Q, m, QQ, mm, hsize, lmemsize, primes_per_bsgs, cache_k;
	bool test, write_state_a_next, nvidia;
	size_t bitmap_bits_per_sequence;
	size_t bitmap_words_per_sequence;
	size_t bitmap_total_words;
	const char *input_file;
	uint32_t *bitmap;
} searchData;

typedef struct {
	cl_mem d_factor;
	cl_mem d_sum;
	cl_mem d_primes, d_primes_full, d_primes_even, d_primes_odd, d_htable, d_hidx;
	cl_mem d_primecount;
	cl_mem d_k, d_k_full, d_k_even, d_k_odd;
	cl_mem d_kcount_full, d_kcount_even, d_kcount_odd;
	cl_mem d_bsgs_count;
	cl_mem d_bitmap;
	sclSoft addsmallprimes, setup, sort, clearn, clearresult, getsegprimes, giantparity, giantfull;
} progData;

void cl_sieve( sclHard hardware, workStatus & st, searchData & sd );

void run_test( sclHard hardware, workStatus & st, searchData & sd );

void read_input(workStatus & st, searchData & sd);

void free_sequences(searchData &sd, workStatus &st);

