//
// sr5sieve-style ABCD reader using one large uint32_t bitmap array.
//
// Parsing rules:
//  - Ignore blank lines and comment lines (# as first non-whitespace char)
//  - If line matches "%u" => delta line (must be after a header):
//        running_n += delta; record running_n
//  - Else must be an ABCD header:
//        ABCD k*b^$a(+/-1) [n] // Sieved to p
//    (" // Sieved to p" optional)
//  - C must be +1 or -1
//  - B must be the same for all sequences
//  - Any other non-empty/non-comment line is a hard error
//
// Bitmap layout:
//  - sd.bitmap is one large uint32_t array.
//  - Every sequence has the same logical N range: st.nmin..st.nmax.
//  - st.nmin is the smallest N read from all sequences, reduced by 1 if odd.
//  - Sequence i starts at: i * sd.bitmap_words_per_sequence.
//  - Bit for N is: (N - st.nmin).
//  - Only N values read from the file are set. All other positions stay zero.
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <inttypes.h>

#include "boinc_api.h"
#include "simpleCL.h"
#include "cl_sieve.h"

// -----------------------------------------------------------------------------
// One-large-array uint32_t bitmap helpers + factor usage helpers
// -----------------------------------------------------------------------------

static uint32_t* init_bitmap_words_zero(size_t nwords) {
    if (nwords == 0) return nullptr;

    if (nwords > ((size_t)-1) / sizeof(uint32_t)) {
        std::fprintf(stderr, "bitmap allocation overflow\n");
        std::exit(EXIT_FAILURE);
    }

    uint32_t *bm = (uint32_t*)std::calloc(nwords, sizeof(uint32_t));
    if (!bm) { std::perror("calloc"); std::exit(EXIT_FAILURE); }
    return bm;
}

static inline size_t bitmap_seq_base_word(const searchData &sd, uint32_t seqno) {
    return (size_t)seqno * sd.bitmap_words_per_sequence;
}

static inline void set_can_use(searchData &sd, const workStatus &st, uint32_t seqno, uint32_t n) {
    if (!sd.bitmap) return;
    if (n < st.nmin || n > st.nmax) return;

    const size_t bit_offset = (size_t)((uint64_t)n - (uint64_t)st.nmin);
    const size_t word_index = bitmap_seq_base_word(sd, seqno) + (bit_offset >> 5);
    const uint32_t mask = UINT32_C(1) << (bit_offset & 31u);

    sd.bitmap[word_index] |= mask;
}

static inline void clear_can_use(searchData &sd, const workStatus &st, uint32_t seqno, uint32_t n) {
    if (!sd.bitmap) return;
    if (n < st.nmin || n > st.nmax) return;

    const size_t bit_offset = (size_t)((uint64_t)n - (uint64_t)st.nmin);
    const size_t word_index = bitmap_seq_base_word(sd, seqno) + (bit_offset >> 5);
    const uint32_t mask = UINT32_C(1) << (bit_offset & 31u);

    sd.bitmap[word_index] &= ~mask;
}

static inline int can_use_N(const searchData &sd, const workStatus &st, uint32_t seqno, uint32_t n) {
    if (!sd.bitmap) return 0;
    if (n < st.nmin || n > st.nmax) return 0;

    const size_t bit_offset = (size_t)((uint64_t)n - (uint64_t)st.nmin);
    const size_t word_index = bitmap_seq_base_word(sd, seqno) + (bit_offset >> 5);
    const uint32_t mask = UINT32_C(1) << (bit_offset & 31u);

    return (sd.bitmap[word_index] & mask) != 0;
}

static int find_sequence(const workStatus &st, uint32_t K, char sign) {
    const int signed_k = (sign == '+') ? (int)K : -(int)K;

    for (int i = 0; i < st.kcount; i++) {
        if (st.klist[i] == signed_k) return i;
    }
    return -1;
}

int factor_can_be_used(const searchData &sd, const workStatus &st, uint32_t K, char sign, uint32_t n) {
    int idx = find_sequence(st, K, sign);
    if (idx < 0) return 0;
    return can_use_N(sd, st, (uint32_t)idx, n);
}

void mark_factor_used(searchData &sd, const workStatus &st, uint32_t K, char sign, uint32_t n) {
    int idx = find_sequence(st, K, sign);
    if (idx >= 0) clear_can_use(sd, st, (uint32_t)idx, n);
}

// Keep the old name if other code already calls free_sequences().
// It now frees the single large bitmap instead of per-Sequence allocations.
void free_sequences(searchData &sd, workStatus &st)
{
    if (sd.bitmap) {
        std::free(sd.bitmap);
        sd.bitmap = nullptr;
    }

    sd.bitmap_bits_per_sequence = 0;
    sd.bitmap_words_per_sequence = 0;
    sd.bitmap_total_words = 0;

    st.kcount = 0;
}

// -----------------------------------------------------------------------------
// sr5sieve-like error reporting
// -----------------------------------------------------------------------------

static void line_error(const char *file_name,
                       uint64_t line_no,
                       const char *msg,
                       const char *line,
                       int current_seq /* -1 if none */)
{
    if (current_seq >= 0) {
        std::fprintf(stderr, "%s:%" PRIu64 ": (seq %d) %s\n",
                     file_name ? file_name : "<input>", line_no, current_seq, msg);
    } else {
        std::fprintf(stderr, "%s:%" PRIu64 ": %s\n",
                     file_name ? file_name : "<input>", line_no, msg);
    }

    if (line && *line) {
        std::fprintf(stderr, "  >> %s", line);
        size_t len = std::strlen(line);
        if (len == 0 || (line[len - 1] != '\n' && line[len - 1] != '\r')) std::fputc('\n', stderr);
    }
    std::exit(EXIT_FAILURE);
}

static void file_error_open(const char *file_name) {
    std::fprintf(stderr, "Error: unable to open input file %s\n", file_name ? file_name : "<null>");
    std::perror("fopen");
    std::exit(EXIT_FAILURE);
}

// -----------------------------------------------------------------------------
// Minimal dynamic list for per-sequence N values
// -----------------------------------------------------------------------------

typedef struct {
    uint32_t *v;
    size_t    n;
    size_t    cap;
} NVec;

static void nvec_push(NVec *a, uint32_t x) {
    if (a->n == a->cap) {
        size_t newcap = (a->cap ? a->cap * 2 : 256);

        if (newcap > ((size_t)-1) / sizeof(uint32_t)) {
            std::fprintf(stderr, "N list allocation overflow\n");
            std::exit(EXIT_FAILURE);
        }

        uint32_t *nv = (uint32_t*)std::realloc(a->v, newcap * sizeof(uint32_t));
        if (!nv) { std::perror("realloc"); std::exit(EXIT_FAILURE); }
        a->v = nv;
        a->cap = newcap;
    }
    a->v[a->n++] = x;
}

// -----------------------------------------------------------------------------
// sr5sieve-style ABCD header parse
// -----------------------------------------------------------------------------

static int parse_abcd_header_sr5(
    const char *line,
    uint32_t *k_out,
    uint32_t *b_out,
    int32_t  *c_out,
    uint32_t *n_out
) {
    // With tail
    {
        uint64_t p_dummy = 0;
        uint32_t k = 0, b = 0, n = 0;
        int32_t  c = 0;

        int got = std::sscanf(line,
            "ABCD %" SCNu32 "*%" SCNu32 "^$a%" SCNd32 " [%" SCNu32 "] // Sieved to %" SCNu64,
            &k, &b, &c, &n, &p_dummy);

        if (got == 5) {
            *k_out = k; *b_out = b; *c_out = c; *n_out = n;
            return 1;
        }
    }

    // Without tail
    {
        uint32_t k = 0, b = 0, n = 0;
        int32_t  c = 0;

        int got = std::sscanf(line,
            "ABCD %" SCNu32 "*%" SCNu32 "^$a%" SCNd32 " [%" SCNu32 "]",
            &k, &b, &c, &n);

        if (got == 4) {
            *k_out = k; *b_out = b; *c_out = c; *n_out = n;
            return 1;
        }
    }

    return 0;
}

// -----------------------------------------------------------------------------
// Safe add for uint32_t with overflow detection
// -----------------------------------------------------------------------------

static inline uint32_t add_u32_checked(uint32_t a, uint32_t b,
                                      const char *file_name, uint64_t line_no,
                                      const char *line, int current_seq,
                                      const char *what)
{
    uint32_t r = a + b;
    if (r < a) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "overflow while %s (uint32_t wrap)", what);
        line_error(file_name, line_no, buf, line, current_seq);
    }
    return r;
}

// -----------------------------------------------------------------------------
// Main reader
// -----------------------------------------------------------------------------

void read_input(workStatus &st, searchData &sd)
{
    if (!sd.input_file || !sd.input_file[0]) {
        std::fprintf(stderr, "Error: input_file is null/empty\n");
        std::exit(EXIT_FAILURE);
    }

    // If the caller reuses searchData, release the previous bitmap first.
    // This assumes searchData is zero-initialized before first use.
    if (sd.bitmap) {
        std::free(sd.bitmap);
        sd.bitmap = nullptr;
    }
    sd.bitmap_bits_per_sequence = 0;
    sd.bitmap_words_per_sequence = 0;
    sd.bitmap_total_words = 0;
    st.kcount = 0;

    char resolved_name[512];
    boinc_resolve_filename(sd.input_file, resolved_name, sizeof(resolved_name));
    FILE *fp = boinc_fopen(resolved_name, "r");

    if (!fp) file_error_open(sd.input_file);

    NVec nlists[MAX_SEQUENCES];
    std::memset(nlists, 0, sizeof(nlists));

    uint32_t firstN[MAX_SEQUENCES];
    uint32_t lastN[MAX_SEQUENCES];
    std::memset(firstN, 0, sizeof(firstN));
    std::memset(lastN,  0, sizeof(lastN));

    char line[2048];
    uint64_t line_no = 0;

    int current_seq = -1;      // no active sequence until first ABCD header
    uint32_t running_n = 0;    // current n for delta accumulation

    int32_t global_B = -1;
    uint32_t NMIN = UINT32_MAX;
    uint32_t NMAX = 0;

    while (std::fgets(line, sizeof(line), fp)) {
        line_no++;

        // Find first non-whitespace char (for blank/comment detection)
        const char *s = line;
        while (*s && std::isspace((unsigned char)*s)) s++;

        if (*s == '\0' || *s == '#')
            continue;

        // Delta line?
        uint32_t delta = 0;
        if (std::sscanf(line, "%" SCNu32, &delta) == 1) {
            if (current_seq < 0) {
                std::fclose(fp);
                line_error(sd.input_file, line_no, "delta line before first ABCD header", line, current_seq);
            }

            running_n = add_u32_checked(running_n, delta,
                                        sd.input_file, line_no, line, current_seq,
                                        "accumulating n += delta");

            lastN[current_seq] = running_n;
            nvec_push(&nlists[current_seq], running_n);
            continue;
        }

        // Must be ABCD header
        uint32_t hk = 0, hb = 0, hn = 0;
        int32_t  hc = 0;

        if (!parse_abcd_header_sr5(line, &hk, &hb, &hc, &hn)) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "invalid or unrecognised line", line, current_seq);
        }

        // Basic sanity: K>0, N0>0 (sr5sieve expects positive values)
        if (hk == 0) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "invalid K (must be > 0)", line, current_seq);
        }
        if (hn == 0) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "invalid N0 (must be > 0)", line, current_seq);
        }

        // Enforce B constant and valid
        if (hb < 2) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "invalid base in ABCD header (B must be >= 2)", line, current_seq);
        }
        if (global_B < 0) global_B = (int32_t)hb;
        else if ((uint32_t)global_B != hb) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "mismatched base (B must match all sequences)", line, current_seq);
        }

        // Standard form: c == +/-1
        if (hc != 1 && hc != -1) {
            std::fclose(fp);
            line_error(sd.input_file, line_no, "not standard form k*b^n+/-1 (c must be +1 or -1)", line, current_seq);
        }

        if (st.kcount >= MAX_SEQUENCES) {
            std::fclose(fp);
            char buf[256];
            std::snprintf(buf, sizeof(buf), "too many sequences (MAX_SEQUENCES=%d)", MAX_SEQUENCES);
            line_error(sd.input_file, line_no, buf, line, current_seq);
        }

        // New sequence starts here.  No per-sequence struct is used anymore.
        current_seq = (int)st.kcount++;

        st.klist[current_seq] = (hc > 0) ? (int)hk : -(int)hk;
        firstN[current_seq] = hn;
        lastN[current_seq] = hn;

        // Reset running_n and record N0
        running_n = hn;
        nvec_push(&nlists[current_seq], running_n);
    }

    if (std::ferror(fp)) {
        std::fclose(fp);
        std::fprintf(stderr, "%s:%" PRIu64 ": read error while reading ABCD file\n", sd.input_file, line_no);
        std::exit(EXIT_FAILURE);
    }
    std::fclose(fp);

    if (st.kcount == 0 || global_B < 0) {
        std::fprintf(stderr, "%s: no sequences found in input\n", sd.input_file);
        std::exit(EXIT_FAILURE);
    }

    // First pass: compute global NMIN/NMAX from all sequences.
    for (int i = 0; i < st.kcount; i++) {
        NVec *nv = &nlists[i];

        if (nv->n == 0) {
            std::fprintf(stderr, "%s: internal error: sequence %d has no N values\n", sd.input_file, i);
            std::exit(EXIT_FAILURE);
        }

        if (lastN[i] < firstN[i]) {
            std::fprintf(stderr, "%s: internal error: lastN < firstN for seq %d\n", sd.input_file, i);
            std::exit(EXIT_FAILURE);
        }

        if (firstN[i] < NMIN) NMIN = firstN[i];
        if (lastN[i]  > NMAX) NMAX = lastN[i];
    }

    // Bitmap NMIN must be even.  This may create one extra zero bit before the
    // smallest file N value, which is exactly what we want.
    if (NMIN & 1u) NMIN--;

    st.base = (uint32_t)global_B;
    st.nmin = NMIN;
    st.nmax = NMAX;

    const uint64_t nbits64 = (uint64_t)st.nmax - (uint64_t)st.nmin + 1u;
    if (nbits64 == 0 || nbits64 > (uint64_t)((size_t)-1)) {
        std::fprintf(stderr, "%s: bitmap range is too large\n", sd.input_file);
        std::exit(EXIT_FAILURE);
    }

    sd.bitmap_bits_per_sequence = (size_t)nbits64;
    sd.bitmap_words_per_sequence = (sd.bitmap_bits_per_sequence + 31u) >> 5;

    if (sd.bitmap_words_per_sequence != 0 &&
        (size_t)st.kcount > ((size_t)-1) / sd.bitmap_words_per_sequence) {
        std::fprintf(stderr, "%s: bitmap total size overflow\n", sd.input_file);
        std::exit(EXIT_FAILURE);
    }

    sd.bitmap_total_words = (size_t)st.kcount * sd.bitmap_words_per_sequence;
    sd.bitmap = init_bitmap_words_zero(sd.bitmap_total_words);

    // Second pass: set only actual N values read from the file.  Everything not
    // read, including positions below/above each sequence's own range, remains 0.
    for (int i = 0; i < st.kcount; i++) {
        NVec *nv = &nlists[i];

        for (size_t j = 0; j < nv->n; j++) {
            uint32_t nn = nv->v[j];
            if (nn < firstN[i] || nn > lastN[i]) {
                std::fprintf(stderr, "%s: internal error: N out of local sequence range while building bitmap (seq %d)\n",
                             sd.input_file, i);
                std::exit(EXIT_FAILURE);
            }
            set_can_use(sd, st, (uint32_t)i, nn);
        }

        std::free(nv->v);
        nv->v = nullptr;
        nv->n = nv->cap = 0;
    }

    std::printf("Sequences read: %d\n", (int)st.kcount);
    std::fprintf(stderr, "Sequences read: %d\n", (int)st.kcount);
    std::printf("Base = %" PRIu32 "\n", st.base);
    std::fprintf(stderr, "Base = %" PRIu32 "\n", st.base);

//    std::printf("NMIN = %" PRIu32 "\n", st.nmin);
//    std::printf("NMAX = %" PRIu32 "\n", st.nmax);
//    std::printf("Bitmap bits/seq = %zu\n", sd.bitmap_bits_per_sequence);
//    std::printf("Bitmap words/seq = %zu\n", sd.bitmap_words_per_sequence);
//    std::printf("Bitmap total words = %zu\n", sd.bitmap_total_words);
//
//    std::printf("klist array: ");
//    for (int i = 0; i < st.kcount; i++) std::printf("%d ", st.klist[i]);
//    std::printf("\n");
}
