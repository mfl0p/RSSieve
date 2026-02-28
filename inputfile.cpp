//
// Even closer to sr5sieve behavior, with extra strictness:
//  - filename + line number on ALL errors
//  - includes current sequence index (when known) in relevant errors
//  - detects and errors on uint32_t overflow when accumulating n += delta
//  - detects and errors on overflow when computing lastN = N0 + sum(deltas)
//  - validates K,B,N0 are nonzero-ish where sensible (B>=2, K>=1, N0>=1)
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
// Output:
//  - Builds seq bitmap bits for exactly the N values present in the file
//  - Sets st.base, st.nmin, st.nmax (and if nmin odd => nmin--)

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <inttypes.h>

#include "simpleCL.h"
#include "cl_sieve.h"

// -----------------------------------------------------------------------------
// Bitmap helpers + factor usage helpers
// -----------------------------------------------------------------------------

static uint8_t* init_bitmap_empty(size_t nbits) {
    size_t nbytes = (nbits + 7) / 8;
    uint8_t *bm = (uint8_t*)std::malloc(nbytes);
    if (!bm) { std::perror("malloc"); std::exit(EXIT_FAILURE); }
    std::memset(bm, 0, nbytes);
    return bm;
}

static inline void set_can_use(uint8_t *bitmap, size_t offset) {
    bitmap[offset / 8] |= (uint8_t)(1u << (offset % 8));
}

static inline void clear_can_use(uint8_t *bitmap, size_t offset) {
    bitmap[offset / 8] &= (uint8_t)~(1u << (offset % 8));
}

static int can_use_N(const Sequence *seq, uint32_t n) {
    if (n < seq->N0 || n > seq->lastN) return 0;
    size_t offset = (size_t)(n - seq->N0);
    if (offset >= seq->nbits) return 0;
    return (seq->bitmap[offset / 8] >> (offset % 8)) & 1;
}

static void mark_N_used(Sequence *seq, uint32_t n) {
    if (n < seq->N0 || n > seq->lastN) return;
    size_t offset = (size_t)(n - seq->N0);
    if (offset >= seq->nbits) return;
    clear_can_use(seq->bitmap, offset);
}

static int find_sequence(const Sequence *sequences, size_t count, uint32_t K, char sign) {
    for (size_t i = 0; i < count; i++) {
        if (sequences[i].K == K && sequences[i].sign == sign) return (int)i;
    }
    return -1;
}

int factor_can_be_used(Sequence *sequences, size_t count, uint32_t K, char sign, uint32_t n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx < 0) return 0;
    return can_use_N(&sequences[idx], n);
}

void mark_factor_used(Sequence *sequences, size_t count, uint32_t K, char sign, uint32_t n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx >= 0) mark_N_used(&sequences[idx], n);
}

void free_sequences(searchData &sd)
{
    for (int i = 0; i < sd.kcount; ++i) {
        if (sd.sequences[i].bitmap) {
            free(sd.sequences[i].bitmap);
            sd.sequences[i].bitmap = nullptr;
            sd.sequences[i].nbits = 0;
        }
    }
    sd.kcount = 0;
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

    FILE *fp = std::fopen(sd.input_file, "r");
    if (!fp) file_error_open(sd.input_file);

    NVec nlists[MAX_SEQUENCES];
    std::memset(nlists, 0, sizeof(nlists));

    sd.kcount = 0;

    char line[2048];
    uint64_t line_no = 0;

    int current_seq = -1;      // no active sequence until first ABCD header
    uint32_t running_n = 0;    // current n for delta accumulation

    int32_t global_B = -1;
    uint32_t NMIN = 0xFFFFFFFFu;
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

        if (sd.kcount >= MAX_SEQUENCES) {
            std::fclose(fp);
            char buf[256];
            std::snprintf(buf, sizeof(buf), "too many sequences (MAX_SEQUENCES=%d)", MAX_SEQUENCES);
            line_error(sd.input_file, line_no, buf, line, current_seq);
        }

        // New sequence starts here
        current_seq = (int)sd.kcount++;
        Sequence *seq = &sd.sequences[current_seq];
        std::memset(seq, 0, sizeof(*seq));

        seq->K    = hk;
        seq->sign = (hc > 0) ? '+' : '-';
        seq->N0   = hn;

        sd.klist[current_seq] = (seq->sign == '+') ? (int)hk : -(int)hk;

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

    if (sd.kcount == 0 || global_B < 0) {
        std::fprintf(stderr, "%s: no sequences found in input\n", sd.input_file);
        std::exit(EXIT_FAILURE);
    }

    // Build bitmaps + compute lastN/NMIN/NMAX
    for (int i = 0; i < sd.kcount; i++) {
        Sequence *seq = &sd.sequences[i];
        NVec *nv = &nlists[i];

        if (nv->n == 0) {
            std::fprintf(stderr, "%s: internal error: sequence %d has no N values\n", sd.input_file, i);
            std::exit(EXIT_FAILURE);
        }

        const uint32_t n0 = seq->N0;
        const uint32_t lastN = nv->v[nv->n - 1];
        seq->lastN = lastN;

        if (n0 < NMIN) NMIN = n0;
        if (lastN > NMAX) NMAX = lastN;

        if (lastN < n0) {
            std::fprintf(stderr, "%s: internal error: lastN < N0 for seq %d\n", sd.input_file, i);
            std::exit(EXIT_FAILURE);
        }

        // nbits = (lastN - n0 + 1) in size_t; lastN>=n0 so safe
        seq->nbits = (size_t)(lastN - n0 + 1);
        if (seq->nbits == 0) {
            std::fprintf(stderr, "%s: internal error: computed nbits == 0 for seq %d\n", sd.input_file, i);
            std::exit(EXIT_FAILURE);
        }

        seq->bitmap = init_bitmap_empty(seq->nbits);

        for (size_t j = 0; j < nv->n; j++) {
            uint32_t nn = nv->v[j];
            if (nn < n0 || nn > lastN) {
                std::fprintf(stderr, "%s: internal error: N out of range while building bitmap (seq %d)\n",
                             sd.input_file, i);
                std::exit(EXIT_FAILURE);
            }
            size_t off = (size_t)(nn - n0);
            if (off >= seq->nbits) {
                std::fprintf(stderr, "%s: internal error: bitmap offset out of range (seq %d)\n",
                             sd.input_file, i);
                std::exit(EXIT_FAILURE);
            }
            set_can_use(seq->bitmap, off);
        }

        std::free(nv->v);
        nv->v = nullptr;
        nv->n = nv->cap = 0;
    }

    // If NMIN is odd, reduce by 1
    if (NMIN & 1u) NMIN--;

    st.base = (uint32_t)global_B;
    st.nmin = NMIN;
    st.nmax = NMAX;

    std::printf("Sequences read: %d\n", (int)sd.kcount);
    std::printf("Base = %" PRIu32 "\n", st.base);
    std::printf("NMIN = %" PRIu32 "\n", st.nmin);
    std::printf("NMAX = %" PRIu32 "\n", st.nmax);

    std::printf("klist array: ");
    for (int i = 0; i < sd.kcount; i++) std::printf("%d ", sd.klist[i]);
    std::printf("\n");
}
