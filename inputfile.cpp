/*
	THIS IS STILL BROKEN

	read ABCD input file

      - Assumes sd.sequences[] and sd.klist[] are allocated and MAX_SEQUENCES exists.
      - Assumes Sequence has fields: bitmap, K, sign, N0, lastN, nbits.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "simpleCL.h"
#include "cl_sieve.h"

// ---- helpers ----

static char* skip_ws(char *p) {
    while (*p && isspace((unsigned char)*p)) p++;
    return p;
}

static void strip_crlf(char *s) {
    // Remove trailing \r and/or \n so Windows CRLF doesn’t leave a '\r' at end.
    size_t n = strlen(s);
    while (n && (s[n-1] == '\n' || s[n-1] == '\r')) {
        s[n-1] = '\0';
        --n;
    }
}

static long long parse_ll(char **p) {
    errno = 0;
    char *end = NULL;
    long long val = strtoll(*p, &end, 10);
    if (end == *p) {
        fprintf(stderr, "Error: expected integer near: '%s'\n", *p);
        exit(EXIT_FAILURE);
    }
    if (errno != 0) {
        perror("strtoll");
        exit(EXIT_FAILURE);
    }
    *p = end;
    return val;
}

static uint8_t* init_bitmap_empty(size_t nbits) {
    size_t nbytes = (nbits + 7) / 8;
    uint8_t *bm = (uint8_t*)malloc(nbytes);
    if (!bm) { perror("malloc"); exit(EXIT_FAILURE); }
    memset(bm, 0, nbytes);
    return bm;
}

// Set a bit as "can be used" (1 = can use)
static void set_can_use(uint8_t *bitmap, size_t offset) {
    bitmap[offset / 8] |= (uint8_t)(1u << (offset % 8));
}

// Clear a bit as "used" (0 = used / cannot use)
static void clear_can_use(uint8_t *bitmap, size_t offset) {
    bitmap[offset / 8] &= (uint8_t)~(1u << (offset % 8));
}

static int can_use_N(const Sequence *seq, long long n) {
    if (n < seq->N0 || n > seq->lastN) return 0;
    size_t offset = (size_t)(n - seq->N0);
    if (offset >= seq->nbits) return 0;
    return (seq->bitmap[offset / 8] >> (offset % 8)) & 1;
}

static void mark_N_used(Sequence *seq, long long n) {
    if (n < seq->N0 || n > seq->lastN) return;
    size_t offset = (size_t)(n - seq->N0);
    if (offset >= seq->nbits) return;
    clear_can_use(seq->bitmap, offset);
}

// Find sequence by K and sign (B is global)
static int find_sequence(const Sequence *sequences, size_t count, long long K, char sign) {
    for (size_t i = 0; i < count; i++) {
        if (sequences[i].K == K && sequences[i].sign == sign) return (int)i;
    }
    return -1;
}

int factor_can_be_used(Sequence *sequences, size_t count, long long K, char sign, long long n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx < 0) return 0;
    return can_use_N(&sequences[idx], n);
}

void mark_factor_used(Sequence *sequences, size_t count, long long K, char sign, long long n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx >= 0) mark_N_used(&sequences[idx], n);
}

// Parse an ABCD header line.
// Expects: "ABCD <K>*<B>^$a<+/-C> [<N0>] ..."
static void parse_abcd_header(
    char *line,
    long long *outK,
    long long *outB,
    char *outSign,
    long long *outC,
    long long *outN0
) {
    strip_crlf(line);

    char *p = skip_ws(line);
    if (strncmp(p, "ABCD", 4) != 0) {
        fprintf(stderr, "Internal error: parse_abcd_header called on non-ABCD line\n");
        exit(EXIT_FAILURE);
    }

    p += 4;
    p = skip_ws(p);

    long long K = parse_ll(&p);
    p = skip_ws(p);

    if (*p != '*') { fprintf(stderr, "Error: Expected '*'\n"); exit(EXIT_FAILURE); }
    p++;
    p = skip_ws(p);

    long long B = parse_ll(&p);
    p = skip_ws(p);

    if (strncmp(p, "^$a", 3) != 0) { fprintf(stderr, "Error: Expected '^$a'\n"); exit(EXIT_FAILURE); }
    p += 3;

    char sign = *p;
    if (sign != '+' && sign != '-') { fprintf(stderr, "Error: Expected + or -\n"); exit(EXIT_FAILURE); }
    p++;
    p = skip_ws(p);

    long long C = parse_ll(&p);
    p = skip_ws(p);

    if (*p != '[') { fprintf(stderr, "Error: N0 is required\n"); exit(EXIT_FAILURE); }
    p++;
    long long N0 = parse_ll(&p);
    if (*p != ']') { fprintf(stderr, "Error: Expected ']'\n"); exit(EXIT_FAILURE); }
    // p++ optional; we don't need to consume further.

    *outK = K;
    *outB = B;
    *outSign = sign;
    *outC = C;
    *outN0 = N0;
}

// Reads offset lines until next ABCD header or EOF.
// Returns lastN. Leaves file positioned at the start of the next ABCD line (or EOF).
static long long scan_offsets_compute_lastN(FILE *fp, long long N0) {
    char line[512];
    long long running_N = N0;

    while (1) {
        fpos_t pos;
        if (fgetpos(fp, &pos) != 0) { perror("fgetpos"); exit(EXIT_FAILURE); }

        if (!fgets(line, sizeof(line), fp)) break;

        strip_crlf(line);
        char *q = skip_ws(line);

        if (*q == '\0' || *q == '#') continue;

        if (strncmp(q, "ABCD", 4) == 0) {
            // Rewind to the start of this ABCD line for the caller.
            if (fsetpos(fp, &pos) != 0) { perror("fsetpos"); exit(EXIT_FAILURE); }
            break;
        }

        errno = 0;
        char *endptr = NULL;
        long long d = strtoll(q, &endptr, 10);
        if (endptr == q) {
            fprintf(stderr, "Error: invalid offset line: '%s'\n", q);
            exit(EXIT_FAILURE);
        }
        if (errno != 0) { perror("strtoll"); exit(EXIT_FAILURE); }

        running_N += d;
    }

    return running_N;
}

// Second pass: mark bitmap bits for N0 and each cumulative offset.
// Leaves file positioned at the start of the next ABCD line (or EOF).
static void scan_offsets_mark_bitmap(FILE *fp, uint8_t *bitmap, size_t nbits, long long N0) {
    char line[512];
    long long running_N = N0;

    // Mark N0 itself as usable.
    if (nbits > 0) set_can_use(bitmap, 0);

    while (1) {
        fpos_t pos;
        if (fgetpos(fp, &pos) != 0) { perror("fgetpos"); exit(EXIT_FAILURE); }

        if (!fgets(line, sizeof(line), fp)) break;

        strip_crlf(line);
        char *q = skip_ws(line);

        if (*q == '\0' || *q == '#') continue;

        if (strncmp(q, "ABCD", 4) == 0) {
            if (fsetpos(fp, &pos) != 0) { perror("fsetpos"); exit(EXIT_FAILURE); }
            break;
        }

        errno = 0;
        char *endptr = NULL;
        long long d = strtoll(q, &endptr, 10);
        if (endptr == q) {
            fprintf(stderr, "Error: invalid offset line: '%s'\n", q);
            exit(EXIT_FAILURE);
        }
        if (errno != 0) { perror("strtoll"); exit(EXIT_FAILURE); }

        running_N += d;

        long long off_ll = running_N - N0;
        if (off_ll < 0) {
            fprintf(stderr, "Error: negative offset encountered (running_N < N0)\n");
            exit(EXIT_FAILURE);
        }
        size_t off = (size_t)off_ll;
        if (off >= nbits) {
            fprintf(stderr, "Error: offset %zu out of range (nbits=%zu)\n", off, nbits);
            exit(EXIT_FAILURE);
        }
        set_can_use(bitmap, off);
    }
}

// ---- main function ----

void read_input(workStatus &st, searchData &sd) {
    FILE *fp = fopen(sd.input_file, "r"); // text mode ok because we use fgetpos/fsetpos
    if (!fp) {
        fprintf(stderr, "Error: unable to open input file %s\n", sd.input_file);
        printf("Error: unable to open input file %s\n", sd.input_file);
        exit(EXIT_FAILURE);
    }

    char line[512];
    long long global_B = -1;
    long long NMIN = -1, NMAX = -1;

    // Ensure kcount starts sane if caller didn’t.
    // (Remove if you intentionally accumulate.)
    // sd.kcount = 0;

    while (1) {
        fpos_t header_line_pos;
        if (fgetpos(fp, &header_line_pos) != 0) { perror("fgetpos"); exit(EXIT_FAILURE); }

        if (!fgets(line, sizeof(line), fp)) break;

        strip_crlf(line);
        char *p = skip_ws(line);
        if (*p == '\0' || *p == '#') continue;

        if (strncmp(p, "ABCD", 4) != 0) continue;

        if (sd.kcount >= MAX_SEQUENCES) {
            fprintf(stderr, "Error: number of sequences exceeds MAX_SEQUENCES (%d)\n", MAX_SEQUENCES);
            printf("Error: number of sequences exceeds MAX_SEQUENCES (%d)\n", MAX_SEQUENCES);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        long long K, B, C, N0;
        char sign;

        parse_abcd_header(p, &K, &B, &sign, &C, &N0);

        if (C != 1) {
            fprintf(stderr, "Error: C must be 1, got %lld\n", C);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        if (global_B == -1) global_B = B;
        else if (B != global_B) {
            fprintf(stderr, "Error: B must be the same for all sequences (found %lld, expected %lld)\n",
                    B, global_B);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Save position right AFTER the ABCD header line
        fpos_t after_header_pos;
        if (fgetpos(fp, &after_header_pos) != 0) { perror("fgetpos"); exit(EXIT_FAILURE); }

        // First pass: compute lastN by scanning offsets until next ABCD
        long long lastN = scan_offsets_compute_lastN(fp, N0);

        if (lastN < N0) {
            fprintf(stderr, "Error: lastN < N0 (corrupt offsets?)\n");
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        size_t nbits = (size_t)(lastN - N0 + 1);
        if (nbits == 0) {
            fprintf(stderr, "Error: computed nbits == 0\n");
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        Sequence *seq = &sd.sequences[sd.kcount];
        seq->bitmap = init_bitmap_empty(nbits);
        seq->K = K;
        seq->sign = sign;
        seq->N0 = N0;
        seq->lastN = lastN;
        seq->nbits = nbits;

        sd.klist[sd.kcount] = (sign == '+') ? (int)K : -(int)K;

        // Second pass: rewind to after-header and mark bitmap
        if (fsetpos(fp, &after_header_pos) != 0) { perror("fsetpos"); exit(EXIT_FAILURE); }
        scan_offsets_mark_bitmap(fp, seq->bitmap, nbits, N0);

        if (NMIN == -1 || N0 < NMIN) NMIN = N0;
        if (NMAX == -1 || lastN > NMAX) NMAX = lastN;

        sd.kcount++;
        // Loop continues from current fp position (already at next ABCD or EOF).
    }

    fclose(fp);

    if (sd.kcount == 0) {
        fprintf(stderr, "Error: no sequences found in input.\n");
        exit(EXIT_FAILURE);
    }
    if (global_B < 0) {
        fprintf(stderr, "Error: global B not found.\n");
        exit(EXIT_FAILURE);
    }

    if (NMIN % 2 != 0) NMIN--;

    st.nmin = (uint32_t)NMIN;
    st.nmax = (uint32_t)NMAX;
    st.base = (uint32_t)global_B;

    printf("Sequences read: %d\n", sd.kcount);
    printf("Base = %lld\n", global_B);
    printf("NMIN = %lld\n", NMIN);
    printf("NMAX = %lld\n", NMAX);

    printf("klist array: ");
    for (int i = 0; i < sd.kcount; i++) printf("%d ", sd.klist[i]);
    printf("\n");

    // Caller can free later:
    // for (int i=0; i<sd.kcount; ++i) free(sd.sequences[i].bitmap);
}
