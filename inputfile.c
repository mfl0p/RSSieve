
/*

	NOTE: this is incomplete and is broken!

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "simpleCL.h"
#include "cl_sieve.h"

// Skip whitespace
char* skip_ws(char *p) {
    while (*p && isspace((unsigned char)*p)) p++;
    return p;
}

// Parse long long
long long parse_ll(char **p) {
    errno = 0;
    char *end;
    long long val = strtoll(*p, &end, 10);
    if (errno != 0) { perror("strtoll"); exit(1); }
    *p = end;
    return val;
}

// Initialize bitmap with all bits cleared
uint8_t* init_bitmap_empty(size_t nbits) {
    size_t nbytes = (nbits + 7) / 8;
    uint8_t *bm = (uint8_t *)malloc(nbytes);
    if (!bm) { perror("malloc"); exit(1); }
    memset(bm, 0, nbytes);
    return bm;
}

// Set a bit as "can be used"
void set_can_use(uint8_t *bitmap, size_t offset) {
    bitmap[offset / 8] |= (1 << (offset % 8));
}

// Check if N can be used
int can_use_N(Sequence *seq, long long n) {
    if (n < seq->N0 || n > seq->lastN) return 0;
    size_t offset = (size_t)(n - seq->N0);
    return (seq->bitmap[offset / 8] >> (offset % 8)) & 1;
}

// Mark N as used
void mark_N_used(Sequence *seq, long long n) {
    if (n < seq->N0 || n > seq->lastN) return;
    size_t offset = (size_t)(n - seq->N0);
    seq->bitmap[offset / 8] &= ~(1 << (offset % 8));
}

// Find sequence by K and sign (B is global)
int find_sequence(Sequence *sequences, size_t count, long long K, char sign) {
    for (size_t i = 0; i < count; i++) {
        if (sequences[i].K == K && sequences[i].sign == sign) return (int)i;
    }
    return -1;
}

// Check if factor can be used for (K, sign, N)
int factor_can_be_used(Sequence *sequences, size_t count, long long K, char sign, long long n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx < 0) return 0;
    return can_use_N(&sequences[idx], n);
}

// Mark factor used for (K, sign, N)
void mark_factor_used(Sequence *sequences, size_t count, long long K, char sign, long long n) {
    int idx = find_sequence(sequences, count, K, sign);
    if (idx >= 0) mark_N_used(&sequences[idx], n);
}

void read_input(workStatus & st, searchData & sd) {

    FILE *fp = fopen(sd.input_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: unable to open input file %s\n", sd.input_file);
	printf("Error: unable to open input file %s\n", sd.input_file);
	exit(EXIT_FAILURE);
    }

    char line[512];
    long long global_B = -1;

    long long NMIN = -1, NMAX = -1;

    while (fgets(line, sizeof(line), fp)) {
        char *p = skip_ws(line);
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        if (strncmp(p, "ABCD", 4) == 0) {
            if (sd.kcount >= MAX_SEQUENCES) {
                fprintf(stderr, "Error: number of sequences exceeds MAX_SEQUENCES (%d)\n", MAX_SEQUENCES);
		printf("Error: number of sequences exceeds MAX_SEQUENCES (%d)\n", MAX_SEQUENCES);
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            p += 4;
            p = skip_ws(p);

            long long K = parse_ll(&p);
            p = skip_ws(p);

            if (*p != '*') { fprintf(stderr, "Expected '*'\n"); exit(1); }
            p++;
            p = skip_ws(p);

            long long B = parse_ll(&p);
            if (global_B == -1) global_B = B;
            else if (B != global_B) {
                fprintf(stderr, "Error: B must be the same for all sequences (found %lld, expected %lld)\n", B, global_B);
                exit(1);
            }

            p = skip_ws(p);
            if (strncmp(p, "^$a", 3) != 0) { fprintf(stderr, "Expected '^$a'\n"); exit(1); }
            p += 3;

            char sign = *p;
            if (sign != '+' && sign != '-') { fprintf(stderr, "Expected + or -\n"); exit(1); }
            p++;
            p = skip_ws(p);

            long long C = parse_ll(&p);
            if (C != 1) { fprintf(stderr, "Error: C must be 1, got %lld\n", C); exit(1); }
            p = skip_ws(p);

            long long N0;
            if (*p == '[') {
                p++;
                N0 = parse_ll(&p);
                if (*p != ']') { fprintf(stderr, "Expected ']'\n"); exit(1); }
                p++;
            } else { fprintf(stderr, "Error: N0 is required\n"); exit(1); }

            long long header_pos = ftell(fp);

            // First pass: compute lastN
            long long lastN = N0;
            long long running_N = N0;
            while (fgets(line, sizeof(line), fp)) {
                char *q = skip_ws(line);
                if (*q == '\0' || *q == '\n' || *q == '#') continue;
                if (strncmp(q, "ABCD", 4) == 0) { fseek(fp, -strlen(line), SEEK_CUR); break; }
                char *endptr;
                errno = 0;
                long long d = strtoll(q, &endptr, 10);
                if (errno != 0) { perror("strtoll"); exit(1); }
                running_N += d;
            }
            lastN = running_N;
            size_t nbits = (size_t)(lastN - N0 + 1);

            sd.sequences[sd.kcount].bitmap = init_bitmap_empty(nbits);
            sd.sequences[sd.kcount].K = K;
            sd.sequences[sd.kcount].sign = sign;
	    sd.sequences[sd.kcount].N0 = N0;
            sd.sequences[sd.kcount].lastN = lastN;
            sd.sequences[sd.kcount].nbits = nbits;

            sd.klist[sd.kcount] = (sign == '+') ? (int)K : -(int)K;

            // Second pass: mark bitmap from offsets
            fseek(fp, header_pos, SEEK_SET);
            running_N = N0;
            while (fgets(line, sizeof(line), fp)) {
                char *q = skip_ws(line);
                if (*q == '\0' || *q == '\n' || *q == '#') continue;
                if (strncmp(q, "ABCD", 4) == 0) { fseek(fp, -strlen(line), SEEK_CUR); break; }
                char *endptr;
                errno = 0;
                long long d = strtoll(q, &endptr, 10);
                if (errno != 0) { perror("strtoll"); exit(1); }
                running_N += d;
                set_can_use(sd.sequences[sd.kcount].bitmap, (size_t)(running_N - N0));
            }

            if (NMIN == -1 || N0 < NMIN) NMIN = N0;
            if (NMAX == -1 || lastN > NMAX) NMAX = lastN;

            sd.kcount++;
        }
    }

    fclose(fp);

    if (NMIN % 2 != 0) NMIN--;

    st.nmin = NMIN;
    st.nmax = NMAX;
    st.base = (uint32_t)global_B;

    printf("Sequences read: %d\n", sd.kcount);
    printf("Base = %lld\n", global_B);
    printf("NMIN = %lld\n", NMIN);
    printf("NMAX = %lld\n", NMAX);

    printf("klist array: ");
    for (int i = 0; i < sd.kcount; i++) printf("%d ", sd.klist[i]);
    printf("\n");


   // for (size_t i = 0; i < sd.kcount; i++) free(sequences[i].bitmap);

}

