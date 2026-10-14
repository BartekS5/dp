#include<stdio.h>
#include <stdlib.h>
/*
All W*H*D cells are stored in one contiguous array ID[0..N-1].
k varies fastest, then j, then i: n = (i*H + j)*D + k = i·(H·D) + j·D + k
Strides: k -> 1, j -> D, i -> H*D elements.
*/
#define W 32
#define H 32
#define D 32
#define N (W * H * D)

static int *ID; 

typedef struct { int i, j, k; } Coord;

static inline int idx(int i, int j, int k)
{
    return (i * H + j) * D + k;
}

static inline Coord coords(int n)
{
    Coord c;
    c.k = n % D;
    c.j = (n / D) % H;
    c.i = n / (H * D);
    return c;
}

int main() {
    printf("Hello, World!\n");
    return 0;
}

static int verify_index_mapping(void)
{
    int errors = 0;

    for (int n = 0; n < N; n++) {
        Coord c = coords(n);
        if (idx(c.i, c.j, c.k) != n) {
            fprintf(stderr, "idx(coords(%d)) != %d\n", n, n);
            errors++;
        }
    }

    for (int i = 0; i < W; i++)
        for (int j = 0; j < H; j++)
            for (int k = 0; k < D; k++) {
                Coord c = coords(idx(i, j, k));
                if (c.i != i || c.j != j || c.k != k) {
                    fprintf(stderr, "coords(idx(%d,%d,%d)) wrong\n", i, j, k);
                    errors++;
                }
            }

    printf("index mapping: %s (%d errors)\n", errors ? "FAILED" : "OK", errors);
    return errors;
}