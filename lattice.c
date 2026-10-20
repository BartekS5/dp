#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define W 256
#define H 256
#define D 256
#define N (W * H * D)

int *lattice;

typedef struct{
    int i, j, k;
} Coord;

int idx(int i, int j, int k)
{
    return i * (H * D) + j * D + k;
}

Coord coords(int n)
{
    Coord c;
    c.k = n % D;
    c.j = (n / D) % H;
    c.i = n / (H * D);
    return c;
}

#define PX (2 * W)
#define PY H
#define PZ D


typedef struct{
    int x, y, z;
} Phys;

Phys to_phys(int i, int j, int k)
{
    Phys r;
    r.x = 2 * i + ((j + k) & 1);
    r.y = j;
    r.z = k;
    return r;
}

Coord to_array(int x, int y, int z)
{
    assert(((x + y + z) & 1) == 0);
    Coord c = { x / 2, y, z };
    return c;
}

int wrap(int v, int L) { return ((v % L) + L) % L; }

const int DIR[12][3] = {
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0},
    { 1, 0, 1}, {-1, 0,-1}, { 0, 1, 1}, { 0,-1,-1},
    { 1, 0,-1}, {-1, 0, 1}, { 0, 1,-1}, { 0,-1, 1},
};

/* Task 4 */
 
int lattice_init(void)
{
    lattice = malloc((size_t)N * sizeof(*lattice));
    if (!lattice) {
        perror("malloc error");
        return -1;
    }
    for (int n = 0; n < N; n++)
        lattice[n] = n;
    return 0;
}
 
void lattice_free(void)
{
    free(lattice);
    lattice = NULL;
}

int neighbor(int n, int d)
{
    Coord c = coords(n);
    Phys  p = to_phys(c.i, c.j, c.k);
    int x = wrap(p.x + DIR[d][0], PX);
    int y = wrap(p.y + DIR[d][1], PY);
    int z = wrap(p.z + DIR[d][2], PZ);
    Coord m = to_array(x, y, z);
    return idx(m.i, m.j, m.k);
}

/* Task 5 */

volatile long long sum_kernel()
{
    long long sum = 0;
    for (int n = 0; n < N; n++)
        sum += lattice[n];
    return sum;
}

/* Task 6 */

long long sum_inner_k()
{
    clock_t start = clock();
    long long sum = 0;
    for (int i = 0; i < W; i++)
        for (int j = 0; j < H; j++)
            for (int k = 0; k < D; k++) {
                sum += lattice[idx(i, j, k)];
            }
    clock_t end = clock();
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("sum_inner_k: %lld, elapsed time: %f s\n", sum, elapsed_time);
    return sum;
}

long long sum_inner_j()
{
    clock_t start = clock();
    long long sum = 0;
    for (int i = 0; i < H; i++)
        for (int k = 0; k < W; k++)
            for (int j = 0; j < D; j++) {
                sum += lattice[idx(i, j, k)];
            }
    clock_t end = clock();
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("sum_inner_j: %lld, elapsed time: %f s\n", sum, elapsed_time);
    return sum;
}

long long sum_inner_i()
{
    clock_t start = clock();
    long long sum = 0;
    for (int k = 0; k < D; k++)
        for (int j = 0; j < H; j++)
            for (int i = 0; i < W; i++) {
                sum += lattice[idx(i, j, k)];
            }
    clock_t end = clock();
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("sum_inner_i: %lld, elapsed time: %f s\n", sum, elapsed_time);
    return sum;
}

/*verification*/

int verify_index_mapping(void)
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

int verify_embedding(void)
{
    int errors = 0;

    for (int n = 0; n < N; n++) {
        Coord c = coords(n);
        Phys p = to_phys(c.i, c.j, c.k);

        if ((p.x + p.y + p.z) & 1) {
            fprintf(stderr, "cell %d maps to non-site (%d,%d,%d)\n", n, p.x, p.y, p.z);
            errors++;
            continue;
        }
        if (p.x < 0 || p.x >= PX || p.y < 0 || p.y >= PY || p.z < 0 || p.z >= PZ) {
            fprintf(stderr, "cell %d maps outside the box\n", n);
            errors++;
            continue;
        }
        Coord b = to_array(p.x, p.y, p.z);
        if (b.i != c.i || b.j != c.j || b.k != c.k) {
            fprintf(stderr, "round trip failed for cell %d\n", n);
            errors++;
        }
    }

    int sites = 0;
    for (int x = 0; x < PX; x++)
        for (int y = 0; y < PY; y++)
            for (int z = 0; z < PZ; z++) {
                if ((x + y + z) & 1) continue;
                sites++;
                Coord c = to_array(x, y, z);
                Phys  p = to_phys(c.i, c.j, c.k);
                if (p.x != x || p.y != y || p.z != z) {
                    fprintf(stderr, "site (%d,%d,%d) round trip failed\n", x, y, z);
                    errors++;
                }
            }
    if (sites != N) {
        fprintf(stderr, "box has %d sites, expected %d\n", sites, N);
        errors++;
    }

    for (int n = 0; n < N; n++) {
        Coord c = coords(n);
        Phys  p = to_phys(c.i, c.j, c.k);
        for (int d = 0; d < 12; d++) {
            int x = wrap(p.x + DIR[d][0], PX);
            int y = wrap(p.y + DIR[d][1], PY);
            int z = wrap(p.z + DIR[d][2], PZ);
            if ((x + y + z) & 1) {
                fprintf(stderr, "neighbor %d of cell %d is not a site\n", d, n);
                errors++;
            }
        }
    }

    printf("embedding: %s (%d errors)\n", errors ? "FAILED" : "OK", errors);
    return errors;
}


int verify_neighbors(void)
{
    int errors = 0;
 
    for (int n = 0; n < N; n++) {
        int nb[12];
 
        for (int d = 0; d < 12; d++) {
            nb[d] = neighbor(n, d);
 
            if (nb[d] < 0 || nb[d] >= N || nb[d] == n) {
                fprintf(stderr, "cell %d: bad neighbor %d in direction %d\n", n, nb[d], d);
                errors++;
                continue;
            }
 
            if (neighbor(nb[d], d ^ 1) != n) {
                fprintf(stderr, "cell %d: direction %d not symmetric\n", n, d);
                errors++;
            }
        }
 
        for (int a = 0; a < 12; a++)
            for (int b = a + 1; b < 12; b++)
                if (nb[a] == nb[b]) {
                    fprintf(stderr, "cell %d: directions %d and %d give the same cell %d\n", n, a, b, nb[a]);
                    errors++;
                }
    }
 
    printf("neighbors:     %s (%d errors)\n", errors ? "FAILED" : "OK", errors);
    return errors;
}


int main(void)
{
    if (lattice_init() != 0) return EXIT_FAILURE;

    int errors = 0;
    errors += verify_index_mapping();
    errors += verify_embedding();
    errors += verify_neighbors();

    long long sum = sum_kernel();
    long long expected = (long long)N * (N - 1) / 2;
    printf("kernel sum: %lld (expected %lld) %s\n", sum, expected, sum == expected ? "OK" : "FAILED");
    if (sum != expected)
        errors++;

    sum_inner_k();
    sum_inner_j();
    sum_inner_i();

    lattice_free();
    return errors ? EXIT_FAILURE : EXIT_SUCCESS;
}