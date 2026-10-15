#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/*
All W*H*D cells are stored in one contiguous array ID[0..N-1].
k varies fastest, then j, then i: n = (i*H + j)*D + k = i·(H·D) + j·D + k
Strides: k -> 1, j -> D, i -> H*D elements.
*/
//atoms will be in places in the lattice where x + y + z is even
//own layer(±1, ±1, 0)
//layer above(±1, 0, +1), (0, ±1, +1)
//layer below(±1, 0, −1), (0, ±1, −1)
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

/* FCC embedding */
/*
 * Physical coordinates are integers in units of a/2 (a = cubic lattice
 * constant). FCC sites are exactly the integer points with x+y+z even.
 * Each z-layer is a checkerboard: in every row (fixed y,z) the sites sit
 * at every second x, shifted by the parity p = (y+z) & 1.
 *
 *   array -> physical:  x = 2*i + ((j+k) & 1),  y = j,  z = k
 *   physical -> array:  i = x / 2,              j = y,  k = z
 *
 * Physical periodic box: PX x PY x PZ = 2W x H x D, containing exactly
 * W*H*D sites = N cells (bijection). H and D must be even so that
 * periodic wrapping preserves the parity of x+y+z.
 */
#define PX (2 * W)
#define PY H
#define PZ D


typedef struct { int x, y, z; } Phys;

/* array (i,j,k) -> physical (x,y,z) */
static inline Phys to_phys(int i, int j, int k)
{
    Phys r;
    r.x = 2 * i + ((j + k) & 1);
    r.y = j;
    r.z = k;
    return r;
}

static inline Coord to_array(int x, int y, int z)
{
    assert(((x + y + z) & 1) == 0);
    Coord c = { x / 2, y, z };
    return c;
}

static inline int wrap(int v, int L) { return ((v % L) + L) % L; }

static const int DIR[12][3] = {
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0},   /* own layer   */
    { 1, 0, 1}, {-1, 0,-1}, { 0, 1, 1}, { 0,-1,-1},   /* above/below */
    { 1, 0,-1}, {-1, 0, 1}, { 0, 1,-1}, { 0,-1, 1},   /* below/above */
};

/*verification*/

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

static int verify_embedding(void)
{
    int errors = 0;

    /*every cell -> a valid site in the box -> same cell */
    for (int n = 0; n < N; n++) {
        Coord c = coords(n);
        Phys  p = to_phys(c.i, c.j, c.k);

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

    /*every site in the box -> a cell -> same site*/
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

    /* closure: every site + every NN vector, wrapped, is again a site */
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

    printf("embedding:     %s (%d errors)\n", errors ? "FAILED" : "OK", errors);
    return errors;
}


static int verify_neighbors(void)
{
    int errors = 0;
 
    for (int n = 0; n < N; n++) {
        int nb[12];
 
        for (int d = 0; d < 12; d++) {
            nb[d] = neighbor(n, d);
 
            /* valid index, and not the cell itself */
            if (nb[d] < 0 || nb[d] >= N || nb[d] == n) {
                fprintf(stderr, "cell %d: bad neighbor %d in direction %d\n", n, nb[d], d);
                errors++;
                continue;
            }
 
            /* symmetry: going back in the opposite direction returns to n */
            if (neighbor(nb[d], d ^ 1) != n) {
                fprintf(stderr, "cell %d: direction %d not symmetric\n", n, d);
                errors++;
            }
        }
 
        /* the 12 neighbours must be distinct */
        for (int a = 0; a < 12; a++)
            for (int b = a + 1; b < 12; b++)
                if (nb[a] == nb[b]) {
                    fprintf(stderr, "cell %d: directions %d and %d give the same cell %d\n",
                            n, a, b, nb[a]);
                    errors++;
                }
    }
 
    printf("neighbors:     %s (%d errors)\n", errors ? "FAILED" : "OK", errors);
    return errors;
}


/* Task 4 */
 
/* allocate the N cells in one contiguous block and set ID[n] = n */
static int lattice_init(void)
{
    ID = malloc((size_t)N * sizeof *ID);
    if (!ID) {
        perror("malloc");
        return -1;
    }
    for (int n = 0; n < N; n++)
        ID[n] = n;
    return 0;
}
 
static void lattice_free(void)
{
    free(ID);
    ID = NULL;
}

static inline int neighbor(int n, int d)
{
    Coord c = coords(n);
    Phys  p = to_phys(c.i, c.j, c.k);
    int x = wrap(p.x + DIR[d][0], PX);
    int y = wrap(p.y + DIR[d][1], PY);
    int z = wrap(p.z + DIR[d][2], PZ);
    Coord m = to_array(x, y, z);
    return idx(m.i, m.j, m.k);
}


int main(void)
{
    if (lattice_init() != 0) return EXIT_FAILURE;

    int errors = 0;
    errors += verify_index_mapping();
    errors += verify_embedding();
    errors += verify_neighbors();

    lattice_free();
    return errors ? EXIT_FAILURE : EXIT_SUCCESS;
}