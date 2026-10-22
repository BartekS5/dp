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

typedef struct{
    int x, y, z;
} Phys;

int idx(int i, int j, int k){
    assert(i >= 0 && i < W);
    assert(j >= 0 && j < H);
    assert(k >= 0 && k < D);

    return i * (H * D) + j * D + k;
}

Coord coords(int n){
    assert(n >= 0 && n < N);

    Coord c;
    c.k = n % D;
    c.j = (n / D) % H;
    c.i = n / (H * D);

    return c;
}

Phys to_phys(int i, int j, int k){
    assert(i >= 0 && i < W);
    assert(j >= 0 && j < H);
    assert(k >= 0 && k < D);

    Phys r;
    r.x = i;
    r.y = j;
    r.z = 2 * k + ((i + j) & 1);

    assert(((r.x + r.y + r.z) & 1) == 0);
    return r;
}

Coord to_array(int x, int y, int z){
    assert(x >= 0 && x < W);
    assert(y >= 0 && y < H);
    assert(z >= 0 && z < 2 * D);
    assert(((x + y + z) & 1) == 0);

    Coord c;
    c.i = x;
    c.j = y;
    c.k = z / 2;

    return c;
}

int wrap(int a, int m){
    assert(m > 0);
    int r = a % m;
    return r < 0 ? r + m : r;
}

const int DIR[12][3] = {
    {1, 1, 0}, {-1, -1, 0}, {1, -1, 0}, { -1, 1, 0},
    {1, 0, 1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, 1},
    {1, 0, -1}, {-1, 0, -1}, {0, 1, -1}, {0, -1, -1},
};

const int OPP[12] = {1, 0, 3, 2, 9, 8, 11, 10, 5, 4, 7, 6};

/* Task 4 */
 
int lattice_init(void){
    lattice = malloc((size_t)N * sizeof(*lattice));
    if (!lattice) {
        perror("malloc error");
        return -1;
    }
    for (int n = 0; n < N; n++){
        lattice[n] = n;
    }
    return 0;
}
 
void lattice_free(void){
    free(lattice);
    lattice = NULL;
}

int neighbor(int n, int d){
    assert(n >= 0 && n < N);
    assert(d >= 0 && d < 12);

    Coord c = coords(n);
    Phys p = to_phys(c.i, c.j, c.k);

    int x = wrap(p.x + DIR[d][0], W);
    int y = wrap(p.y + DIR[d][1], H);
    int z = wrap(p.z + DIR[d][2], 2 * D);

    Coord m = to_array(x, y, z);

    return idx(m.i, m.j, m.k);
}

/* Task 5 */

long long sum_kernel(){
    long long sum = 0;
    for (int n = 0; n < N; n++){
        sum += lattice[n];
    }
    return sum;
}

/* Task 6 */

#define RUNS 10

void sum_inner_k(long long expected){
    long long sum = 0;
    double best = 1e30;
    for (int r = 0; r < RUNS; r++) {
        clock_t start = clock();
        sum = 0;
        for (int i = 0; i < W; i++){
            for (int j = 0; j < H; j++){
                for (int k = 0; k < D; k++) {
                    sum += lattice[idx(i, j, k)];
                }
            }
        }
        double t = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (t < best) best = t;
    }
    printf("sum_inner_k: %lld, best of %d runs: %f s\n", sum, RUNS, best);
    assert(sum == expected);
}

void sum_inner_j(long long expected){
    long long sum = 0;
    double best = 1e30;
    for (int r = 0; r < RUNS; r++) {
        clock_t start = clock();
        sum = 0;
        for (int i = 0; i < W; i++){
            for (int k = 0; k < D; k++){
                for (int j = 0; j < H; j++) {
                    sum += lattice[idx(i, j, k)];
                }
            }
        }
        double t = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (t < best) best = t;
    }
    printf("sum_inner_j: %lld, best of %d runs: %f s\n", sum, RUNS, best);
    assert(sum == expected);
}

void sum_inner_i(long long expected){
    long long sum = 0;
    double best = 1e30;
    for (int r = 0; r < RUNS; r++) {
        clock_t start = clock();
        sum = 0;
        for (int k = 0; k < D; k++){
            for (int j = 0; j < H; j++){
                for (int i = 0; i < W; i++) {
                    sum += lattice[idx(i, j, k)];
                }
            }
        }
        double t = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (t < best) best = t;
    }
    printf("sum_inner_i: %lld, best of %d runs: %f s\n", sum, RUNS, best);
    assert(sum == expected);
}

/*verification*/

void test_wrap(void){
    assert(wrap(0, 5) == 0);
    assert(wrap(7, 5) == 2);
    assert(wrap(-1, 5) == 4);
    assert(wrap(-5, 5) == 0);
    assert(wrap(-6, 5) == 4);
}

void test_idx_coords(void){
    for (int n = 0; n < N; n++) {
        Coord c = coords(n);
        assert(idx(c.i, c.j, c.k) == n);
        assert(lattice[idx(c.i, c.j, c.k)] == n);
    }
 
    for (int i = 0; i < W; i++){
        for (int j = 0; j < H; j++){
            for (int k = 0; k < D; k++) {
                int n = idx(i, j, k);
                assert(n >= 0 && n < N);
                Coord c = coords(n);
                assert(c.i == i && c.j == j && c.k == k);
            }
        }
    }
    printf("test_idx_coords passed\n");
}

void test_embedding(void){
    for (int i = 0; i < W; i++){
        for (int j = 0; j < H; j++){
            for (int k = 0; k < D; k++) {
                Phys p = to_phys(i, j, k);
                assert(p.z >= 0 && p.z < 2 * D);
                assert(((p.x + p.y + p.z) & 1) == 0);
                Coord a = to_array(p.x, p.y, p.z);
                assert(a.i == i && a.j == j && a.k == k);
            }
        }
    }
 
    for (int x = 0; x < W; x++){
        for (int y = 0; y < H; y++){
            for (int z = 0; z < 2 * D; z++) {
                if ((x + y + z) & 1) continue;
                Coord a = to_array(x, y, z);
                Phys p = to_phys(a.i, a.j, a.k);
                assert(p.x == x && p.y == y && p.z == z);
            }
        }
    }
    printf("test_embedding passed\n");
}

void test_neighbor(void){
    for (int n = 0; n < N; n++) {
        Coord a = coords(n);
        Phys p = to_phys(a.i, a.j, a.k);
        int nb[12];
 
        for (int d = 0; d < 12; d++) {
            int m = neighbor(n, d);
            assert(m >= 0 && m < N);
            assert(m != n);
            nb[d] = m;
 
            Coord am = coords(m);
            Phys pm = to_phys(am.i, am.j, am.k);
            assert(wrap(pm.x - p.x, W) == wrap(DIR[d][0], W));
            assert(wrap(pm.y - p.y, H) == wrap(DIR[d][1], H));
            assert(wrap(pm.z - p.z, 2 * D) == wrap(DIR[d][2], 2 * D));
            assert(neighbor(m, OPP[d]) == n);
        }
    }
    printf("test_neighbor passed\n");
}

int main(void){
    if (lattice_init() != 0) return 1;

    long long expected = (long long)N * (N - 1) / 2;
    long long sum = sum_kernel();
    assert(sum == expected);

    test_wrap();
    test_idx_coords();
    test_embedding();
    test_neighbor();

    sum_inner_k(expected);
    sum_inner_j(expected);
    sum_inner_i(expected);

    lattice_free();
    return 0;
}