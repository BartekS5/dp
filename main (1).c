#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define W (32)
#define H (32)
#define D (32)
#define N (W * H * D)
#define NDIR (12)

int *ID;

const int dir[NDIR][3] = {{1, 1, 0},  {-1, 1, 0}, {-1, -1, 0}, {1, -1, 0},
                          {1, 0, 1},  {0, 1, 1},  {-1, 0, 1},  {0, -1, 1},
                          {1, 0, -1}, {0, 1, -1}, {-1, 0, -1}, {0, -1, -1}};

typedef struct {
  int i;
  int j;
  int k;
} coords_t;

void init() {
  ID = malloc(N * sizeof(int));

  if (ID == NULL) {
    perror("malloc");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < N; i++) {
    ID[i] = i;
  }
}

int mod(int a, int m) {
  assert(m > 0);
  int r = a % m;
  return r < 0 ? r + m : r;
}

int idx(int i, int j, int k) {
  assert(i >= 0 && i < W);
  assert(j >= 0 && j < H);
  assert(k >= 0 && k < D);
  return (i * H + j) * D + k;
}

coords_t coords(int n) {
  assert(n >= 0 && n < N);
  coords_t c;
  c.k = n % D;
  n /= D;
  c.j = n % H;
  n /= H;
  c.i = n;
  return c;
}

coords_t arr_phys(int i, int j, int k) {
  assert(i >= 0 && i < W);
  assert(j >= 0 && j < H);
  assert(k >= 0 && k < D);
  coords_t c;
  c.i = i;
  c.j = j;
  c.k = (i + j) % 2 == 0 ? k * 2 : k * 2 + 1;
  assert((c.i + c.j + c.k) % 2 == 0);
  return c;
}

coords_t phys_arr(int x, int y, int z) {
  assert(x >= 0 && x < W);
  assert(y >= 0 && y < H);
  assert(z >= 0 && z < 2 * D);
  assert((x + y + z) % 2 == 0);
  coords_t c;
  c.i = x;
  c.j = y;
  c.k = z / 2;
  return c;
}

int neighbor(int n, int d) {
  assert(n >= 0 && n < N);
  assert(d >= 0 && d < NDIR);
  coords_t arr = coords(n);
  coords_t phys = arr_phys(arr.i, arr.j, arr.k);
  coords_t neigh =
      phys_arr(mod(phys.i + dir[d][0], W), mod(phys.j + dir[d][1], H),
               mod(phys.k + dir[d][2], 2 * D));
  return idx(neigh.i, neigh.j, neigh.k);
}

volatile int sum_k() {
  assert(ID != NULL);
  int sum = 0;
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < H; j++) {
      for (int k = 0; k < D; k++) {
        sum += ID[idx(i, j, k)];
      }
    }
  }
  return sum;
}

volatile int sum_j() {
  assert(ID != NULL);
  int sum = 0;
  for (int k = 0; k < D; k++) {
    for (int i = 0; i < W; i++) {
      for (int j = 0; j < H; j++) {
        sum += ID[idx(i, j, k)];
      }
    }
  }
  return sum;
}

volatile int sum_i() {
  assert(ID != NULL);
  int sum = 0;
  for (int j = 0; j < H; j++) {
    for (int k = 0; k < D; k++) {
      for (int i = 0; i < W; i++) {
        sum += ID[idx(i, j, k)];
      }
    }
  }
  return sum;
}

int opposite(int d) {
  for (int e = 0; e < NDIR; e++) {
    if (dir[e][0] == -dir[d][0] && dir[e][1] == -dir[d][1] &&
        dir[e][2] == -dir[d][2]) {
      return e;
    }
  }
  assert(0 && "direction has no opposite");
  return -1;
}

void test_mod() {
  assert(mod(0, 5) == 0);
  assert(mod(7, 5) == 2);
  assert(mod(-1, 5) == 4);
  assert(mod(-5, 5) == 0);
  assert(mod(-6, 5) == 4);
}

void test_idx_coords() {
  // n -> (i,j,k) -> n
  for (int n = 0; n < N; n++) {
    coords_t c = coords(n);
    assert(idx(c.i, c.j, c.k) == n);
    assert(ID[idx(c.i, c.j, c.k)] == n);
  }

  // (i,j,k) -> n -> (i,j,k)
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < H; j++) {
      for (int k = 0; k < D; k++) {
        int n = idx(i, j, k);
        assert(n >= 0 && n < N);
        coords_t c = coords(n);
        assert(c.i == i && c.j == j && c.k == k);
      }
    }
  }
}

void test_arr_phys() {
  // (i,j,k) -> (x,y,z) -> (i,j,k), every site is a valid FCC site
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < H; j++) {
      for (int k = 0; k < D; k++) {
        coords_t p = arr_phys(i, j, k);
        assert(p.k >= 0 && p.k < 2 * D);
        assert((p.i + p.j + p.k) % 2 == 0);
        coords_t a = phys_arr(p.i, p.j, p.k);
        assert(a.i == i && a.j == j && a.k == k);
      }
    }
  }

  // (x,y,z) -> (i,j,k) -> (x,y,z) for every valid FCC site
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
      for (int z = 0; z < 2 * D; z++) {
        if ((x + y + z) % 2 != 0) {
          continue;
        }
        coords_t a = phys_arr(x, y, z);
        coords_t p = arr_phys(a.i, a.j, a.k);
        assert(p.i == x && p.j == y && p.k == z);
      }
    }
  }
}

void test_dir() {
  for (int d = 0; d < NDIR; d++) {
    // FCC nearest-neighbor vector: two components +-1, one 0
    int zeros = 0;
    for (int a = 0; a < 3; a++) {
      assert(dir[d][a] >= -1 && dir[d][a] <= 1);
      zeros += dir[d][a] == 0;
    }
    assert(zeros == 1);
    // directions are unique
    for (int e = d + 1; e < NDIR; e++) {
      assert(!(dir[e][0] == dir[d][0] && dir[e][1] == dir[d][1] &&
               dir[e][2] == dir[d][2]));
    }
    int o = opposite(d);
    assert(o != d && opposite(o) == d);
  }
}

void test_neighbor() {
  for (int n = 0; n < N; n++) {
    coords_t a = coords(n);
    coords_t p = arr_phys(a.i, a.j, a.k);
    int nb[NDIR];

    for (int d = 0; d < NDIR; d++) {
      int m = neighbor(n, d);
      assert(m >= 0 && m < N);
      assert(m != n);
      nb[d] = m;

      // physical displacement equals dir[d] (modulo the periodic box)
      coords_t am = coords(m);
      coords_t pm = arr_phys(am.i, am.j, am.k);
      assert(mod(pm.i - p.i, W) == mod(dir[d][0], W));
      assert(mod(pm.j - p.j, H) == mod(dir[d][1], H));
      assert(mod(pm.k - p.k, 2 * D) == mod(dir[d][2], 2 * D));

      // stepping back returns to the start
      assert(neighbor(m, opposite(d)) == n);
    }

    // all 12 neighbors are distinct
    for (int d = 0; d < NDIR; d++) {
      for (int e = d + 1; e < NDIR; e++) {
        assert(nb[d] != nb[e]);
      }
    }
  }

  // explicit boundary cases
  // (0,0,0) + (-1,-1,0) wraps to (W-1,H-1,0)
  assert(neighbor(idx(0, 0, 0), 2) == idx(W - 1, H - 1, 0));
  // (0,0,0) + (0,-1,-1) wraps to (0,H-1,2D-1) -> k = D-1
  assert(neighbor(idx(0, 0, 0), 11) == idx(0, H - 1, D - 1));
  // (W-1,H-1,z=2D-2) + (1,1,0) wraps to (0,0,2D-2) -> k = D-1
  assert(neighbor(idx(W - 1, H - 1, D - 1), 0) == idx(0, 0, D - 1));
  // (0,1,z=2D-1) + (1,0,1) wraps to (1,1,0) -> k = 0
  assert(neighbor(idx(0, 1, D - 1), 4) == idx(1, 1, 0));
}

void test_sums() {
  int expected = (int)((long long)N * (N - 1) / 2);
  assert(sum_i() == expected);
  assert(sum_j() == expected);
  assert(sum_k() == expected);
}

#define WARMUP (10)
#define REPS (100) // timed kernel calls

// results are written here so the compiler cannot drop the kernel calls
volatile int sink;

double now() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void time_sum(const char *name, volatile int (*sum)()) {
  // calling through a volatile pointer keeps the kernel from being inlined
  // and hoisted out of the timing loop
  volatile int (*volatile kernel)() = sum;
  int expected = (int)((long long)N * (N - 1) / 2);

  for (int r = 0; r < WARMUP; r++) {
    sink = kernel();
    assert(sink == expected);
  }

  double total = 0.0;
  double best = 1e30;
  for (int r = 0; r < REPS; r++) {
    double t0 = now();
    sink = kernel();
    double t1 = now();
    assert(sink == expected);
    double t = t1 - t0;
    total += t;
    if (t < best) {
      best = t;
    }
  }

  double mean = total / REPS;
  printf("%s: mean %8.3f us, min %8.3f us, %6.3f ns/elem\n", name, mean * 1e6,
         best * 1e6, best * 1e9 / N);
}

void test_timing() {
  printf("Traversal timing (%dx%dx%d, %d calls):\n", W, H, D, REPS);
  time_sum("k innermost", sum_k);
  time_sum("j innermost", sum_j);
  time_sum("i innermost", sum_i);
}

void test() {
  test_mod();
  test_idx_coords();
  test_arr_phys();
  test_dir();
  test_neighbor();
  test_sums();

  printf("All checks passed (N = %d)\n", N);

  test_timing();
}

int main() {
#ifdef NDEBUG
  fprintf(stderr, "warning: built with NDEBUG, asserts are disabled\n");
#endif
  init();
  test();
  free(ID);
  return 0;
}
