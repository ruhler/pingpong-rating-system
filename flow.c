
#include <assert.h>   // for assert
#include <limits.h>   // for INT_MAX
#include <stdlib.h>   // for malloc
#include <stdio.h>    // for printf
#include <string.h>   // for memset

// Min --
//   Return the smaller of two numbers.
//
// Inputs:
//   a - the first number.
//   b - the second number.
//
// Result:
//   The smaller of the two numbers.
//
// Side effects:
//   None.
static int Min(int a, int b)
{
  return a < b ? a : b;
}

// Flow --
//   Compute the maximum flow from src to dst in the graph.
//
// Inputs:
//   n - The number of nodes in the graph.
//   edges - edges[i][j] is the capacity from node i to j in the graph.
//   src - The source node.
//   dst - The destination node.
//   cap - The maximum allowed flow capacity.
//
// Result:
//   The maximum flow from src to dst in the graph.
//
// Side effects:
//  reduces capacity on the edges as necessary to sustain the maximum flow,
//  producing the residual graph.
static int Flow(int n, int** edges, int cap, int src, int dst)
{
  if (src == dst) {
    return cap;
  }

  int total = 0;
  for (int i = 0; cap > 0 && i < n; ++i) {
    int c = edges[src][i];
    edges[src][i] = 0;
    int flow = Flow(n, edges, Min(cap, c), i, dst);
    total += flow;
    cap -= flow;
    edges[src][i] = c - flow;
  }
  return total;
}

// Copy --
//   Copy a graph from src to dst.
//
// Inputs:
//   n - The number of nodes.
//   src - The n*n edges of the source graph.
//   dst - The n*n edges to copy the source graph to.
//
// Results:
//   none.
//
// Side effects:
//   dst[i][j] = src[i][j] for all i, j less than n.
static void Copy(int n, int** src, int** dst)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dst[i][j] = src[i][j];
    }
  }
}

// Alloc --
//   Allocate space for the edges of a graph with n players.
//
// Inputs:
//   n - The number of nodes in the graph.
//
// Result:
//   Newly allocated n*n edges initialized to 0.
//
// Side effects:
//   Allocates memory that should be freed using the Free function when no
//   longer needed.
static int** Alloc(int n)
{
  int** edges = malloc(n * sizeof(int*));
  for (int i = 0; i < n; ++i) {
    edges[i] = malloc(n * sizeof(int));
    memset(edges[i], 0, sizeof(int));
  }
  return edges;
}

// Free --
//   Free memory for edges allocated using Alloc.
//
// Inputs:
//   n - The number of of nodes in the graph.
//   edges - The edges to free.
//
// Result:
//   None.
//
// Side effects:
//   Frees memory for the given edges of the graph.
static void Free(int n, int** edges)
{
  for (int i = 0; i < n; ++i) {
    free(edges[i]);
  }
  free(edges);
}

// FlowAll --
//   Compute the max flow between all pairs of nodes in the given graph.
//
// Inputs:
//   n - The number of nodes in the graph.
//   edges - edges[i][j] is the capacity from node i to j in the graph.
//
// Result:
//   For each pair of nodes (i, j), the max flow from i to j in the graph, or
//   0 if i == j.
//
// Side effects:
//   Allocates an edges structure using Alloc that must be free with Free when
//   no longer needed.
int** FlowAll(int n, int** edges)
{
  int** flows = Alloc(n);
  int** copy = Alloc(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        Copy(n, edges, copy);
        flows[i][j] = Flow(n, copy, INT_MAX, i, j);
      }
    }
  }
  Free(n, copy);
  return flows;
}

// Rate --
//   Rate a player based on the given flows.
//
// Inputs:
//   n - The number of players.
//   flows - The pre-computed max flows between all players.
//   p - The player to rate.
//
// Results:
//   The players rating, as a number between 0 and 1.
//
// Side effects:
//   None.
static double Rate(int n, int** flows, int p)
{
  double rating = 0;
  for (int i = 0; i < n; ++i) {
    double wins = (double)flows[p][i];
    double losses = (double)flows[i][p];
    rating += wins / (wins + losses + 1.0);
  }
  return rating / (double)n;
}

static void TestSimple()
{
  // a -- 2 --> b
  int a[2] = {0, 2};
  int b[2] = {0, 0};
  int* edges[2] = {a, b};
  assert(Flow(2, edges, 5, 0, 1) == 2);
}

static void TestMultiPath()
{
  // a -- 1 --> b -- 1 --> c
  //  \------->----- 1 -->/
  int a[3] = {0, 2, 1};
  int b[3] = {0, 0, 1};
  int c[3] = {0, 0, 0};
  int* edges[3] = {a, b, c};
  assert(Flow(3, edges, 5, 0, 2) == 2);
}

static void TestCapped()
{
  // a -- 1 --> b -- 1 --> c -- 1 --> d
  //             \------->----- 1 -->/
  int a[4] = {0, 1, 0, 0};
  int b[4] = {0, 0, 1, 1};
  int c[4] = {0, 0, 0, 1};
  int d[4] = {0, 0, 0, 0};
  int* edges[4] = {a, b, c, d};
  assert(Flow(4, edges, 10, 0, 3) == 1);
}

static void TestFlowAll()
{
  // a -- 1 --> b -- 1 --> c -- 1 --> d
  //             \------->----- 1 -->/
  int a[4] = {0, 1, 0, 0};
  int b[4] = {0, 0, 1, 1};
  int c[4] = {0, 0, 0, 1};
  int d[4] = {0, 0, 0, 0};
  int* edges[4] = {a, b, c, d};

  int aw[4] = {0, 1, 1, 1};
  int bw[4] = {0, 0, 1, 2};
  int cw[4] = {0, 0, 0, 1};
  int dw[4] = {0, 0, 0, 0};
  int* flows_want[4] = {aw, bw, cw, dw};

  int** flows_got = FlowAll(4, edges);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      assert(flows_want[i][j] == flows_got[i][j]);
    }
  }
  Free(4, flows_got);
}

static void TestRate()
{
  int aw[4] = {0, 1, 1, 1};
  int bw[4] = {0, 0, 1, 3};
  int cw[4] = {0, 0, 0, 1};
  int dw[4] = {0, 0, 0, 0};
  int* flows[4] = {aw, bw, cw, dw};
  assert(Rate(4, flows, 0) == 0.375);
  assert(Rate(4, flows, 1) == 0.3125);
  assert(Rate(4, flows, 2) == 0.125);
  assert(Rate(4, flows, 3) == 0.0);
}

int main() {
  TestSimple();
  TestMultiPath();
  TestCapped();
  TestFlowAll();
  TestRate();
  printf("passed\n");
  return 0;
}

