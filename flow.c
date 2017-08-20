
#include <assert.h>   // for assert
#include <limits.h>   // for INT_MAX
#include <stdlib.h>   // for malloc
#include <stdio.h>    // for printf
#include <string.h>   // for memset

#include "rate.h"

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
static size_t Min(size_t a, size_t b)
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
static size_t Flow(size_t n, size_t** edges, size_t cap, size_t src, size_t dst)
{
  if (src == dst) {
    return cap;
  }

  size_t total = 0;
  for (size_t i = 0; cap > 0 && i < n; ++i) {
    size_t c = edges[src][i];
    edges[src][i] = 0;
    size_t flow = Flow(n, edges, Min(cap, c), i, dst);
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
static void Copy(size_t n, size_t** src, size_t** dst)
{
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
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
static size_t** Alloc(size_t n)
{
  size_t** edges = malloc(n * sizeof(size_t*));
  for (size_t i = 0; i < n; ++i) {
    edges[i] = malloc(n * sizeof(size_t));
    memset(edges[i], 0, sizeof(size_t));
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
static void Free(size_t n, size_t** edges)
{
  for (size_t i = 0; i < n; ++i) {
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
size_t** FlowAll(size_t n, size_t** edges)
{
  size_t** flows = Alloc(n);
  size_t** copy = Alloc(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
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
static double Rate(size_t n, size_t** flows, size_t p)
{
  double rating = 0;
  for (size_t i = 0; i < n; ++i) {
    double wins = (double)flows[p][i];
    double losses = (double)flows[i][p];
    rating += wins / (wins + losses + 1.0);
  }
  return rating / (double)n;
}

static void TestSimple()
{
  // a -- 2 --> b
  size_t a[2] = {0, 2};
  size_t b[2] = {0, 0};
  size_t* edges[2] = {a, b};
  assert(Flow(2, edges, 5, 0, 1) == 2);
}

static void TestMultiPath()
{
  // a -- 1 --> b -- 1 --> c
  //  \------->----- 1 -->/
  size_t a[3] = {0, 2, 1};
  size_t b[3] = {0, 0, 1};
  size_t c[3] = {0, 0, 0};
  size_t* edges[3] = {a, b, c};
  assert(Flow(3, edges, 5, 0, 2) == 2);
}

static void TestCapped()
{
  // a -- 1 --> b -- 1 --> c -- 1 --> d
  //             \------->----- 1 -->/
  size_t a[4] = {0, 1, 0, 0};
  size_t b[4] = {0, 0, 1, 1};
  size_t c[4] = {0, 0, 0, 1};
  size_t d[4] = {0, 0, 0, 0};
  size_t* edges[4] = {a, b, c, d};
  assert(Flow(4, edges, 10, 0, 3) == 1);
}

static void TestFlowAll()
{
  // a -- 1 --> b -- 1 --> c -- 1 --> d
  //             \------->----- 1 -->/
  size_t a[4] = {0, 1, 0, 0};
  size_t b[4] = {0, 0, 1, 1};
  size_t c[4] = {0, 0, 0, 1};
  size_t d[4] = {0, 0, 0, 0};
  size_t* edges[4] = {a, b, c, d};

  size_t aw[4] = {0, 1, 1, 1};
  size_t bw[4] = {0, 0, 1, 2};
  size_t cw[4] = {0, 0, 0, 1};
  size_t dw[4] = {0, 0, 0, 0};
  size_t* flows_want[4] = {aw, bw, cw, dw};

  size_t** flows_got = FlowAll(4, edges);
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      assert(flows_want[i][j] == flows_got[i][j]);
    }
  }
  Free(4, flows_got);
}

static void TestRate()
{
  size_t aw[4] = {0, 1, 1, 1};
  size_t bw[4] = {0, 0, 1, 3};
  size_t cw[4] = {0, 0, 0, 1};
  size_t dw[4] = {0, 0, 0, 0};
  size_t* flows[4] = {aw, bw, cw, dw};
  assert(Rate(4, flows, 0) == 0.375);
  assert(Rate(4, flows, 1) == 0.3125);
  assert(Rate(4, flows, 2) == 0.125);
  assert(Rate(4, flows, 3) == 0.0);
}

static void TestAll()
{
  TestSimple();
  TestMultiPath();
  TestCapped();
  TestFlowAll();
  TestRate();
  printf("(Flow Tests Passed)\n");
}

// FlowRate -- see documentation in rate.h
void FlowRate(Data* data, double ratings[])
{
  // Always run the tests for now.
  TestAll();

  size_t** flows = FlowAll(data->n, data->wins);
  for (size_t i = 0; i < data->n; ++i) {
    ratings[i] = Rate(data->n, flows, i);
  }
  Free(data->n, flows);
}

