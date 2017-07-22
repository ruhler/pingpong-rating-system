
#include <assert.h>   // for assert
#include <stdio.h>    // for printf

static int Min(int a, int b) {
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

int main() {
  TestSimple();
  TestMultiPath();
  TestCapped();
  printf("passed\n");
  return 0;
}

