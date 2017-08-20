
#include <assert.h>   // for assert
#include <math.h>     // for fabs
#include <stdio.h>    // for printf
#include <stdlib.h>   // for malloc

#include "rate.h"

// Edge --
//   r - The rate of games the src player has won over the dst player.
//   c - The weight of the edge.
typedef struct {
  double r;
  double c;
} Edge;

#define MAX_ERR_THRESHOLD 0.0001

// RatePlayer --
//   Determine the rating of the player with id p.
//
// Inputs:
//   n - The total number of players.
//   edges - The win percentages and edge weights of all players
//   p - The player to compute the rating for.
static double RatePlayer(size_t n, Edge** edges, size_t p)
{
  double rs[n];
  for (size_t k = 0; k < n; ++k) {
    rs[k] = 0;
  }
  rs[p] = 0.5;

  double max_err = 1.0;
  for (size_t i = 0; i < 1000 && max_err > MAX_ERR_THRESHOLD; ++i) {
    max_err = 0;
    for (size_t j = 0; j < n; ++j) {
      if (j != p) {
        double rj = edges[p][j].c * edges[p][j].r;
        for (size_t k = 0; k < n; ++k) {
          if (k != j && k != p) {
            double x = rs[k];
            double y = edges[k][j].r;
            assert(x * y + (1 - x) * (1 - y) > 0);
            rj += edges[k][j].c * (x * y / (x * y + (1 - x) * (1 - y)));
          }
        }
        max_err = fmax(max_err, fabs(rj - rs[j]));
        rs[j] = rj;
      }
    }
  }

  if (max_err > MAX_ERR_THRESHOLD) {
    printf("(Fixed rate error threshold not achieved)\n");
  }

  double rating = 0;
  for (size_t i = 0; i < n; ++i) {
    rating += rs[i];
  }
  return rating / n;
}

void FixedRate(Data* data, double ratings[])
{
  // Convert wins and losses to win rates and edge weights.
  Edge* edges[data->n];
  for (size_t j = 0; j < data->n; ++j) {
    edges[j] = malloc(data->n * sizeof(Edge));
  }

  for (size_t j = 0; j < data->n; ++j) {
    size_t total_games = 0;
    for (size_t k = 0; k < data->n; ++k) {
      total_games += data->wins[k][j] + data->wins[j][k];
    }

    for (size_t k = 0; k < data->n; ++k) {
      size_t games = data->wins[k][j] + data->wins[j][k];
      if (games == 0) {
        edges[k][j].r = 0.5;
        edges[k][j].c = 0;
      } else {
        assert(total_games != 0);
        double w = (double)data->wins[k][j];
        double g = (double)games;
        edges[k][j].r = (w + 0.5) / (g + 1);
        edges[k][j].c = g / (double)total_games;
      }
    }
  }
  
  for (size_t p = 0; p < data->n; ++p) {
    ratings[p] = RatePlayer(data->n, edges, p);
  }
}
