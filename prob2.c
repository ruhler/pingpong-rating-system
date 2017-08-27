
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>

#include "rate.h"

// We assume the following:
//  - The probability a player with rating x beats a player with rating y is 
//    pwin(x, y) = 1 / (1 + exp(y - x)).
//  Note that the player ratings can be translated without affecting the
//  probabilities.
//
//  - Player ratings are normally distributed with some variance S2 and mean
//  0. It is fine to take the mean to be 0 because translating the ratings
//  does not affect the probabilities.
//    prate(x) = exp(-x^2/(2*S2))/sqrt(2*pi*S2)
//  The parameter S2 is unknown, but will turn out not to matter less and less
//  as players play more games against each other.
//
// The ratings are determined by maximizing the probability of seeing the
// match results we saw. That is, maximizing:
//    (Product of prate(x_i) for all players i)
//  * (Product of pwin(x_i, x_j) for all matches where player i beat player j)
//
// This turns out to be equivalent to minimizing the sum of the squares of
// error E(x_i) for each player i, where
//   E(x_i) = w_i - (x/(m_i*S2) + sum(g(i,j) * pwin(x_i, x_j)))
// Here 'm_i' is the total number of matches player i has played.
//      'w_i' is the fraction out of m_i of matches player i has won.
//      'g(i,j)' is the fraction out of m_i of matches player i has played
//               against player j.
// In other words, the error is the difference in the observed win percentage
// and the expected win percentage. The expected win percentage is the average
// of wins against other players, weighted by number of matches played against
// each other player. And there is an additional adjustment of x/(m_i * S2),
// which goes to zero as the player players more matches, representing the
// matches the player won by fluke, as it were.

#define S2 1.0

typedef struct {
  size_t n;     // the total number of players.
  double* m;    // m_i for each player i.
  double* w;    // w_i for each player i.
  double** g;   // g(i, j) for all pairs of players (i, j).
} Params;

// sse_f - the sum of the squared error.
static double sse_f(const gsl_vector* v, void* params) {
    Params* data = (Params*)params;
    double f = 0;
    for (size_t i = 0; i < data->n; i++) {
      double x_i = gsl_vector_get(v, i);
      double w = data->w[i];
      double u = x_i / (data->m[i] * S2);
      for (size_t j = 0; j < data->n; j++) {
        double x_j = gsl_vector_get(v, j);
        u += data->g[i][j] / (1.0 + exp(x_j - x_i));
      }
      f += (w - u) * (w - u);
    }
    return f;
}

// sse_df - the gradiant of the sum of the squared error.
static void sse_df(const gsl_vector* v, void* params, gsl_vector* df) {
  Params* data = (Params*)params;
  for (size_t i = 0; i < data->n; i++) {
    double x_i = gsl_vector_get(v, i);
    double w = data->w[i];
    double u = x_i / (data->m[i] * S2);
    double du = 1.0 / (data->m[i] * S2);
    for (size_t j = 0; j < data->n; j++) {
      double x_j = gsl_vector_get(v, j);
      double k = exp(x_j - x_i);
      double kp1 = 1.0 + k;
      u += data->g[i][j] / kp1;
      du += data->g[i][j] * k / (kp1 * kp1);
    }
    *(gsl_vector_ptr(df, i)) = -2.0 * (w - u) * du;
  }
}

// sse_fdf - both sse_f and sse_df together
static void sse_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
  Params* data = (Params*)params;
  *f = 0;
  for (size_t i = 0; i < data->n; i++) {
    double x_i = gsl_vector_get(v, i);
    double w = data->w[i];
    double u = x_i / (data->m[i] * S2);
    double du = 1.0 / (data->m[i] * S2);
    for (size_t j = 0; j < data->n; j++) {
      double x_j = gsl_vector_get(v, j);
      double k = exp(x_j - x_i);
      double kp1 = 1.0 + k;
      u += data->g[i][j] / kp1;
      du += data->g[i][j] * k / (kp1 * kp1);
    }
    *f += (w - u) * (w - u);
    *(gsl_vector_ptr(df, i)) = -2.0 * (w - u) * du;
  }
}


// Prob2Rate --
//   See documentation in rate.h.
void Prob2Rate(Data* data, double ratings[])
{
  // Compute the parameters.
  Params p;
  p.n = data->n;
  p.m = malloc(data->n * sizeof(double));
  p.w = malloc(data->n * sizeof(double));
  p.g = malloc(data->n * sizeof(double*));
  for (size_t i = 0; i < data->n; ++i) {
    p.g[i] = malloc(data->n * sizeof(double));
    size_t total = 0;
    size_t wins = 0;
    for (size_t j = 0; j < data->n; ++j) {
      total += data->wins[i][j] + data->wins[j][i];
      wins += data->wins[i][j];
    }
    p.m[i] = (double)total;
    p.w[i] = (double)wins / p.m[i];
    for (size_t j = 0; j < data->n; ++j) {
      p.g[i][j] = (double)(data->wins[i][j] + data->wins[j][i]) / p.m[i];
    }
  }

  // Compute the ratings.
  gsl_multimin_function_fdf sse;
  sse.n = data->n;
  sse.f = &sse_f;
  sse.df = &sse_df;
  sse.fdf = &sse_fdf;
  sse.params = (void*)&p;

  gsl_vector* x = gsl_vector_calloc(data->n);
  const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, data->n);
  gsl_multimin_fdfminimizer_set(s, &sse, x, 0.01, 1e-5);

  int status = GSL_CONTINUE;
  for (int i = 0; status == GSL_CONTINUE && i < 1000; i++) {
    if (gsl_multimin_fdfminimizer_iterate(s)) {
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, 1e-3);
  }

  if (status != GSL_SUCCESS) {
    printf("(Prob2Flow failed to find minimum)\n");
  }

  for (int i = 0; i < data->n; i++) {
    ratings[i] = 1.0 / (1.0 + exp(-gsl_vector_get(s->x, i)));
    free(p.g[i]);
  }
  free(p.g);
  free(p.w);
  free(p.m);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
}
