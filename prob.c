
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>

#include "rate.h"

typedef double RATING;

// If player 'a' beat player 'b', we say player 'a' has a K < 1
// probability of beating player 'b', because in practice a player can't
// have a 100% probability of beating some other player. This ensures the
// existence of some minimum error in the RATING choices.
#define K 0.9545

// Estimated probability of a player with rating 'rwinner' beating a
// player with rating 'rloser'.
static double pwin(RATING rwinner, RATING rloser) {
    return 1 / (1 + exp(rloser - rwinner));
}

// Compute the sum of the squared error in estimated probability for each
// match result based on the given player ratings.
static double my_f(const gsl_vector* v, void* params) {
    Data* data = (Data*)params;
    double f = 0;
    for (int i = 0; i < data->n; i++) {
      for (int j = 0; j < data->n; j++) {
        for (int k = 0; k < data->wins[i][j]; k++) {
          RATING rwinner = gsl_vector_get(v, i);
          RATING rloser = gsl_vector_get(v, j);
          double p = pwin(rwinner, rloser);
          double e = K - p;
          f += (e * e);
        }
      }
    }
    return f;
}

// The gradient of f.
static void my_df(const gsl_vector* v, void* params, gsl_vector* df) {
    Data* data = (Data*)params;
    gsl_vector_set_zero(df);
    for (int i = 0; i < data->n; i++) {
      for (int j = 0; j < data->n; j++) {
        for (int k = 0; k < data->wins[i][j]; k++) {
          RATING rwinner = gsl_vector_get(v, i);
          RATING rloser = gsl_vector_get(v, j);
          double p = pwin(rwinner, rloser);
          double e = K - p;
          double half_gradient = e * p * (p - 1);
          *(gsl_vector_ptr(df, i)) += half_gradient;
          *(gsl_vector_ptr(df, j)) -= half_gradient;
        }
      }
    }
}

// Compute both f and df together
static void my_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
    Data* data = (Data*)params;
    *f = 0;
    gsl_vector_set_zero(df);
    for (int i = 0; i < data->n; i++) {
      for (int j = 0; j < data->n; j++) {
        for (int k = 0; k < data->wins[i][j]; k++) {
          RATING rwinner = gsl_vector_get(v, i);
          RATING rloser = gsl_vector_get(v, j);
          double p = pwin(rwinner, rloser);
          double e = K - p;
          *f += (e * e);
          double half_gradient = e * p * (p - 1);
          *(gsl_vector_ptr(df, i)) += half_gradient;
          *(gsl_vector_ptr(df, j)) -= half_gradient;
        }
      }
    }
}


// ProbRate --
//   See documentation in rate.h.
void ProbRate(Data* data, double ratings[])
{
    // Compute the ratings.
    gsl_multimin_function_fdf my_func;
    my_func.n = data->n;
    my_func.f = &my_f;
    my_func.df = &my_df;
    my_func.fdf = &my_fdf;
    my_func.params = (void*)data;

    gsl_vector* x = gsl_vector_calloc(data->n);

    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr;
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, data->n);
    gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-4);

    int status = GSL_CONTINUE;
    for (int i = 0; status == GSL_CONTINUE && i < 1000; i++) {
        if (gsl_multimin_fdfminimizer_iterate(s)) {
            break;
        }
        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
    }

    if (status != GSL_SUCCESS) {
        printf("(ProbFlow failed to find minimum)\n");
    }

    // Scale the ratings so that the mean is 0.5 and the stddev is .125.
    double mean = gsl_stats_mean(s->x->data, s->x->stride, s->x->size);
    double stddev = gsl_stats_sd_m(s->x->data, s->x->stride, s->x->size, mean);
    if (stddev == 0.0) {
        stddev = 0.125;
    }
    gsl_vector_add_constant(s->x, -mean);
    gsl_vector_scale(s->x, 0.125/stddev);
    gsl_vector_add_constant(s->x, 0.5);

    for (int i = 0; i < data->n; i++) {
      ratings[i] = gsl_vector_get(s->x, i);
    }

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
}
