
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>

typedef int PLAYER_ID;
typedef double RATING;

typedef struct {
    PLAYER_ID winner;
    PLAYER_ID loser;
} MatchResult;

typedef struct {
    int playerc;
    const char** playerv;
    int gamec;
    MatchResult* gamev;
} MatchResults;

// If player 'a' beat player 'b', we say player 'a' has a K < 1
// probability of beating player 'b', because in practice a player can't
// have a 100% probability of beating some other player. This ensures the
// existence of some minimum error in the RATING choices.
#define K 0.9545

// Estimated probability of a player with rating 'rwinner' beating a
// player with rating 'rloser'.
double pwin(RATING rwinner, RATING rloser) {
    return 1 / (1 + exp(rloser - rwinner));
}

// Compute the sum of the squared error in estimated probability for each
// match result based on the given player ratings.
double my_f(const gsl_vector* v, void* params) {
    MatchResults* matches = (MatchResults*)params;
    double f = 0;
    for (int i = 0; i < matches->gamec; i++) {
        PLAYER_ID winner = matches->gamev[i].winner;
        PLAYER_ID loser = matches->gamev[i].loser;
        RATING rwinner = gsl_vector_get(v, winner);
        RATING rloser = gsl_vector_get(v, loser);
        double p = pwin(rwinner, rloser);
        double e = K - p;
        f += (e * e);
    }
    return f;
}

// The gradient of f.
void my_df(const gsl_vector* v, void* params, gsl_vector* df) {
    MatchResults* matches = (MatchResults*)params;
    gsl_vector_set_zero(df);
    for (int i = 0; i < matches->gamec; i++) {
        PLAYER_ID winner = matches->gamev[i].winner;
        PLAYER_ID loser = matches->gamev[i].loser;
        RATING rwinner = gsl_vector_get(v, winner);
        RATING rloser = gsl_vector_get(v, loser);
        double p = pwin(rwinner, rloser);
        double e = K - p;
        double half_gradient = e * p * (p - 1);
        *(gsl_vector_ptr(df, winner)) += half_gradient;
        *(gsl_vector_ptr(df, loser)) -= half_gradient;
    }
}

// Compute both f and df together
void my_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
    MatchResults* matches = (MatchResults*)params;
    *f = 0;
    gsl_vector_set_zero(df);
    for (int i = 0; i < matches->gamec; i++) {
        PLAYER_ID winner = matches->gamev[i].winner;
        PLAYER_ID loser = matches->gamev[i].loser;
        RATING rwinner = gsl_vector_get(v, winner);
        RATING rloser = gsl_vector_get(v, loser);
        double p = pwin(rwinner, rloser);
        double e = K - p;
        *f += (e * e);
        double half_gradient = e * p * (p - 1);
        *(gsl_vector_ptr(df, winner)) += half_gradient;
        *(gsl_vector_ptr(df, loser)) -= half_gradient;
    }
}

int main(int argc, char* argv) {
    const char* players[] = {"richard", "jessie", "x"};
    MatchResult games[] = {
        {0, 1}, {1, 2}, {2, 1},
    };

    MatchResults matches;
    matches.playerv = players;
    matches.playerc = sizeof(players)/sizeof(players[0]);
    matches.gamev = games;
    matches.gamec = sizeof(games)/sizeof(games[0]);

    gsl_multimin_function_fdf my_func;
    my_func.n = matches.playerc;
    my_func.f = &my_f;
    my_func.df = &my_df;
    my_func.fdf = &my_fdf;
    my_func.params = (void*) &matches;

    gsl_vector* x = gsl_vector_calloc(matches.playerc);

    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr;
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, matches.playerc);
    gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-4);

    int status = GSL_CONTINUE;
    for (int i = 0; status == GSL_CONTINUE && i < 1000; i++) {
        if (gsl_multimin_fdfminimizer_iterate(s)) {
            break;
        }
        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
    }

    if (status != GSL_SUCCESS) {
        printf("failed to find minimum\n");
    }

    // Scale the ratings so that the mean is 1000 and the stddev is 250.
    double mean = gsl_stats_mean(s->x->data, s->x->stride, s->x->size);
    double stddev = gsl_stats_sd_m(s->x->data, s->x->stride, s->x->size, mean);
    if (stddev == 0.0) {
        stddev = 250.0;
    }
    gsl_vector_add_constant(s->x, -mean);
    gsl_vector_scale(s->x, 250.0/stddev);
    gsl_vector_add_constant(s->x, 1000.0);

    for (int i = 0; i < matches.playerc; i++) {
        printf("%s %f\n", matches.playerv[i], gsl_vector_get(s->x, i));
    }

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    return 0;
}
