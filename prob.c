
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
    int playera;                // allocated number of players in playerv.
    int playerc;                // count of the number of players in playerv.
    char** playerv;
    int gamea;                  // allocated number of games in gamev.
    int gamec;                  // count of number of games in gamev.
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

// Returns the player id associated with the given player. Adds the player to
// the players list if it is not already present.
//
// The name will be either stored in the playerv vector, or freed by this
// function.
int player_id(MatchResults* matches, char* name) {
    for (int i = 0; i < matches->playerc; i++) {
        if (strcmp(matches->playerv[i], name) == 0) {
            free(name);
            return i;
        }
    }

    if (matches->playerc == matches->playera) {
        matches->playera *= 2;
        matches->playerv = realloc(matches->playerv, matches->playera);
    }

    int id = matches->playerc;
    matches->playerv[id] = name;
    matches->playerc++;
    return id;
}

int main(int argc, char* argv) {
    MatchResults matches;
    matches.playera = 1;
    matches.playerc = 0;
    matches.playerv = calloc(matches.playera, sizeof(char*));
    matches.gamea = 1;
    matches.gamec = 0;
    matches.gamev = malloc(matches.gamea * sizeof(MatchResult));

    // Read match data from stdin.
    char* winner;
    char* loser;
    while (scanf(" %ms %ms", &winner, &loser) == 2) {
        if (matches.gamec == matches.gamea) {
            matches.gamea *= 2;
            matches.gamev = realloc(matches.gamev, matches.gamea);
        }

        matches.gamev[matches.gamec].winner = player_id(&matches, winner);
        matches.gamev[matches.gamec].loser = player_id(&matches, loser);
        matches.gamec++;
    }

    // Compute the ratings.
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

    for (int i = 0; i < matches.playerc; i++) {
        free(matches.playerv[i]);
    }
    free(matches.playerv);
    free(matches.gamev);
    return 0;
}
