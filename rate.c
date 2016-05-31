
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef int PLAYER_ID;

typedef struct {
    PLAYER_ID winner;
    PLAYER_ID loser;
} MatchResult;

double f(double x) {
    // return fmax(0.0, fmin(1.0, x + 0.5));
    return 1 / (1 + exp(-x));
}

double mse(int matchc, MatchResult* matchv, double* ratingv) {
    double sum = 0;
    for (int i = 0; i < matchc; i++) {
        double e = 1.0 - f(ratingv[matchv[i].winner] - ratingv[matchv[i].loser]);
        sum += e*e;
    }
    return sum;
}

int main(int argc, char* argv) {
    const char* players[] = {"richard", "jessie", "x"};
    double ratings[] = {0.5, 0.5, 0.5, 0.5, 0.5};
    double nratings[] = {0.5, 0.5, 0.5, 0.5, 0.5};
    MatchResult matches[] = {
        {0, 1}, {0, 1},
        {0, 1}, {0, 1},
        {0, 1}, {0, 1}, {1, 2},
    };
    int num_players = sizeof(players)/sizeof(players[0]);
    int num_matches = sizeof(matches)/sizeof(matches[0]);

    double err = mse(num_matches, matches, ratings);
    double nerr = err;
    do {
        // Commit the new ratings.
        for (int j = 0; j < num_players; j++) {
            ratings[j] = nratings[j];
            printf("%f ", ratings[j]);
        }
        err = nerr;
        printf(" ERR %f\n", err);

        // Refine the ratings.
        for (int j = 0; j < num_matches; j++) {
            double exchange = f(ratings[matches[j].loser] - ratings[matches[j].winner]);
            nratings[matches[j].winner] += exchange;
            nratings[matches[j].loser] -= exchange;
        }
        for (int j = 0; j < num_players; j++) {
            nratings[j] *= 0.99;
        }
        nerr = mse(num_matches, matches, nratings);
    } while (nerr < err);
    return 0;
}
