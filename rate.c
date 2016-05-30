
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef int PLAYER_ID;

typedef struct {
    PLAYER_ID winner;
    PLAYER_ID loser;
} MatchResult;

double logistic(double x) {
    return 1 / (1 + exp(-x));
}

double likelihood(int matchc, MatchResult* matchv, double* ratingv) {
    double x = 1.0;
    for (int i = 0; i < matchc; i++) {
        x *= logistic(ratingv[matchv[i].winner] - ratingv[matchv[i].loser]);
    }
    return x;
}

// Prerequistit: ratein[x] = rateout[x] for all x.
void refine(int matchc, MatchResult* matchv, double* ratein, double* rateout) {
    for (int i = 0; i < matchc; i++) {
        double exchange = logistic(ratein[matchv[i].loser] - ratein[matchv[i].winner]);
        rateout[matchv[i].winner] += exchange;
        rateout[matchv[i].loser] -= exchange;
    }
}

int main(int argc, char* argv) {
    const char* players[] = {"loser", "winner", "richard", "jessie"};
    double MIN_RATING = 0.0;
    double MAX_RATING = 1.0;
    double ratings0[] = {MIN_RATING, MAX_RATING, 0.5, 0.5};
    double ratings1[] = {MIN_RATING, MAX_RATING, 0.5, 0.5};
    double* ratings[2] = {ratings0, ratings1};
    MatchResult matches[] = {{2, 0}, {3, 0}, {1, 2}, {1, 3}, {2, 3}};
    int num_players = sizeof(players)/sizeof(players[0]);
    int num_matches = sizeof(matches)/sizeof(matches[0]);

    double like = 0.0;
    double nlike = likelihood(num_matches, matches, ratings[0]);
    for (int i = 0; nlike > like ; i++) {
        like = nlike;
        double* old_ratings = ratings[i%2];
        double* new_ratings = ratings[(i+1)%2];

        for (int j = 0; j < num_players; j++) {
            printf("%f ", old_ratings[j]);
            new_ratings[j] = old_ratings[j];
        }
        printf("\n");

        refine(num_matches, matches, old_ratings, new_ratings);
        new_ratings[0] = MIN_RATING;
        new_ratings[1] = MAX_RATING;
        nlike = likelihood(num_matches, matches, new_ratings);
    }

    return 0;
}
