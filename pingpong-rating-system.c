
#include <math.h>     // for exp, sqrt
#include <stdio.h>    // for scanf, printf
#include <stdlib.h>   // for malloc, free
#include <string.h>   // for strcmp

// MatchHistory --
//   Information about players and match history.
//
// Fields:
//   n - The total number of players.
//   players - The names of the players by player id.
//   wins - wins[i][j] is the number of matches player i won against player j.
//   total_wins - total_wins[i] is the total number of wins by player i.
//   total_losses - total_losses[i] is the total number of losses by player i.
typedef struct {
  size_t n;
  char** players;
  size_t** wins;
  size_t* total_wins;
  size_t* total_losses;
} MatchHistory;

// We assume the following:
//  - The probability a player with rating x beats a player with rating y is
//    f(x, y) = 1 / (1 + exp(y - x)).
//  Note that the player ratings can be translated without affecting the
//  probabilities.
//
//  - Player ratings are normally distributed with some variance SIGMA_SQUARED and mean
//  0. It is fine to take the mean to be 0 because translating the ratings
//  does not affect the probabilities.
//    g(x) = exp(-x^2/(2*SIGMA_SQUARED))/sqrt(2*pi*SIGMA_SQUARED)
//  The parameter SIGMA_SQUARED is unknown, but will turn out to matter less and less
//  as players play more games against each other.
//
// The ratings are determined by maximizing the probability of seeing the
// match results we saw. That is, maximizing:
//    (Product of g(x_i) for all players i)
//  * (Product of f(x_i, x_j) for all matches where player i beat player j)
//
// This turns out to be equivalent to minimizing the sum of the squares of
// error E(x_i) for each player i, where
//   E(x_i) = w_i - (x/(m_i*SIGMA_SQUARED) + sum(g(i,j) * f(x_i, x_j)))
// Here 'm_i' is the total number of matches player i has played.
//      'w_i' is the fraction out of m_i of matches player i has won.
//      'g(i,j)' is the fraction out of m_i of matches player i has played
//               against player j.
// In other words, the error is the difference in the observed win percentage
// and the expected win percentage. The expected win percentage is the average
// of wins against other players, weighted by number of matches played against
// each other player. And there is an additional adjustment of x/(m_i *
// SIGMA_SQUARED),
// which goes to zero as the player players more matches, representing the
// matches the player won by fluke, as it were.

#define SIGMA_SQUARED 1.0f

// The target mean and standard deviation of ratings after normalization.
#define NORMAL_MEAN 1000.0
#define NORMAL_STDDEV 250.0

// A small delta value used for precision of prenormalized computations
// leading to about 1 post-normalized rating point.
#define D (1.0 / (NORMAL_STDDEV * 4.0))

// Params
//   Parameters of the sum of squared error function.
typedef struct {
  size_t n;     // the total number of players.
  double* m;    // m_i for each player i.
  double* w;    // w_i for each player i.
  double** g;   // g(i, j) for all pairs of players (i, j).
} Params;

static size_t PlayerId(char* name, MatchHistory* history, size_t* capacity);
static MatchHistory* ReadMatchHistory();
static void FreeMatchHistory(MatchHistory* history);

static void ComputeGradients(double* v, Params* data, double* df);

static void Rate(MatchHistory* history, double ratings[]);

static size_t* Min(size_t* a, size_t* b);
static void SortPlayers(size_t n, double* ratings, size_t* sorted);


// PlayerId --
//   Look up or create a player id for the given player.
//
// Inputs:
//   name - The name of the player to look up the id for. This name should be
//          dynamically allocated using malloc. This function takes ownership
//          of the name, using it or freeing it as appropriate.
//   history - The existing match history.
//   capacity - The number of players that space has been allocated for in
//              data.
//
// Results:
//   The id of the player with the given name.
//
// Side effects:
//   If the player is not found in the data, the player will be added to the
//   data. If there is not enough capacity to add the player, the capacity of
//   the data is expanded first. Takes ownership of the memory for the name,
//   using it or freeing it as appropriate.
static size_t PlayerId(char* name, MatchHistory* history, size_t* capacity)
{
  for (size_t i = 0; i < history->n; ++i) {
    if (strcmp(name, history->players[i]) == 0) {
      free(name);
      return i;
    }
  }

  if (*capacity == history->n) {
    *capacity = (*capacity == 0) ? 1 : 2 * (*capacity);
    history->players = realloc(history->players, *capacity * sizeof(char*));
    history->wins = realloc(history->wins, *capacity * sizeof(size_t*));
    for (size_t i = 0; i < history->n; ++i) {
      history->wins[i] = realloc(history->wins[i], *capacity * sizeof(size_t));
    }
    history->total_wins = realloc(history->total_wins, *capacity * sizeof(size_t));
    history->total_losses = realloc(history->total_losses, *capacity * sizeof(size_t));
  }

  history->players[history->n] = name;
  history->wins[history->n] = malloc(*capacity * sizeof(size_t));
  for (size_t i = 0; i < history->n; ++i) {
    history->wins[i][history->n] = 0;
    history->wins[history->n][i] = 0;
  }
  history->total_wins[history->n] = 0;
  history->total_losses[history->n] = 0;
  return history->n++;
}

// ReadMatchHistory --
//   Read match history from stdin.
//
// Input:
//   Each line of input data is in the form "winner loser [n]", where
//   winner is the name of the winning player, loser is the name of the losing
//   player, and n is a positive integer indicating the number of matches in
//   which the winner beat the loser. If n is not provided it defaults to 1.
//   There may be multiple lines with the same winner loser pairs, in which
//   case the sum of the counts for all lines are used to determine the total
//   number of matches for which the winner beat the loser.
//
// Results:
//   The parsed match data, or null if there was an error parsing the data.
//
// Side effects:
//   Reads input from stdin and allocates a MatchHistory structure using
//   malloc. The caller is responsible for calling FreeMatchHistory on the
//   returned data to reclaim the memory for the data once it is no longer
//   needed.
static MatchHistory* ReadMatchHistory()
{
  MatchHistory* history = malloc(sizeof(MatchHistory));
  history->n = 0;
  history->players = NULL;
  history->wins = NULL;
  history->total_wins = NULL;
  history->total_losses = NULL;

  size_t capacity = 0;
  char* winner;
  char* loser;
  size_t wins = 1;
  while (scanf(" %ms %ms %zu", &winner, &loser, &wins) >= 2) {
    size_t w = PlayerId(winner, history, &capacity);
    size_t l = PlayerId(loser, history, &capacity);
    history->wins[w][l] += wins;
    history->total_wins[w] += wins;
    history->total_losses[l] += wins;
    wins = 1;
  }
  return history;
}

// FreeMatchHistory --
//   Free memory used by the given match history.
//
// Inputs:
//   history - the match history to free.
//
// Results:
//   none.
//
// Side effects:
//   Memory used by the given match history is freed. Future accesses to the
//   match history are undefined.
static void FreeMatchHistory(MatchHistory* history)
{
  for (size_t i = 0; i < history->n; ++i) {
    free(history->wins[i]);
  }
  free(history->players);
  free(history->wins);
  free(history->total_wins);
  free(history->total_losses);
  free(history);
}

// ComputeGradients - The gradiant of the sum of the squared error (SSE).
static void ComputeGradients(double* v, Params* data, double* df) {
  for (size_t i = 0; i < data->n; i++) {
    double x_i = v[i];
    double w = data->w[i];
    double u = x_i / (data->m[i] * SIGMA_SQUARED);
    double du = 1.0 / (data->m[i] * SIGMA_SQUARED);
    for (size_t j = 0; j < data->n; j++) {
      double x_j = v[j];
      double k = exp(x_j - x_i);
      double kp1 = 1.0 + k;
      u += data->g[i][j] / kp1;
      du += data->g[i][j] * k / (kp1 * kp1);
    }
    df[i] = -2.0 * (w - u) * du;
  }
}

// Rate --
//   Rate players.
//
// Inputs:
//   history - The match history.
//   rating - Output array of data->n ratings corresponding to player ratings.
//
// Results:
//   none.
//
// Side effects:
//   Sets rating[i] to the rating of the ith player.
static void Rate(MatchHistory* history, double ratings[])
{
  // Compute the parameters.
  Params p;
  p.n = history->n;
  p.m = malloc(history->n * sizeof(double));
  p.w = malloc(history->n * sizeof(double));
  p.g = malloc(history->n * sizeof(double*));
  for (size_t i = 0; i < history->n; ++i) {
    p.g[i] = malloc(history->n * sizeof(double));
    p.m[i] = (double)(history->total_wins[i] + history->total_losses[i]);
    p.w[i] = (double)history->total_wins[i] / p.m[i];
    for (size_t j = 0; j < history->n; ++j) {
      p.g[i][j] = (double)(history->wins[i][j] + history->wins[j][i]) / p.m[i];
    }
  }

  // Initialize the ratings.
  for (size_t i = 0; i < history->n; ++i) {
    ratings[i] = 0.0;
  }

  double max_gradient = 2 * D;
  double gradients[history->n];

  size_t progress = 0;
  while (max_gradient > D) {
    ComputeGradients(ratings, &p, gradients);

    max_gradient = 0.0;
    for (size_t i = 0; i < history->n; ++i) {
      max_gradient = fmax(max_gradient, gradients[i]);
      ratings[i] -= D * gradients[i];
    }

    size_t nprogress = (size_t)(100.0 * D / max_gradient);
    if (nprogress > progress) {
      progress = nprogress;
      printf("\r%zi%% done", progress);
      fflush(stdout);
    }
  }
  printf("\n");

  for (size_t i = 0; i < history->n; i++) {
    free(p.g[i]);
  }
  free(p.g);
  free(p.w);
  free(p.m);
}

static size_t* Min(size_t* a, size_t* b)
{
  return a < b ? a : b;
}

// SortPlayers --
//   Sort players by rating.
//
// Inputs:
//   n - The number of players.
//   ratings - The rating for each player.
//   sorted - Output array to fill with player ids in order of rating.
//
// Results:
//   None.
//
// Side effects:
//   Sets sorted to be a list of player ids in decreasing order of rating.
static void SortPlayers(size_t n, double* ratings, size_t* sorted)
{
  size_t temp[n];
  size_t* src = sorted;
  size_t* dst = temp;

  for (size_t i = 0; i < n; ++i) {
    src[i] = i;
  }

  for (size_t w = 1; w < n; w *= 2) {
    // Merge subarrays of width w from src to dst.
    for (size_t i = 0; i < n; i += 2*w) {
      size_t* a = src + i;
      size_t* b = a + w;
      size_t* aend = Min(a + w, src + n);
      size_t* bend = Min(b + w, src + n);
      size_t* z = dst + i;
      while (a < aend || b < bend) {
        if (b >= bend || (a < aend && ratings[*a] >= ratings[*b])) {
          *(z++) = *(a++);
        } else {
          *(z++) = *(b++);
        }
      }
    }
    size_t* tmp = src;
    src = dst;
    dst = tmp;
  }

  if (src != sorted) {
    memcpy(sorted, src, n);
  }
}

// main --
//   Main entry point. Reads match history from stdin, computes rating
//   results, sorts players by rating, and outputs the rating results.
//
// Inputs:
//   None.
//
// Returns:
//   zero on success, non-zero on error.
//
// Side effects:
//   Reads match data from stdin and outputs rating results to stdout.
int main()
{
  MatchHistory* history = ReadMatchHistory();
  double ratings[history->n];
  Rate(history, ratings);

  double var = 0;
  for (size_t i = 0; i < history->n; ++i) {
    var += ratings[i] * ratings[i];
  }
  var /= history->n;
  double a = NORMAL_STDDEV / sqrt(var);

  size_t sorted[history->n];
  SortPlayers(history->n, ratings, sorted);

  printf("player raw normal matches wins losses\n");
  for (size_t i = 0; i < history->n; ++i) {
    size_t p = sorted[i];
    double raw = ratings[p];
    double normal = a * raw + NORMAL_MEAN;
    printf("%10s %+1.4f %4.0f %zi %zi %zi\n", history->players[p], raw, normal,
        history->total_wins[p] + history->total_losses[p],
        history->total_wins[p], history->total_losses[p]);
  }

  FreeMatchHistory(history);
  return 0;
}
