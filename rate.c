
#include <stdio.h>    // for scanf, printf
#include <stdlib.h>   // for malloc, free
#include <string.h>   // for strcmp

#include <gsl/gsl_multimin.h> // for gsl_*

// Data --
//   Information about players and match history.
//
// Fields:
//   n - The total number of players.
//   players - The names of the players by player id.
//   wins - wins[i][j] is the number of matches player i won against player j.
typedef struct {
  size_t n;
  char** players;
  size_t** wins;
} Data;

static size_t PlayerId(char* name, Data* data, size_t* capacity);
static Data* ReadData();
static void FreeData(Data* data);

static double sse_f(const gsl_vector* v, void* params);
static void sse_df(const gsl_vector* v, void* params, gsl_vector* df);
static void sse_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df);
static void Rate(Data* data, double ratings[]);


// PlayerId --
//   Look up or create a player id for the given player.
//
// Inputs:
//   name - The name of the player to look up the id for. This name should be
//          dynamically allocated using malloc. This function takes ownership
//          of the name, using it or freeing it as appropriate.
//   data - The existing player data.
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
static size_t PlayerId(char* name, Data* data, size_t* capacity)
{
  for (size_t i = 0; i < data->n; ++i) {
    if (strcmp(name, data->players[i]) == 0) {
      free(name);
      return i;
    }
  }

  if (*capacity == data->n) {
    *capacity = (*capacity == 0) ? 1 : 2 * (*capacity);
    data->players = realloc(data->players, *capacity * sizeof(char*));
    data->wins = realloc(data->wins, *capacity * sizeof(size_t*));
    for (size_t i = 0; i < data->n; ++i) {
      data->wins[i] = realloc(data->wins[i], *capacity * sizeof(size_t));
    }
  }

  data->players[data->n] = name;
  data->wins[data->n] = malloc(*capacity * sizeof(size_t));
  for (size_t i = 0; i < data->n; ++i) {
    data->wins[i][data->n] = 0;
    data->wins[data->n][i] = 0;
  }
  return data->n++;
}

// ReadData --
//   Read player data from stdin.
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
//   Reads input from stdin and allocates a Data structure using malloc. The
//   caller is responsible for calling FreeData on the returned data to
//   reclaim the memory for the data once it is no longer needed.
static Data* ReadData()
{
  Data* data = malloc(sizeof(Data));
  data->n = 0;
  data->players = NULL;
  data->wins = NULL;

  size_t capacity = 0;
  char* winner;
  char* loser;
  size_t wins = 1;
  while (scanf(" %ms %ms %zd", &winner, &loser, &wins) >= 2) {
    size_t w = PlayerId(winner, data, &capacity);
    size_t l = PlayerId(loser, data, &capacity);
    data->wins[w][l] += wins;
    wins = 1;
  }
  return data;
}

// FreeData --
//   Free memory used by the given data.
//
// Inputs:
//   data - the data to free.
//
// Results:
//   none.
//
// Side effects:
//   Memory used by the given data is freed. Future accesses to the data
//   are undefined.
static void FreeData(Data* data)
{
  for (size_t i = 0; i < data->n; ++i) {
    free(data->wins[i]);
  }
  free(data->players);
  free(data->wins);
  free(data);
}

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

// Rate --
//   Rate players.
// 
// Inputs:
//   data - The player data.
//   rating - Output array of data->n ratings corresponding to player ratings.
//
// Results:
//   none.
//
// Side effects:
//   Sets rating[i] to the rating of the ith player.
static void Rate(Data* data, double ratings[])
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
  for (int i = 0; status == GSL_CONTINUE && i < 10000; i++) {
    if (gsl_multimin_fdfminimizer_iterate(s)) {
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, 1e-5);
  }

  if (status != GSL_SUCCESS) {
    fprintf(stderr, "Rate failed to find minimum\n");
    abort();
  }

  for (int i = 0; i < data->n; i++) {
    ratings[i] = gsl_vector_get(s->x, i);
    free(p.g[i]);
  }
  free(p.g);
  free(p.w);
  free(p.m);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
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

int main() {
  Data* data = ReadData();
  double ratings[data->n];
  Rate(data, ratings);

  double var = 0;
  for (size_t i = 0; i < data->n; ++i) {
    var += ratings[i] * ratings[i];
  }
  var /= data->n;
  double a = 250.0 / sqrt(var);

  size_t sorted[data->n];
  SortPlayers(data->n, ratings, sorted);
  for (size_t i = 0; i < data->n; ++i) {
    size_t p = sorted[i];
    double raw = ratings[p];
    double normal = a * raw + 1000.0;
    printf("%10s %1.4f %.0f\n", data->players[p], raw, normal);
  }

  FreeData(data);
  return 0;
}
