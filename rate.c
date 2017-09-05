
#include <stdio.h>    // for scanf, printf
#include <stdlib.h>   // for malloc, free
#include <string.h>   // for strcmp

#include "rate.h"

static size_t PlayerId(char* name, Data* data, size_t* capacity);
static Data* ReadData();
static void FreeData(Data* data);
static void PrintData(Data* data);

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

// PrintData --
//   Print the given data in human readable form for debugging purposes.
//
// Inputs:
//   data - the data to print.
//
// Result: 
//   none.
//
// Side effects:
//   Prints the data to stdout.
static void PrintData(Data* data)
{
  printf("%10s", "");
  for (size_t i = 0; i < data->n; ++i) {
    printf(" %10s", data->players[i]);
  }
  printf("\n");

  for (size_t i = 0; i < data->n; ++i) {
    printf("%10s", data->players[i]);
    for (size_t j = 0; j < data->n; ++j) {
      if (data->wins[i][j] > 0) {
        printf(" %10zd", data->wins[i][j]);
      } else {
        printf(" %10s", "");
      }
    }
    printf("\n");
  }
}

// PrintRatings --
//   Print the given ratings in human readable form.
//
// Inputs:
//   data - the match data.
//   ratings - the ratings to print.
//
// Result: 
//   none.
//
// Side effects:
//   Prints the ratings to stdout.
static void PrintRatings(Data* data, double* ratings)
{
  printf("== Ratings ==\n");
  for (size_t i = 0; i < data->n; ++i) {
    printf("%10s %1.4f\n", data->players[i], ratings[i]);
  }
}

int main() {
  Data* data = ReadData();
  PrintData(data);

  double ratings[data->n];

  printf("Prob2Rate...\n");
  Prob2Rate(data, ratings);
  PrintRatings(data, ratings);

  FreeData(data);
  return 0;
}

