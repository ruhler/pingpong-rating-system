
#ifndef RATE_H_
#define RATE_H_

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

// XXXRate --
//   Rate players using the XXX algorithm.
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
void Prob2Rate(Data* data, double ratings[]);

#endif//RATE_H_
