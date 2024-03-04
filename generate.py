
import math
import random
import sys

# Number of players playing matches
n = int(sys.argv[1])

# Number of matches played
m = int(sys.argv[2])

def prob(a, b):
  return 1 / (1 + math.exp(-(a-b)))

# Randomly pick ratings for n players, normally distributed.
sigma = 1.0
ratings = []
for x in range(0, n):
  ratings.append(random.gauss(0.0, sigma))
  sys.stderr.write("p%i %1.5f\n" % (x, ratings[x]))

# Randomly generate match results for those players
for x in range(0, m):
  winner = random.randint(0, n-1)
  loser = winner
  while loser == winner:
    loser = random.randint(0, n-1)

  if random.uniform(0, 1) > prob(ratings[winner], ratings[loser]):
    winner, loser = loser, winner

  print("p%i p%i" % (winner, loser))

