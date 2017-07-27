
import sys
import random

# Number of players playing matches
n = int(sys.argv[1])

# Number of matches played
m = int(sys.argv[2])

# Generate random match results, where the probability of player i beating
# player j is i / (i+j)
for x in range(0, m):
  winner = random.randint(1, n)
  loser = winner
  while loser == winner:
    loser = random.randint(1, n)

  if random.uniform(0, 1) > (winner / float(winner + loser)):
    winner, loser = loser, winner

  print "p%i p%i" % (winner, loser)

