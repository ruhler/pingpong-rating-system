# Ping Pong Rating System README

A ping pong rating system for rating players in a local league. Includes a
specification and implementation of the rating system.

## Source

The source code is available from github at
<https://github.com/ruhler/pingpong-rating-system>.

## Build Requirements

To install required dependencies on a debian based system:

    apt install gcc texlive-latex-base texlive-pictures

## Build

To build, type `make`.

This generates two items:

* **pingpong-rating-system.pdf**: A description of the rating system.

* **pingpong-rating-system**: An implementation of the rating system.

## Run

The ping pong rating system takes a file as input listing match results. Each
line lists a match winner, followed by a match loser, followed by an optional
number of matches played. See `m1.txt` and `m2.txt` for example.

To run the rating system:

    ./pingpong-rating-system < m1.txt

It should output something like:

    player raw normal matches wins losses
    richard +0.3515 1341 1 1 0
    x -0.0865 916 2 1 1
    jessie -0.2615 746 3 1 2

### Seeding

You can provide seed rating values for players by providing a file with the
`--seed` option to the program listing initial player seeds. The format of the
seed file is one entry per line, with the player name followed by a space
followed by the raw rating value for that player. You could, for example,
provide a previously output ratings result as the seed.

The purpose of seeding is to speed up the rating computation when you have a
good idea of the player ratings already.

For example, generate a seed rating from m2.txt using:

   ./pingpong-rating-system < m2.txt > m2.ratings.txt

Then rerun the computation using the result as the seed:

   ./pingpong-rating-system --seed m2.ratings.txt < m2.txt

It runs much faster the second time.

