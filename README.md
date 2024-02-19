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
number of matches played. See `matches.txt` for example.

To run the rating system:

    ./pingpong-rating-system < matches.txt

It should output something like:

    player raw normal matches wins losses
       richard +0.3515 1341 1 1 0
             x -0.0865  916 2 1 1
        jessie -0.2615  746 3 1 2
