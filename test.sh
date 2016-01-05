#!/bin/sh

rm tweetcube
g++ -lm -pthread -Ofast -march=native -Wall -funroll-loops -ffast-math -Wno-unused-result tweetcube.cpp -o tweetcube -lgsl -lm -lgslcblas
./tweetcube -workdir data -input cube -dim 4 -metric 2 -depth 5 -seedsize 3 -seeds 3 -threads 3 -iter 30
