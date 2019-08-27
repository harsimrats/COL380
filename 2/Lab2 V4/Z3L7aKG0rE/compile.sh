#!/bin/sh

g++ -o run.out -fopenmp graph_partition.cpp main.cpp -std=c++11 -O3 -g