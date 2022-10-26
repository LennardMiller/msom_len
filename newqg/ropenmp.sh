#!/bin/bash
export OMP_NUM_THREADS=6
qcc -O3 qg.c -o qg.e -fopenmp -lm -lnetcdf
./qg.e

