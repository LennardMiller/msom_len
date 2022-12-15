CC99='mpicc -std=c99' qcc -O3 -D_MPI=4 qg.c -o qg.e -lm -lnetcdf
mpirun -np 4 ./qg.e


