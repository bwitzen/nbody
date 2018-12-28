#
# Makefile to compile my three NBody programs
# Ben Witzen <b.a.witzen@uva.nl>
# December 5th, 2015
#
# nbody-seq: Sequential
# nbody-par: Parallelized with OpenMPI
# nbody-omp: Parallelized with OpenMPI and OpenMP
#

all: nbody-seq nbody-par nbody-omp

nbody-seq:
	gcc -O3 -o nbody/nbody-seq nbody/nbody-seq.c -lm

nbody-par:
	mpicc -std=c99 -O3 -o nbody/nbody-par nbody/nbody-par.c nbody/ppm6.c -lm

nbody-omp:
	mpicc -std=c99 -O3 -o nbody/nbody-omp nbody/nbody-par.c nbody/ppm6.c -fopenmp -lm

clean:
	rm -f nbody/nbody-seq nbody/nbody-par nbody/nbody-omp
