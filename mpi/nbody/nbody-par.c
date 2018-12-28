/*******************************************************************************
 *
 * nbody-par.c
 * by Ben Witzen <b.a.witzen@uva.nl>
 * November 12th, 2015
 *
 * Simulates the N-Body problem using the MPI library for parallel execution.
 *
 * As a bonus assignment, I have also added OpenMP directives (pragmas) that can
 * be activated by compiling with the -fopenmp flag. On the DAS-4 cluster, it is
 * easiest to use gcc/4.4.6 for this purpose.
 *
 * This code was based on the provided sequential implementation, but I have
 * rewritten it from scratch to eliminate the use of global variables, remove
 * other obsolete functionalities, and to improve code readability.
 *
 * Usage: ./nbody-par bodies update-freq ppm steps
 *
 * Note that update-freq is not used, as this version does not draw updates back
 * to the PPM file. This argument is accounted for only to maintain backwards
 * compatibility with the sequential version and the automatic tests.
 *
 ******************************************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include "ppm6.h"

/**
 * Simulation constants.
 */
#define GRAVITY       1.1
#define FRICTION      0.01
#define DELTA_T       (0.025/5000)
#define	BOUNCE        -0.9
#define	DEFAULT_SEED  27102015
#define MPI_TAG       42

/**
 * Define a body structure and some macros to easily access elements within.
 */
typedef struct body {
  double x, y, xf, yf, xv, yv, mass, radius;
} body;

#define	X(B)          bodies[B].x
#define	Y(B)          bodies[B].y
#define	XF(B)         bodies[B].xf
#define	YF(B)         bodies[B].yf
#define	XV(B)         bodies[B].xv
#define	YV(B)         bodies[B].yv
#define	R(B)          bodies[B].radius
#define	M(B)          bodies[B].mass

/**
 * Gets width and height from a PPM image. This is a legacy function to provide
 * backward compatibility with the sequential code, which relied on this method
 * to obtain the dimensions of the space. Only PPM 6 is supported.
 *
 * char *infile : (IN) Filename to read from
 * int *xdim    : (OUT) X dimension of space
 * int *ydim    : (OUT) Y dimension of space
 */
bool get_dimensions_from_ppm(char *infile, int *xdim, int *ydim) {
  if (!map_P6(infile, xdim, ydim))
    return false;
  return true;
}

/**
 * Broadcast out my own data to the other nodes and receive broadcasts of their
 * data.
 *
 * bodies   : (IN/OUT) Pointer to start of entire bodies array.
 * mpi_rank : (IN) This node's mpi rank.
 * mpi_size : (IN) Total amount of nodes.
 * count    : (IN) Array of assigned amount of bodies to each node.
 */
void broadcast_positions(body *bodies, int mpi_rank, int mpi_size, int *count) {
  // each iteration either sends or receives multiple coordinate pairs
  for (int i = 0, offset = 0; i < mpi_size; i++) {
    // set up the coordinate pairs for this sender or recipient
    double data[count[i] * 2];
    
    // loop to "package" the data
    if (i == mpi_rank) {
      for (int b = offset, d = 0; b < offset + count[i]; b++) {
        data[d++] = X(b);
        data[d++] = Y(b);
      }
    }

    MPI_Bcast(data, count[i] * 2, MPI_DOUBLE, i, MPI_COMM_WORLD);

    // loop to "unpackage" the data
    if (i != mpi_rank) {
      for (int b = offset, d = 0; b < offset + count[i]; b++) {
        X(b) = data[d++];
        Y(b) = data[d++];
      }
    }
    
    offset += count[i];
  }
}

/**
 * Computes new values for all bodies local to this node. This function fuses
 * the simulation steps "clear forces", "compute forces", "compute velocities",
 * and "compute positions" together in one loop, which allows for better OpenMP
 * optimization.
 *
 * bodies : (IN/OUT) Pointer to start of entire bodies array.
 * len    : (IN) Total amount of bodies.
 * offset : (IN) Offset to bodies' array to point at this node's first body.
 * my_len : (IN) Amount of bodies this node has to work on.
 * xdim   : (IN) X dimension of space.
 * ydim   : (IN) Y dimension of space.
 */
void compute_timestep(body *bodies, int len, int offset, int my_len, int xdim, int ydim) {
  // storage for position updating
  double new_vals[2][len];

  #pragma omp parallel for schedule(static) if(my_len > 8)
  for (int b = offset; b < my_len + offset; b++) {
    // clear force accumulation
    YF(b) = XF(b) = 0;
    
    // compute forces
    for (int c = 0; c < len; c++) {
      if (b != c) {
        double dx = X(c) - X(b);
        double dy = Y(c) - Y(b);
        double angle = atan2(dy, dx);
        double dsqr = dx*dx + dy*dy;
        double mindist = R(b) + R(c);
        double mindsqr = mindist*mindist;
        double forced = ((dsqr < mindsqr) ? mindsqr : dsqr);
        double force = M(b) * M(c) * GRAVITY / forced;
        double xf = force * cos(angle);
        double yf = force * sin(angle);
        XF(b) += xf;
        YF(b) += yf;
      }
    }
    
    // compute velocities
    double xv = XV(b);
    double yv = YV(b);
    double force = sqrt(xv*xv + yv*yv) * FRICTION;
    double angle = atan2(yv, xv);
    double xf = XF(b) - (force * cos(angle));
    double yf = YF(b) - (force * sin(angle));
    XV(b) += (xf / M(b)) * DELTA_T;
    YV(b) += (yf / M(b)) * DELTA_T;
    
    // compute positions    
    double xn = X(b) + (XV(b) * DELTA_T);
    double yn = Y(b) + (YV(b) * DELTA_T);

    // bounce of walls
    if (xn < 0) {
      xn = 0;
      XV(b) = -XV(b);
    }
    else if (xn >= xdim) {
      xn = xdim - 1;
      XV(b) = -XV(b);
    }

    if (yn < 0) {
      yn = 0;
      YV(b) = -YV(b);
    }
    else if (yn >= ydim) {
      yn = ydim - 1;
      YV(b) = -YV(b);
    }

    // store new positions
    new_vals[0][b] = xn;
    new_vals[1][b] = yn;
  }
  
  // flip old positions for new ones
  for (int b = offset, n = my_len + offset; b < n; b++) {
    X(b) = new_vals[0][b];
    Y(b) = new_vals[1][b];
  }
}

/**
 * Main. Handles arguments, sets up MPI, and contains the code for slaves and
 * master to execute.
 */
int main(int argc, char **argv) {
  // initialize the MPI environment
  int mpi_rank, mpi_size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // runtime parameters
  int body_count, steps, xdim, ydim;
  int seed = DEFAULT_SEED;

  // macro to exit cleanly if invalid arguments are encountered
  #define ARG_FAIL(x) { if (mpi_rank == 0) \
                          fprintf(stderr, "\e[31mArg error:\e[0m " x "\n"); \
                        MPI_Finalize(); \
                        exit(EXIT_SUCCESS); }

  // handle arguments
  if (argc == 5) {
    body_count = atoi(argv[1]);
    steps = atoi(argv[4]);
    if (!get_dimensions_from_ppm(argv[3], &xdim, &ydim))
      ARG_FAIL("Cannot read your PPM #6 file.");
  }
  else
    ARG_FAIL("Incorrect usage: ./nbody-par bodies update-freq ppm steps");

  // perform a very simple sanity check on the arguments
  if (body_count <= 0 || steps < 0 || xdim <= 0 || ydim <= 0)
    ARG_FAIL("One or more parameters is invalid.");

  if (mpi_size > body_count)
    ARG_FAIL("It is not allowed to have more nodes than bodies.");

  #undef ARG_FAIL

  if (mpi_size < 2)
    fprintf(stderr, "\e[33mNote:\e[0m Only one node is active!\n");

  // compute how many bodies each participant will be calculating
  int local_body_count[mpi_size];
  memset(local_body_count, 0, mpi_size * sizeof(int));

  for (int i = 0; i < body_count; i++)
    local_body_count[i % mpi_size]++;

  // allocate space for all bodies on all nodes
  // this is needed for initialization and to store received broadcasting info
  body *bodies = malloc(body_count * sizeof(body));
  if (!bodies) {
    fprintf(stderr, "\e[31mError:\e[0m Malloc failure at node %d.\n", mpi_rank);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  // define timer moments
  struct timeval start, end;

  // master code
  if (mpi_rank == 0) {
    // initialize bodies
    srand(seed);
    for (int b = 0; b < body_count; b++) {
      X(b) = (rand() % xdim);
      Y(b) = (rand() % ydim);
      R(b) = 1 + ((b*b + 1.0) * sqrt(1.0 * ((xdim * xdim) + (ydim * ydim)))) /
             (25.0 * (body_count * body_count + 1.0));
      M(b) = R(b) * R(b) * R(b);
      XV(b) = ((rand() % 20000) - 10000) / 2000.0;
      YV(b) = ((rand() % 20000) - 10000) / 2000.0;
    }

    // send the entire initial state to other nodes
    for (int dest = 1; dest < mpi_size; dest++)
      if (local_body_count[dest] > 0)
        MPI_Send(bodies, body_count * sizeof(body), MPI_BYTE, dest, MPI_TAG, MPI_COMM_WORLD);

    // prepare my length and data offset
    int length = local_body_count[mpi_rank];
    int offset = 0;

    if (gettimeofday(&start, 0) != 0) {
      fprintf(stderr, "\e[31mError:\e[0m Could not do timing.\n");
      exit(EXIT_FAILURE);
    }

    // perform local computation
    while (steps--) {
      broadcast_positions(bodies, mpi_rank, mpi_size, local_body_count);
      compute_timestep(bodies, body_count, offset, length, xdim, ydim);
    }

    // stop timer and calculate time spent -- since the master is guaranteed to
    // have the most nodes, stopping the timer here is justified
    if (gettimeofday(&end, 0) != 0) {
      fprintf(stderr, "\e[31mError:\e[0m Could not do timing.\n");
      exit(EXIT_FAILURE);
    }

    // receive final results
    for (int src = 1, offset = 0; src < mpi_size; src++) {
      offset += local_body_count[src - 1];
      if (local_body_count[src] > 0) {
        MPI_Recv(bodies + offset,
                 local_body_count[src] * sizeof(body),
                 MPI_BYTE,
                 src,
                 MPI_TAG,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    }

    // print final body state
    for (int b = 0; b < body_count; b++)
      printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
             X(b), Y(b), XF(b), YF(b), XV(b), YV(b));
  
    // print time used
    double time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
    fprintf(stderr, "\e[32mSuccess:\e[0m Worker %d took %10.3f seconds.\n", mpi_rank, time);
  }

  // slave code
  else {
    // prepare my length and data offset
    int length = local_body_count[mpi_rank];
    int offset = 0;
    for (int i = 0; i < mpi_rank; i++)
      offset += local_body_count[i];

    // receive initial data for my bodies
    MPI_Recv(bodies, body_count * sizeof(body), MPI_BYTE, 0, MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (gettimeofday(&start, 0) != 0) {
      fprintf(stderr, "\e[31mError:\e[0m Could not do timing.\n");
      exit(EXIT_FAILURE);
    }

    // perform local computation
    while (steps--) {
      broadcast_positions(bodies, mpi_rank, mpi_size, local_body_count);
      compute_timestep(bodies, body_count, offset, length, xdim, ydim);
    }

    // stop timer and calculate time spent -- since the master is guaranteed to
    // have the most nodes, stopping the timer here is justified
    if (gettimeofday(&end, 0) != 0) {
      fprintf(stderr, "\e[31mError:\e[0m Could not do timing.\n");
      exit(EXIT_FAILURE);
    }

    // send final results to master
    MPI_Send(bodies + offset, local_body_count[mpi_rank] * sizeof(body), MPI_BYTE, 0, MPI_TAG, MPI_COMM_WORLD);
    
    // print time used
    double time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
    fprintf(stderr, "\e[32mSuccess:\e[0m Worker %d took %10.3f seconds.\n", mpi_rank, time);
  }

  // cleanup
  free(bodies);
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
