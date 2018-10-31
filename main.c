/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "include/pgmio.h"

#define M 192
#define N 128

#define MAXITER   5000
#define PRINTFREQ  200

#define MASTER 0

double boundaryval(int i, int m);

void allocate_tables(double ***buf, double ***old, double ***new, double ***edge, int m, int n) {
  size_t size_of_rows = (m+2)*sizeof(double*);
  // size_t size_of_rows_buf = m*sizeof(double*);
  size_t size_of_cols = (n+2)*sizeof(double);
  // size_t size_of_cols_buf = n*sizeof(double);

  // *buf = (double**)malloc(size_of_rows_buf);
  *old = (double**)malloc(size_of_rows);
  *new = (double**)malloc(size_of_rows);
  *edge = (double**)malloc(size_of_rows);

    for (int i = 0; i < m+2; ++i)
    {
      // if(i<m)
      //   (*buf)[i] = (double*)malloc(size_of_cols_buf);
      (*old)[i] = (double*)malloc(size_of_cols);
      (*new)[i] = (double*)malloc(size_of_cols);
      (*edge)[i] = (double*)malloc(size_of_cols);
    }
}

int main (void) {
  // double old[M+2][N+2], new[M+2][N+2], edge[M+2][N+2];
  double **old, **new, **edge;
  // **buf, 
  double buf[M][N];

  int i, j, iter;
  char *filename;
  double val;

  int world_rank, world_size;
  double start_time, end_time;
  int m, n;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  printf("Process rank %d of %d processes\n", world_rank, world_size);
  printf("Number of iterations = %d\n", MAXITER);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    start_time = MPI_Wtime();

  filename = "resources/edgenew192x128.pgm";

  pgmsize (filename, &m, &n);
  printf("Processing %d x %d image\n", m, n);

  allocate_tables(&buf, &old, &new, &edge, m, n);


  printf("\nReading <%s>\n", filename);
  pgmread(filename, buf, m, n);
  printf("\n");

  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
      edge[i][j]=buf[i-1][j-1];
    }
  }

  for (i=0; i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      old[i][j]=255.0;
    }
  }

  /* Set fixed boundary conditions on the left and right sides */
  for (j=1; j < n+1; j++) {
    /* compute sawtooth value */
    val = boundaryval(j, n);

    old[0][j]   = (int)(255.0*(1.0-val));
    old[m+1][j] = (int)(255.0*val);
  }

  for (iter=1;iter<=MAXITER; iter++) {
    if(iter%PRINTFREQ==0) {
      printf("Iteration %d\n", iter);
    }

    /* Implement periodic boundary conditions on bottom and top sides */
    for (i=1; i < m+1; i++) {
      old[i][0]   = old[i][n];
      old[i][n+1] = old[i][1];
    }

    for (i=1;i<m+1;i++) {
      for (j=1;j<n+1;j++) {
        new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1] - edge[i][j]);
      }
    }

    for (i=1;i<m+1;i++) {
      for (j=1;j<n+1;j++) {
        old[i][j]=new[i][j];
      }
    }
  }

  printf("\nFinished %d iterations\n", iter-1);

  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
      buf[i-1][j-1]=old[i][j];
    }
  }

  filename="output/imagenew192x128.pgm";
  pgmwrite(filename, buf, m, n);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    end_time = MPI_Wtime();
    printf("Running time = %f\n", end_time - start_time);
  }

  MPI_Finalize();

  return 0;
} 

double boundaryval(int i, int m) {
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  return val;
}
