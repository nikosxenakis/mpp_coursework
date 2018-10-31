/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "include/pgmio.h"

#define M 192
#define N 128

#define MAXITER   50000
#define PRINTFREQ  200

double boundaryval(int i, int m);

int main (void) {
  double old[M+2][N+2], new[M+2][N+2], edge[M+2][N+2];
  double buf[M][N];

  int i, j, iter;
  char *filename;
  double val;

  printf("Processing %d x %d image\n", M, N);
  printf("Number of iterations = %d\n", MAXITER);

  filename = "resources/edgenew192x128.pgm";

  printf("\nReading <%s>\n", filename);
  pgmread(filename, buf, M, N);
  printf("\n");

  for (i=1;i<M+1;i++) {
    for (j=1;j<N+1;j++) {
      edge[i][j]=buf[i-1][j-1];
    }
  }

  for (i=0; i<M+2;i++) {
    for (j=0;j<N+2;j++) {
      old[i][j]=255.0;
    }
  }

  /* Set fixed boundary conditions on the left and right sides */
  for (j=1; j < N+1; j++) {
    /* compute sawtooth value */
    val = boundaryval(j, N);

    old[0][j]   = (int)(255.0*(1.0-val));
    old[M+1][j] = (int)(255.0*val);
  }

  for (iter=1;iter<=MAXITER; iter++) {
    if(iter%PRINTFREQ==0) {
      printf("Iteration %d\n", iter);
    }

    /* Implement periodic boundary conditions on bottom and top sides */
    for (i=1; i < M+1; i++) {
      old[i][0]   = old[i][N];
      old[i][N+1] = old[i][1];
    }

    for (i=1;i<M+1;i++) {
      for (j=1;j<N+1;j++) {
        new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1] - edge[i][j]);
      }
    }

    for (i=1;i<M+1;i++) {
      for (j=1;j<N+1;j++) {
        old[i][j]=new[i][j];
      }
    }
  }

  printf("\nFinished %d iterations\n", iter-1);

  for (i=1;i<M+1;i++) {
    for (j=1;j<N+1;j++) {
      buf[i-1][j-1]=old[i][j];
    }
  }

  filename="output/imagenew192x128.pgm";
  printf("\nWriting <%s>\n", filename); 
  pgmwrite(filename, buf, M, N);

  return 0;
} 

double boundaryval(int i, int m) {
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  return val;
}
