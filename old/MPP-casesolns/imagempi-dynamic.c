/*
 * Simple solution to Case Study exercise from the EPCC MPI course.
 * Note that this uses dynamic memory allocation for the 2D arrays.
 *
 * Communications is done using the sendrecv routine; a proper
 * solution would use non-blocking communications (eg some combination
 * of issend/recv and ssend/irecv).
 *
 * Note that the special rank of MPI_PROC_NULL is a "black hole" for
 * communications (similar to /dev/null in Unix). Using this value for
 * processes off the edges of the image means we do not need any
 * additional logic to ensure processes at the edges do not attempt to
 * send to or receive from invalid ranks (ie rank = -1 and rank = P).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "arralloc.h"
#include "pgmio.h"

#define M 192
#define N 128

#define P 4

#define MP M/P
#define NP N

#define MAXITER   1500
#define PRINTFREQ  200

int main (int argc, char **argv)
{
  double **old, **new, **edge, **masterbuf, **buf;

  int i, j, iter, maxiter;
  char *filename;

  int rank, size, next, prev;
  MPI_Status status;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(size != P)
    {
      if (rank == 0) printf("ERROR: size = %d, P = %d\n", size, P);
      MPI_Finalize();
      exit(-1);
    }

  /*
   * Allocate arrays dynamically. In a real code you would enquire the
   * image dimensions at runtime using "pgmsize" and then compute
   * local array sizes; this would remove the requirement for P to be
   * defined at compile time, so the executable would run for any P.
   */

  masterbuf = (double **) arralloc(sizeof(double), 2, M , N );
  buf       = (double **) arralloc(sizeof(double), 2, MP, NP);

  new  = (double **) arralloc(sizeof(double), 2, MP+2, NP+2);
  old  = (double **) arralloc(sizeof(double), 2, MP+2, NP+2);
  edge = (double **) arralloc(sizeof(double), 2, MP+2, NP+2);

  next = rank + 1;
  prev = rank - 1;

  if (next >= size)
    {
      next = MPI_PROC_NULL;
    }

  if (prev < 0)
    {
      prev = MPI_PROC_NULL;
    }

  if(rank == 0)
    {
      printf("Processing %d x %d image on %d processes\n", M, N, P);
      printf("Number of iterations = %d\n", MAXITER);

      filename = "edge192x128.pgm";

      printf("\nReading <%s>\n", filename);
      pgmread(filename, &masterbuf[0][0], M, N);
      printf("\n");

    }

  MPI_Scatter(&masterbuf[0][0], MP*NP, MPI_DOUBLE,
	      &buf[0][0],       MP*NP, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  for (i=1;i<MP+1;i++)
    {
      for (j=1;j<NP+1;j++)
	{
	  edge[i][j]=buf[i-1][j-1];
	}
    }

  for (i=0;i<MP+2;i++)
    {
      for (j=0;j<NP+2;j++)
	{
	  old[i][j]=255.0;
	}
    }

  for (iter=1;iter<=MAXITER; iter++)
    {
      if(iter%PRINTFREQ==0)
	{
	  if(rank==0)
	    {
	      printf("Iteration %d\n", iter);
	    }
	}

      MPI_Sendrecv(&old[MP][1], NP, MPI_DOUBLE, next, 1, 
		   &old[0][1],  NP, MPI_DOUBLE, prev, 1,
		   MPI_COMM_WORLD, &status);

      MPI_Sendrecv(&old[1][1],    NP, MPI_DOUBLE, prev, 2, 
		   &old[MP+1][1], NP, MPI_DOUBLE, next, 2,
		   MPI_COMM_WORLD, &status);

      for (i=1;i<MP+1;i++)
	{
	  for (j=1;j<NP+1;j++)
	    {
	      new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
			      - edge[i][j]);
	    }
	}
	
      for (i=1;i<MP+1;i++)
	{
	  for (j=1;j<NP+1;j++)
	    {
	      old[i][j]=new[i][j];
	    }
	}
    }

  if (rank==0)
    {
      printf("\nFinished %d iterations\n", iter-1);
    }

  for (i=1;i<MP+1;i++)
    {
      for (j=1;j<NP+1;j++)
	{
	  buf[i-1][j-1]=old[i][j];
	}
    }

  MPI_Gather(&buf[0][0],       MP*NP, MPI_DOUBLE,
	     &masterbuf[0][0], MP*NP, MPI_DOUBLE,
	     0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      filename="image192x128.pgm";
      printf("\nWriting <%s>\n", filename); 
      pgmwrite(filename, &masterbuf[0][0], M, N);
    }

  free(masterbuf);
  free(buf);
  free(new);
  free(old);
  free(edge);

  MPI_Finalize();
} 
