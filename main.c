/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "include/arralloc.h"
#include "include/pgmio.h"

#define M 192
#define N 128

#define MAXITER   1000
#define PRINTFREQ  200

#define MASTER 0

double boundaryval(int i, int m) {
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  return val;
}

void allocate_tables(double ***old, double ***new, double ***edge, int m, int n) {
  *old = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *new = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *edge = (double **) arralloc(sizeof(double), 2, m+2, n+2);
}

void print_table(double **t, int m, int n) {
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%.1f ", t[i][j]);
    }
    printf("\n");
  }
}

void scatter_masterbuf(double **masterbuf, double **buf, int m, int n, int mp, int np, int max_mp, int max_np, int world_rank, int world_size, MPI_Comm comm) {
  int curr_coord[2];
  int send_rank;
  MPI_Status status;

  if(world_rank == MASTER)
    printf("MPI_Scatter\n");

  if(world_rank == MASTER) {
    double ***sendBuffs = (double ***) arralloc(sizeof(double), 3, world_size, max_mp, max_np);

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        curr_coord[0] = i/mp;
        curr_coord[1] = j/np;
        MPI_Cart_rank(comm, curr_coord, &send_rank);
        sendBuffs[send_rank][i%mp][j%np] = masterbuf[i][j];
      }
    }

    for (int i = 0; i < world_size; ++i) {
      MPI_Send(&(sendBuffs[i])[0][0], max_mp*max_np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }

  }

  MPI_Recv(&buf[0][0], max_mp*max_np, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
}

void gather_masterbuf(double **masterbuf, double **buf, int m, int n, int mp, int np, int max_mp, int max_np, int world_rank, int world_size, MPI_Comm comm) {
  int curr_coord[2];
  int send_rank;
  MPI_Status status;

  if(world_rank == MASTER)
    printf("MPI_Gather\n");

  MPI_Send(&buf[0][0], max_mp*max_np, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    double ***recvBuffs = (double ***) arralloc(sizeof(double), 3, world_size, max_mp, max_np);

    for (int i = 0; i < world_size; ++i) {
      MPI_Recv(&(recvBuffs[i])[0][0], max_mp*max_np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
    }

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        curr_coord[0] = i/mp;
        curr_coord[1] = j/np;
        MPI_Cart_rank(comm, curr_coord, &send_rank);
        masterbuf[i][j] = recvBuffs[send_rank][i%mp][j%np];
      }
    }
  }
}

void calculate(double **buf, double **old, double **new, double **edge, int m, int n) {
  int i, j, iter;
  double val;

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
      // printf("Iteration %d\n", iter);
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

  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
      buf[i-1][j-1]=old[i][j];
    }
  }
}

void free_tables(double **old, double **new, double **edge) {
  free(old);
  free(new);
  free(edge);
}

int main (void) {
    // double buf[M][N];
  // double old[M+2][N+2], new[M+2][N+2], edge[M+2][N+2];
  double **masterbuf = NULL, **buf, **old, **new, **edge;

  char *filename;

  int world_rank, world_size;
  double start_time, end_time;
  int m, n, mp, np, max_mp, max_np;

  MPI_Comm comm;
  int dim[2], period[2], reorder;
  int coord[2], id;
  
  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  // cart.c
  dim[0]=sqrt(world_size);
  dim[1]=world_size/dim[0];
  // what about world_size%dim[0] ???
  period[0]=1; period[1]=0;
  reorder=1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
  if(world_rank == MASTER)
  {
    printf("Number of iterations = %d\n", MAXITER);
    printf("Result topology: %d x %d\n", dim[0], dim[1]);
  }

  MPI_Cart_coords(comm, world_rank, 2, coord);
  printf("Rank %d coordinates are %d %d of %d processes\n", world_rank, coord[0], coord[1], world_size);
  MPI_Cart_rank(comm, coord, &id);
  printf("The processor at position (%d, %d) has rank %d\n", coord[0], coord[1], id);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    start_time = MPI_Wtime();

  filename = "resources/edgenew192x128.pgm";
  pgmsize (filename, &m, &n);

  //new sizes of tables
  mp = m/dim[0];
  if(coord[0] + 1 == dim[0]) {
    mp += m%dim[0];
  }      
  max_mp = m/dim[0] + m%dim[0];

  np = n/dim[1];
  if(coord[1] + 1 == dim[1]) {
    np += n%dim[1];
  }
  max_np = n/dim[1] + n%dim[1];

  printf("Rank: %d (%d, %d) Processing %d x %d = %d fraction of image %d x %d = %d\n", world_rank, coord[0], coord[1], mp, np, mp*np, m, n, m*n);
  if(world_rank == MASTER)
  {
    printf("Fractions: %d x %d x %d = %d of image %d\n", world_size, mp, np, world_size*mp*np, m*n);
  }
  allocate_tables(&old, &new, &edge, mp, np);

  MPI_Barrier(MPI_COMM_WORLD);

  masterbuf = (double **) arralloc(sizeof(double), 2, m, n);
  buf = (double **) arralloc(sizeof(double), 2, max_mp, max_np);

  if(world_rank == MASTER)
  {
    printf("\nReading <%s>\n", filename);
    pgmread(filename, &masterbuf[0][0], m, n);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  scatter_masterbuf(masterbuf, buf, m, n, mp, np, max_mp, max_np, world_rank, world_size, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    printf("calculate\n");
  calculate(buf, old, new, edge, mp, np);
  if(world_rank == MASTER)
    printf("Finished %d iterations\n", MAXITER);

  MPI_Barrier(MPI_COMM_WORLD);

  gather_masterbuf(masterbuf, buf, m, n, mp, np, max_mp, max_np, world_rank, world_size, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
  {
    filename="output/imagenew192x128.pgm";
    pgmwrite(filename, &masterbuf[0][0], m, n);
  }

  free(masterbuf);
  free(buf);

  free_tables(old, new, edge);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    end_time = MPI_Wtime();
    printf("Running time = %f\n", end_time - start_time);
  }

  MPI_Finalize();

  return 0;
} 
