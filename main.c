#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "include/arralloc.h"
#include "include/pgmio.h"
#include "include/mpi_types.h"
#include "include/cart_info.h"
#include "include/calculate.h"

#define FILENAME_SIZE 128
#define MASTER 0

double boundaryval(int i, int m) {
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  return val;
}

void decomposition(int world_size, int m, int n, int dim[2], int *mp, int *np, int *max_mp, int *max_np) {

  MPI_Dims_create(world_size, 2, dim);

  //new sizes of tables
  *mp = m/dim[0];
  *np = n/dim[1];
  *max_mp = m/dim[0] + m%dim[0];
  *max_np = n/dim[1] + n%dim[1];

  if(n%dim[1] != 0) {
    dim[0] = 0;
    dim[1] = 1;
    decomposition(world_size, m, n, dim, mp, np, max_mp, max_np);
  }
}

void allocate_tables(double ***edge, double ***old, double ***new, int m, int n) {
  *edge = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *old = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *new = (double **) arralloc(sizeof(double), 2, m+2, n+2);
}

void scatter_masterbuf(double **masterbuf, double **edge, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes) {
  int curr_coord[2], send_rank, i, j;
  MPI_Status recv_status, status[cart_info.world_size];
  MPI_Request recv_request, request[cart_info.world_size];

  if(!has_right(cart_info))
    MPI_Irecv(&edge[1][1], 1, mpi_Datatypes.max_cont_table, MASTER, 0, cart_info.comm, &recv_request);
  else
    MPI_Irecv(&edge[1][1], 1, mpi_Datatypes.cont_table, MASTER, 0, cart_info.comm, &recv_request);

  if(cart_info.id == MASTER) {

    for (i = 0; i < cart_info.dim[0]; i++) {
      curr_coord[0] = i;
      for (j = 0; j < cart_info.dim[1]; j++) {
        curr_coord[1] = j;
        MPI_Cart_rank(cart_info.comm, curr_coord, &send_rank);
        if(i + 1 == cart_info.dim[0])
          MPI_Isend(&(masterbuf[i*mp][j*np]), 1, mpi_Datatypes.max_table, send_rank, 0, cart_info.comm, &request[send_rank]);
        else
          MPI_Isend(&(masterbuf[i*mp][j*np]), 1, mpi_Datatypes.table, send_rank, 0, cart_info.comm, &request[send_rank]);
      }
    }
    MPI_Waitall(cart_info.world_size, request, status);
  }

  MPI_Wait(&recv_request, &recv_status);
}

void gather_masterbuf(double **masterbuf, double **old, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes) {
  int curr_coord[2], recv_rank, i, j;
  MPI_Status send_status, status[cart_info.world_size];
  MPI_Request send_request, request[cart_info.world_size];

  if(!has_right(cart_info))
    MPI_Isend(&old[1][1], 1, mpi_Datatypes.max_cont_table, MASTER, 0, cart_info.comm, &send_request);
  else
    MPI_Isend(&old[1][1], 1, mpi_Datatypes.cont_table, MASTER, 0, cart_info.comm, &send_request);

  if(cart_info.id == MASTER) {

    for (i = 0; i < cart_info.dim[0]; i++) {
      curr_coord[0] = i;
      for (j = 0; j < cart_info.dim[1]; j++) {
        curr_coord[1] = j;
        MPI_Cart_rank(cart_info.comm, curr_coord, &recv_rank);
        if(i + 1 == cart_info.dim[0])
          MPI_Irecv(&(masterbuf[i*mp][j*np]), 1, mpi_Datatypes.max_table, recv_rank, 0, cart_info.comm, &request[recv_rank]);
        else
          MPI_Irecv(&(masterbuf[i*mp][j*np]), 1, mpi_Datatypes.table, recv_rank, 0, cart_info.comm, &request[recv_rank]);
      }
    }

    MPI_Waitall(cart_info.world_size, request, status);
  }

  MPI_Wait(&send_request, &send_status);
}

void initialize_tables(double **old, int m, int n, Cart_info cart_info) {
  int i, j;
  double val;

  for (i=0; i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      old[i][j]=255.0;
    }
  }

  /* compute sawtooth value */
  for (j=0; j<n+2; j++) {
    val = boundaryval(cart_info.coord[1]*n+j, cart_info.dim[1]*n);

    if(!has_left(cart_info))
      old[0][j] = (int)(255.0*(1.0-val));
    if(!has_right(cart_info))
      old[m+1][j] = (int)(255.0*val);
  }

}

void free_tables(double **edge, double **old, double **new) {
  free(edge);
  free(old);
  free(new);
}

int main (int argc, char** argv) {
  double **masterbuf = NULL, **edge, **old, **new;
  FILE * fp;
  char filename[FILENAME_SIZE];
  int world_rank, world_size, m, n, mp, np, max_mp, max_np;
  double start_time = 0, end_time = 0;
  int dim[2] = {0, 0}, period[2] = {0, 1}, reorder = 1;
  MPI_Comm comm;
  Cart_info cart_info;
  Mpi_Datatypes mpi_Datatypes;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(argc == 2) {
    strcpy(filename, argv[1]);
  }
  else {
    strcpy(filename, "./resources/edgenew192x128.pgm");
  }

  pgmsize (filename, &m, &n);

  decomposition(world_size, m, n, dim, &mp, &np, &max_mp, &max_np);

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

  cart_info = discoverCart(world_rank, comm, world_size, dim);

  mpi_Datatypes = init_mpi_datatypes(n, mp, np, max_mp, max_np);

  if(!has_right(cart_info))
    mp = max_mp;
  if(is_top(cart_info))
    np = max_np;

  if(world_rank == MASTER) {
    masterbuf = (double **) arralloc(sizeof(double), 2, m, n);
    pgmread(filename, &masterbuf[0][0], m, n);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    start_time = MPI_Wtime();

  allocate_tables(&edge, &old, &new, mp, np);

  initialize_tables(old, mp, np, cart_info);

  scatter_masterbuf(masterbuf, edge, mp, np, cart_info, mpi_Datatypes);

  calculate(edge, old, new, m, n, mp, np, cart_info, &mpi_Datatypes, &filename[12]);

  gather_masterbuf(masterbuf, old, mp, np, cart_info, mpi_Datatypes);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    end_time = MPI_Wtime();
    fp = fopen ("./data/results.tsv", "a");
    fprintf(fp, "%s\t%d\t%f\n", filename, world_size, end_time - start_time);
    fclose(fp);

    if(argc == 2) {
      strcpy(filename, "./output/image");
      strcat(filename, &(argv[1][16]));
    }
    else {
      strcpy(filename, "./output/imagenew192x128.pgm");
    }
    pgmwrite(filename, &masterbuf[0][0], m, n);
    free(masterbuf);
  }

  free_tables(edge, old, new);

  MPI_Finalize();

  return 0;
}
