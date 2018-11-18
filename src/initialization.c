#include "../include/initialization.h"

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

void initialize_sizes(Cart_info cart_info, int *mp, int *np, int max_mp, int max_np) {
  if(!has_right(cart_info))
    *mp = max_mp;
  if(is_top(cart_info))
    *np = max_np;
}

void allocate_tables(double ***edge, double ***old, double ***new, int m, int n) {
  *edge = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *old = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *new = (double **) arralloc(sizeof(double), 2, m+2, n+2);
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


void initialization(int *world_rank, int *world_size, char *filename, int argc, char **argv, int *m, int *n, double ***masterbuf) {

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, world_size);

  if(argc == 2) {
    strcpy(filename, argv[1]);
  }
  else {
    strcpy(filename, "./resources/edgenew192x128.pgm");
  }

  pgmsize (filename, m, n);

  if(*world_rank == MASTER) {
    *masterbuf = (double **) arralloc(sizeof(double), 2, *m, *n);
    pgmread(filename, &(*masterbuf)[0][0], *m, *n);
  }

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
