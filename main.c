#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <assert.h>

#include "include/arralloc.h"
#include "include/pgmio.h"
#include "include/mpi_types.h"
#include "include/cart_info.h"

#define MAXITER   1000000
#define PRINTFREQ  100
#define FILENAME_SIZE 128

#define MIN_DIFF 0.075

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

}

void allocate_tables(double ***edge, double ***old, double ***new, int m, int n) {
  *edge = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *old = (double **) arralloc(sizeof(double), 2, m+2, n+2);
  *new = (double **) arralloc(sizeof(double), 2, m+2, n+2);
}

void scatter_masterbuf(double **masterbuf, double **edge, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes) {
  int curr_coord[2], send_rank;
  MPI_Status recv_status, status[cart_info.world_size];
  MPI_Request recv_request, request[cart_info.world_size];

  if(!has_right(cart_info))
    MPI_Irecv(&edge[1][1], 1, mpi_Datatypes.max_cont_table, MASTER, 0, cart_info.comm, &recv_request);
  else
    MPI_Irecv(&edge[1][1], 1, mpi_Datatypes.cont_table, MASTER, 0, cart_info.comm, &recv_request);

  if(cart_info.id == MASTER) {

    for (int i = 0; i < cart_info.dim[0]; i++) {
      curr_coord[0] = i;
      for (int j = 0; j < cart_info.dim[1]; j++) {
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

void gather_masterbuf(double **masterbuf, double **edge, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes) {
  int curr_coord[2], recv_rank;
  MPI_Status send_status, status[cart_info.world_size];
  MPI_Request send_request, request[cart_info.world_size];

  if(!has_right(cart_info))
    MPI_Isend(&edge[1][1], 1, mpi_Datatypes.max_cont_table, MASTER, 0, cart_info.comm, &send_request);
  else
    MPI_Isend(&edge[1][1], 1, mpi_Datatypes.cont_table, MASTER, 0, cart_info.comm, &send_request);


  if(cart_info.id == MASTER) {

    for (int i = 0; i < cart_info.dim[0]; i++) {
      curr_coord[0] = i;
      for (int j = 0; j < cart_info.dim[1]; j++) {
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

void initialize_tables(double **edge, double **old, int m, int n, Cart_info cart_info) {
  int i, j;
  double val;

  for (i=0; i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      old[i][j]=255.0;
    }
  }

  /* compute sawtooth value */
  if(!has_left(cart_info))
    i = 0;
  if(!has_right(cart_info))
    i = m+1;

  for (j=0; j<n+2; j++) {
    val = boundaryval(cart_info.coord[1]*n+j, cart_info.dim[1]*n);

    if(i == 0)
      old[i][j] = (int)(255.0*(1.0-val));
    if(i == m+1)
      old[i][j] = (int)(255.0*val);
  }

}

double get_average_pixel(double **table, MPI_Comm comm, int m, int n, int mp, int np) {
  double sum = 0, global_sum = 0;

  for (int i=1; i<mp+1; i++) {
    for (int j=1; j<np+1; j++) {
      sum += table[i][j];
    }
  }
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  return global_sum/(m*n);
}

void print_average_pixel(double **table, int iter, Cart_info cart_info, int m, int n, int mp, int np, char *filename) {
  double average_pixel;
  FILE * fp;
  char average_pixel_filename[FILENAME_SIZE];

  average_pixel = get_average_pixel(table, cart_info.comm, m, n, mp, np);
  if(cart_info.id == MASTER) {
    strcpy(average_pixel_filename, "./data/");
    strncat(average_pixel_filename, filename, 14);
    strcat(average_pixel_filename, "_average_pixel.tsv");
    if(iter == 1) {
      fp = fopen (average_pixel_filename, "w");
      fprintf(fp, "iteration\taverage pixel\n");
      fclose(fp);
    }
    fp = fopen (average_pixel_filename, "a");
    fprintf(fp, "%d\t%f\n", iter, average_pixel);
    fclose(fp);
  }
}

void halo_swaps(double **old, int m, int n, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, MPI_Request request[8], MPI_Status status[8]) {

  MPI_Irecv(&(old[1][n+1]), 1, mpi_Datatypes->column, cart_info.up, TOP_TO_BOTTOM, cart_info.comm, &request[0]);
  MPI_Isend(&(old[1][n]), 1, mpi_Datatypes->column, cart_info.up, BOTTOM_TO_TOP, cart_info.comm, &request[1]);

  MPI_Irecv(&(old[1][0]), 1, mpi_Datatypes->column, cart_info.down, BOTTOM_TO_TOP, cart_info.comm, &request[2]);
  MPI_Isend(&(old[1][1]), 1, mpi_Datatypes->column, cart_info.down, TOP_TO_BOTTOM, cart_info.comm, &request[3]);

  MPI_Irecv(&(old[0][0]), 1, mpi_Datatypes->row, cart_info.left, LEFT_TO_RIGHT, cart_info.comm, &request[4]);
  MPI_Isend(&(old[1][0]), 1, mpi_Datatypes->row, cart_info.left, RIGHT_TO_LEFT, cart_info.comm, &request[5]);

  MPI_Irecv(&(old[m+1][0]), 1, mpi_Datatypes->row, cart_info.right, RIGHT_TO_LEFT, cart_info.comm, &request[6]);
  MPI_Isend(&(old[m][0]), 1, mpi_Datatypes->row, cart_info.right, LEFT_TO_RIGHT, cart_info.comm, &request[7]);

}

void caclulate_halo(double **old, double **new, double **edge, int m, int n, Cart_info cart_info, MPI_Request request[8], MPI_Status status[8]) {
  int i, j;

  MPI_Waitall(8, request, status);

  // calculate halo
  for (i=1;i<m+1;i+=m-1) {
    for (j=1;j<n+1;j++) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }

  for (j=1;j<n+1;j+=n-1) {
    for (i=1;i<m+1;i++) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }
}

double calculate_max_diff(double **old, double **new, int m, int n) {
  double max_diff = -1, diff;

  for (int i=1;i<m+1;i++) {
    for (int j=1;j<n+1;j++) {
      diff = fabs(old[i][j] - new[i][j]);
      if(diff > max_diff)  max_diff = diff;
      if(max_diff >= MIN_DIFF) break;
    }
  }

  return max_diff;
}

int can_terminate(double **old, double **new, int m, int n, MPI_Comm comm) {
  double max_diff, global_max_diff;

  max_diff = calculate_max_diff(old, new, m, n);

  // all reduce max_diff
  MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, comm);

  if(global_max_diff < MIN_DIFF && global_max_diff >= 0)
    return 1;

  return 0;
}

void calculate(double **edge, double **old, double **new, int m, int n, int mp, int np, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, char *filename) {
  int i, j, iter;
  MPI_Request request[8];
  MPI_Status status[8];

  init_mpi_datatypes_row_col(mpi_Datatypes, mp, np);

  for (iter=1; iter<=MAXITER; iter++) {

    // Send and Receive halo swaps
    halo_swaps(old, mp, np, cart_info, mpi_Datatypes, request, status);

    // calculate inner table
    for (i=2;i<mp;i++) {
      for (j=2;j<np;j++) {
        new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
      }
    }

    // calculate borders with halo
    caclulate_halo(old, new, edge, mp, np, cart_info, request, status);

    // if can_terminate every PRINTFREQ
    if(iter%PRINTFREQ==0 || iter == 1) {
      print_average_pixel(new, iter, cart_info, m, n, mp, np, filename);

      if(can_terminate(old, new, mp, np, cart_info.comm))
        break;
    }

    for (i=1;i<mp+1;i++) {
      for (j=1;j<np+1;j++) {
        old[i][j]=new[i][j];
      }
    }

  }

  for (i=1;i<mp+1;i++) {
    for (j=1;j<np+1;j++) {
      edge[i][j]=old[i][j];
    }
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
  double start_time, end_time;
  int dim[2], period[2] = {0, 1}, reorder = 1;
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

  scatter_masterbuf(masterbuf, edge, mp, np, cart_info, mpi_Datatypes);

  initialize_tables(edge, old, mp, np, cart_info);

  calculate(edge, old, new, m, n, mp, np, cart_info, &mpi_Datatypes, &filename[12]);

  gather_masterbuf(masterbuf, edge, mp, np, cart_info, mpi_Datatypes);

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
