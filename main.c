#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "include/arralloc.h"
#include "include/pgmio.h"
#include "include/cart_info.h"

#define MAXITER   1500
#define PRINTFREQ  100
#define FILENAME_SIZE 128

#define MIN_DIFF 0.075

#define MASTER 0

#define TOP_TO_BOTTOM 0
#define BOTTOM_TO_TOP 1
#define LEFT_TO_RIGHT 3
#define RIGHT_TO_LEFT 4

double boundaryval(int i, int m) {
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  return val;
}

void allocate_tables(double ***buf, double ***old, double ***new, double ***edge, int m, int n) {
  *buf = (double **) arralloc(sizeof(double), 2, m, n);
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

void scatter_masterbuf(double **masterbuf, double **buf, int m, int n, int mp, int np, int world_rank, int world_size, MPI_Comm comm) {
  int curr_coord[2];
  int send_rank;
  MPI_Status recv_status;
  MPI_Request recv_request;
  MPI_Status status[world_size];
  MPI_Request request[world_size];
  double ***sendBuffs = NULL;
  MPI_Datatype table;

  MPI_Type_contiguous(mp*np, MPI_DOUBLE, &table);
  MPI_Type_commit(&table);

  // if(world_rank == MASTER)
  //   printf("MPI_Scatter\n");

  MPI_Irecv(&buf[0][0], 1, table, MASTER, 0, comm, &recv_request);

  if(world_rank == MASTER) {
    sendBuffs = (double ***) arralloc(sizeof(double), 3, world_size, mp, np);

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {

        curr_coord[0] = i/mp;
        curr_coord[1] = j/np;
        MPI_Cart_rank(comm, curr_coord, &send_rank);

        sendBuffs[send_rank][i%mp][j%np] = masterbuf[i][j];
      }
    }

    for (int i = 0; i < world_size; ++i) {
      MPI_Isend(&(sendBuffs[i])[0][0], 1, table, i, 0, comm, &request[i]);
    }
    MPI_Waitall(world_size, request, status);
  }

  MPI_Wait(&recv_request, &recv_status);

  if(world_rank == MASTER)
    free(sendBuffs);
}

void gather_masterbuf(double **masterbuf, double **buf, int m, int n, int mp, int np, int world_rank, int world_size, MPI_Comm comm) {
  int curr_coord[2];
  int send_rank;
  MPI_Status send_status;
  MPI_Request send_request;
  MPI_Status status[world_size];
  MPI_Request request[world_size];
  double ***recvBuffs = NULL;
  MPI_Datatype table;

  MPI_Type_contiguous(mp*np, MPI_DOUBLE, &table);
  MPI_Type_commit(&table);
  // if(world_rank == MASTER)
  //   printf("MPI_Gather\n");

  MPI_Isend(&buf[0][0], 1, table, MASTER, 0, comm, &send_request);

  if(world_rank == MASTER) {
    recvBuffs = (double ***) arralloc(sizeof(double), 3, world_size, mp, np);

    for (int i = 0; i < world_size; ++i) {
      MPI_Irecv(&(recvBuffs[i])[0][0], 1, table, i, 0, comm, &request[i]);
    }

    MPI_Waitall(world_size, request, status);

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        curr_coord[0] = i/mp;
        curr_coord[1] = j/np;
        MPI_Cart_rank(comm, curr_coord, &send_rank);
        masterbuf[i][j] = recvBuffs[send_rank][i%mp][j%np];
      }
    }
  }

  MPI_Wait(&send_request, &send_status);

  if(world_rank == MASTER)
    free(recvBuffs);
}

void initialize_tables(double **buf, double **old, double **edge, int m, int n, Cart_info cart_info, int dim[2]) {
  int i, j;
  double val;

  // printf("initialize_tables: The processor at position (%d, %d) has rank %d\n", coord[0], coord[1], world_rank);

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

  for (j=0; j < n+2; j++) {
    /* compute sawtooth value */
    val = boundaryval(cart_info.coord[1]*n+j, dim[1]*n);

    if( (cart_info.coord[1] == 0 && j == 0) || (cart_info.coord[1] + 1 == dim[1] && j == n+1) ) {
        continue;
      }

      if(!has_left(cart_info))
        old[0][j]   = (int)(255.0*(1.0-val));
      if(!has_right(cart_info))
        old[m+1][j] = (int)(255.0*val);
  }

}

void halo_swaps(double **old, int m, int n, Cart_info cart_info, MPI_Datatype row, MPI_Datatype column, MPI_Request request[8], MPI_Status status[8]) {

  MPI_Irecv(&(old[1][n+1]), 1, column, cart_info.up, TOP_TO_BOTTOM, cart_info.comm, &request[0]);
  MPI_Isend(&(old[1][n]), 1, column, cart_info.up, BOTTOM_TO_TOP, cart_info.comm, &request[1]);

  MPI_Irecv(&(old[1][0]), 1, column, cart_info.down, BOTTOM_TO_TOP, cart_info.comm, &request[2]);
  MPI_Isend(&(old[1][1]), 1, column, cart_info.down, TOP_TO_BOTTOM, cart_info.comm, &request[3]);

  if(has_left(cart_info)) {
    MPI_Irecv(&(old[0][0]), 1, row, cart_info.left, LEFT_TO_RIGHT, cart_info.comm, &request[4]);
    MPI_Isend(&(old[1][0]), 1, row, cart_info.left, RIGHT_TO_LEFT, cart_info.comm, &request[5]);
  }

  if(has_right(cart_info)) {
    MPI_Irecv(&(old[m+1][0]), 1, row, cart_info.right, RIGHT_TO_LEFT, cart_info.comm, &request[6]);
    MPI_Isend(&(old[m][0]), 1, row, cart_info.right, LEFT_TO_RIGHT, cart_info.comm, &request[7]);
  }

}

void caclulate_halo(double **old, double **new, double **edge, int m, int n, Cart_info cart_info, MPI_Request request[8], MPI_Status status[8]) {
  int i, j;

  MPI_Waitall(4, request, status);

  if(has_left(cart_info))
    MPI_Waitall(2, &(request[4]), &(status[4]));

  if(has_right(cart_info))
    MPI_Waitall(2, &(request[6]), &(status[6]));

  // calculate swaps
  for (i=1;i<m+1;i+=m-1) {
    for (j=1;j<n+1;j++) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }

  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j+=n-1) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }
}

double calculate_max_diff(double **old, double **new, int m, int n) {
  double max_diff = -1, diff;

  for (int i=1;i<m+1;i++) {
    for (int j=1;j<n+1;j++) {
      if(new[i][j] - old[i][j] > 0) {
        diff = new[i][j] - old[i][j];
      }
      else {
        diff = old[i][j] - new[i][j];
      }

      if(max_diff == -1) max_diff = diff;
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

void calculate(double **buf, double **old, double **new, double **edge, int m, int n, Cart_info cart_info) {
  int i, j, iter;
  MPI_Datatype column, row;
  MPI_Request request[8];
  MPI_Status status[8];

  MPI_Type_contiguous(n+2, MPI_DOUBLE, &row);
  MPI_Type_commit(&row);
  MPI_Type_vector(m, 1, n+2, MPI_DOUBLE, &column);
  MPI_Type_commit(&column);

  // if(world_rank == MASTER)
  //   printf("calculate\n");

  for (iter=1; iter<=MAXITER; iter++) {
    // if(iter%PRINTFREQ==0) {
    //   printf("Iteration %d\n", iter);
    // }

    // Send and Receive halo swaps
    halo_swaps(old, m, n, cart_info, row, column, request, status);

    // calculate inner table
    for (i=2;i<m;i++) {
      for (j=2;j<n;j++) {
        new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
      }
    }

    // calculate borders with halo
    caclulate_halo(old, new, edge, m, n, cart_info, request, status);

    // if can_terminate every PRINTFREQ
    if(iter%PRINTFREQ==0) {
      if(can_terminate(old, new, m, n, cart_info.comm)) {
        // printf("ITER = %d\n", iter);
        // break;
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

  // if(world_rank == MASTER)
  //   printf("Finished %d iterations\n", MAXITER);
}

void free_tables(double **buf, double **old, double **new, double **edge) {
  free(buf);
  free(old);
  free(new);
  free(edge);
}

int main (int argc, char** argv) {
  double **masterbuf = NULL, **buf, **old, **new, **edge;

  char filename[FILENAME_SIZE];

  int world_rank, world_size, m, n, mp, np, dim[2];
  double start_time, end_time;
  int period[2] = {0, 1};
  int reorder = 1;
  MPI_Comm comm;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  dim[0]=sqrt(world_size);
  dim[1]=world_size/dim[0];

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

  Cart_info cart_info = discoverCart(world_rank, comm, dim);

  // print_cart_info(cart_info);

  if(argc == 2) {
    strcpy(filename, argv[1]);
  }
  else {
    strcpy(filename, "resources/");
    strcat(filename, "edgenew192x128.pgm");
  }

  pgmsize (filename, &m, &n);

  //new sizes of tables
  mp = m/dim[0];
  np = n/dim[1];

  if(world_rank == MASTER) {
    masterbuf = (double **) arralloc(sizeof(double), 2, m, n);
    // printf("\nReading <%s>\n", filename);
    pgmread(filename, &masterbuf[0][0], m, n);
  }

  allocate_tables(&buf, &old, &new, &edge, mp, np);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    start_time = MPI_Wtime();

  scatter_masterbuf(masterbuf, buf, m, n, mp, np, world_rank, world_size, comm);

  initialize_tables(buf, old, edge, mp, np, cart_info, dim);

  calculate(buf, old, new, edge, mp, np, cart_info);

  gather_masterbuf(masterbuf, buf, m, n, mp, np, world_rank, world_size, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    end_time = MPI_Wtime();
    // printf("Running time = %f\n", end_time - start_time);
  }

  if(world_rank == MASTER) {
    if(argc == 2) {
      strcpy(filename, "./output/image");
      strcat(filename, &(argv[1][16]));
    }
    else {
      strcpy(filename, "./output/");
      strcat(filename, "imagenew192x128.pgm");
    }
    pgmwrite(filename, &masterbuf[0][0], m, n);
    free(masterbuf);
  }

  free_tables(buf, old, new, edge);

  MPI_Finalize();

  return 0;
}
