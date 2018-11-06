#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "include/arralloc.h"
#include "include/pgmio.h"

#define MAXITER   1500
#define PRINTFREQ  MAXITER/100
#define FILENAME_SIZE 128

#define MASTER 0

typedef struct Cart_info {
  int id;
  MPI_Comm comm;
  int coord[2];
  int up;
  int down;
  int left;
  int right;
} Cart_info;

int has_top(Cart_info cart_info) {
  return 1;
}

int has_bottom(Cart_info cart_info) {
  return 1;
}

int has_left(Cart_info cart_info) {
  return cart_info.left != -1;
}

int has_right(Cart_info cart_info) {
  return cart_info.right != -1;
}

Cart_info discoverCart(int id, MPI_Comm comm, int dim[2]) {
  Cart_info cart_info;

  cart_info.id = id;
  cart_info.comm = comm;

  MPI_Cart_coords(comm, id, 2, cart_info.coord);

  cart_info.coord[1] = cart_info.coord[1] + 1;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.up);
  cart_info.coord[1] = cart_info.coord[1] - 2;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.down);
  cart_info.coord[1] = cart_info.coord[1] + 1;

  cart_info.left = -1;
  cart_info.right = -1;

  if(cart_info.coord[0] != 0) {
    cart_info.coord[0] = cart_info.coord[0] - 1;
    MPI_Cart_rank(comm, cart_info.coord, &cart_info.left);
    cart_info.coord[0] = cart_info.coord[0] + 1;
  }

  if(cart_info.coord[0] + 1 != dim[0]) {
    cart_info.coord[0] = cart_info.coord[0] + 1;
    MPI_Cart_rank(comm, cart_info.coord, &cart_info.right);
    cart_info.coord[0] = cart_info.coord[0] - 1;
  }

  return cart_info;
}

void print_cart_info(Cart_info cart_info) {
  printf("id = %d has up = %d, down = %d, left = %d, right = %d\n", cart_info.id, cart_info.up, cart_info.down, cart_info.left, cart_info.right);
}


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
  MPI_Status status;
  MPI_Request request;
  double ***sendBuffs = NULL;

  // if(world_rank == MASTER)
  //   printf("MPI_Scatter\n");

  MPI_Irecv(&buf[0][0], mp*np, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &request);

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
      MPI_Send(&(sendBuffs[i])[0][0], mp*np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }

  }

  MPI_Wait(&request, &status);

  if(world_rank == MASTER)
    free(sendBuffs);
}

void gather_masterbuf(double **masterbuf, double **buf, int m, int n, int mp, int np, int world_rank, int world_size, MPI_Comm comm) {
  int curr_coord[2];
  int send_rank;
  MPI_Status status;
  MPI_Request request;
  double ***recvBuffs = NULL;

  // if(world_rank == MASTER)
  //   printf("MPI_Gather\n");

  MPI_Isend(&buf[0][0], mp*np, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &request);

  if(world_rank == MASTER) {
    recvBuffs = (double ***) arralloc(sizeof(double), 3, world_size, mp, np);

    for (int i = 0; i < world_size; ++i) {
      MPI_Recv(&(recvBuffs[i])[0][0], mp*np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
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

void calculate_halo_swaps(double **old, int m, int n, Cart_info cart_info) {

  MPI_Request request[4], send_request[4];
  MPI_Status status[4], send_status[4];
  double *top_send_buff, *bottom_send_buff, *left_send_buff, *right_send_buff;
  double *top_recv_buff, *bottom_recv_buff, *left_recv_buff, *right_recv_buff;

  top_send_buff = (double *) arralloc(sizeof(double), 1, m);
  bottom_send_buff = (double *) arralloc(sizeof(double), 1, m);
  left_send_buff = (double *) arralloc(sizeof(double), 1, n+2);
  right_send_buff = (double *) arralloc(sizeof(double), 1, n+2);

  top_recv_buff = (double *) arralloc(sizeof(double), 1, m);
  bottom_recv_buff = (double *) arralloc(sizeof(double), 1, m);
  left_recv_buff = (double *) arralloc(sizeof(double), 1, n+2);
  right_recv_buff = (double *) arralloc(sizeof(double), 1, n+2);

  for (int i = 1; i < m+1; ++i) {
    top_send_buff[i-1] = old[i][n];
    bottom_send_buff[i-1] = old[i][1];
  }
  for (int i = 0; i < n+2; ++i) {
    left_send_buff[i] = old[1][i];
    right_send_buff[i] = old[m][i];
  }

  if(has_top(cart_info))
    MPI_Irecv(top_recv_buff, m, MPI_DOUBLE, cart_info.up, 0, cart_info.comm, &request[0]);
  if(has_bottom(cart_info))
    MPI_Irecv(bottom_recv_buff, m, MPI_DOUBLE, cart_info.down, 1, cart_info.comm, &request[1]);
  if(has_left(cart_info))
    MPI_Irecv(left_recv_buff, n+2, MPI_DOUBLE, cart_info.left, 2, cart_info.comm, &request[2]);
  if(has_right(cart_info))
    MPI_Irecv(right_recv_buff, n+2, MPI_DOUBLE, cart_info.right, 3, cart_info.comm, &request[3]);

  if(has_top(cart_info))
    MPI_Isend(top_send_buff, m, MPI_DOUBLE, cart_info.up, 1, cart_info.comm, &send_request[0]);
  if(has_bottom(cart_info))
    MPI_Isend(bottom_send_buff, m, MPI_DOUBLE, cart_info.down, 0, cart_info.comm, &send_request[1]);
  if(has_left(cart_info))
    MPI_Isend(left_send_buff, n+2, MPI_DOUBLE, cart_info.left, 3, cart_info.comm, &send_request[2]);
  if(has_right(cart_info))
    MPI_Isend(right_send_buff, n+2, MPI_DOUBLE, cart_info.right, 2, cart_info.comm, &send_request[3]);

  if(has_top(cart_info)) {
    MPI_Wait(&request[0], &status[0]);
    MPI_Wait(&send_request[0], &send_status[0]);
  }
  if(has_bottom(cart_info)) {
    MPI_Wait(&request[1], &status[1]);
    MPI_Wait(&send_request[1], &send_status[1]);
  }
  if(has_left(cart_info)) {
    MPI_Wait(&request[2], &status[2]);
    MPI_Wait(&send_request[2], &send_status[2]);
  }
  if(has_right(cart_info)){
    MPI_Wait(&request[3], &status[3]);
    MPI_Wait(&send_request[3], &send_status[3]);
  }

  for (int i = 1; i < m+1; ++i) {
    if(has_top(cart_info))
      old[i][n+1] = top_recv_buff[i-1];
    if(has_bottom(cart_info))
      old[i][0] = bottom_recv_buff[i-1];
  }

  for (int i = 0; i < n+2; ++i) {
    if(has_left(cart_info))
      old[0][i] = left_recv_buff[i];
    if(has_right(cart_info))
      old[m+1][i] = right_recv_buff[i];
  }

  free(top_recv_buff);
  free(bottom_recv_buff);
  free(left_recv_buff);
  free(right_recv_buff);
  free(top_send_buff);
  free(bottom_send_buff);
  free(left_send_buff);
  free(right_send_buff);
}

void calculate(double **buf, double **old, double **new, double **edge, int m, int n, Cart_info cart_info) {
  int i, j, iter;

  // if(world_rank == MASTER)
  //   printf("calculate\n");

  for (iter=1;iter<=MAXITER; iter++) {
    // if(iter%PRINTFREQ==0) {
    //   printf("Iteration %d\n", iter);
    // }

    /* Implement periodic boundary conditions on bottom and top sides */
    /* Implement halo swaps */
    calculate_halo_swaps(old, m, n, cart_info);

    MPI_Barrier(cart_info.comm);

    for (i=1;i<m+1;i++) {
      for (j=1;j<n+1;j++) {
        new[i][j]=0.25*(old[i-1][j  ]
                       +old[i+1][j  ]
                       +old[i  ][j-1]
                       +old[i  ][j+1] - edge[i][j]);
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

  int world_rank, world_size;
  double start_time, end_time;
  int m, n, mp, np;

  MPI_Comm comm;
  int dim[2], period[2], reorder;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  dim[0]=sqrt(world_size);
  dim[1]=world_size/dim[0];

  // doesn't work for
  // n=3
  // dim[1]=world_size;
  // dim[0]=1;

  period[0]=0;
  period[1]=1;
  reorder=1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
  // if(world_rank == MASTER)
  // {
  //   printf("Number of iterations = %d\n", MAXITER);
  //   printf("Result topology: %d x %d\n", dim[0], dim[1]);
  // }

  Cart_info cart_info = discoverCart(world_rank, comm, dim);

  print_cart_info(cart_info);

  strcpy(filename, "resources/");
  if(argc == 2) {
    strcat(filename, argv[1]);
  }
  else {
    strcat(filename, "edgenew192x128.pgm");
  }

  pgmsize (filename, &m, &n);

  //new sizes of tables
  mp = m/dim[0];
  np = n/dim[1];

  // printf("Rank: %d (%d, %d) Processing %d x %d = %d fraction of image %d x %d = %d\n", world_rank, coord[0], coord[1], mp, np, mp*np, m, n, m*n);

  MPI_Barrier(MPI_COMM_WORLD);

  // if(world_rank == MASTER)
  // {
  //   printf("Fractions: %d x %d x %d = %d of image %d\n", world_size, mp, np, world_size*mp*np, m*n);
  // }

  if(world_rank == MASTER)
  {
    masterbuf = (double **) arralloc(sizeof(double), 2, m, n);
    // printf("\nReading <%s>\n", filename);
    pgmread(filename, &masterbuf[0][0], m, n);
  }

  allocate_tables(&buf, &old, &new, &edge, mp, np);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
    start_time = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  scatter_masterbuf(masterbuf, buf, m, n, mp, np, world_rank, world_size, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  initialize_tables(buf, old, edge, mp, np, cart_info, dim);

  MPI_Barrier(MPI_COMM_WORLD);

  calculate(buf, old, new, edge, mp, np, cart_info);

  MPI_Barrier(MPI_COMM_WORLD);

  gather_masterbuf(masterbuf, buf, m, n, mp, np, world_rank, world_size, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER) {
    end_time = MPI_Wtime();
    // printf("Running time = %f\n", end_time - start_time);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == MASTER)
  {
    strcpy(filename, "output/");
    if(argc == 2) {
      strcat(filename, argv[1]);
    }
    else {
      strcat(filename, "imagenew192x128.pgm");
    }
    pgmwrite(filename, &masterbuf[0][0], m, n);
    free(masterbuf);
  }

  free_tables(buf, old, new, edge);

  MPI_Finalize();

  return 0;
}
