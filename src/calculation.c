#include "../include/calculation.h"

//finds the average pixel
double get_average_pixel(double **table, MPI_Comm comm, int m, int n, int mp, int np) {
  double sum = 0, global_sum = 0;
  int i, j;

  for (i=1; i<mp+1; i++) {
    for (j=1; j<np+1; j++) {
      sum += table[i][j];
    }
  }
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  return global_sum/(m*n);
}

//prints the average pixel
void print_average_pixel(double **table, int iter, Cart_info cart_info, int m, int n, int mp, int np, char *filename) {
  double average_pixel;
  FILE * fp;
  char average_pixel_filename[FILENAME_SIZE], tmp_name[FILENAME_SIZE], c;
  int i=0;

  average_pixel = get_average_pixel(table, cart_info.comm, m, n, mp, np);
  if(cart_info.id == MASTER) {
    strcpy(average_pixel_filename, AVERAGE_PIXEL_FOLDER);
    
    while(filename[i] != '.') {
      tmp_name[i] = filename[i];
      i++;
    }
    tmp_name[i] = '\0';

    strcat(average_pixel_filename, tmp_name);

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

//send and receives the halo pixels from and to the neighbours
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

//calculates the borders depending on the halo pixels
void caclulate_halo(double **old, double **new, double **edge, int m, int n, Cart_info cart_info, MPI_Request request[8], MPI_Status status[8]) {
  int i, j;

  MPI_Waitall(8, request, status);

  // calculate halo for columns
  for (i=1;i<m+1;i+=m-1) {
    for (j=1;j<n+1;j++) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }

  // calculate halo for rows
  for (j=1;j<n+1;j+=n-1) {
    for (i=1;i<m+1;i++) {
      new[i][j]=0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
    }
  }
}

//calculates the max difference in the pixels
double calculate_max_diff(double **old, double **new, int m, int n) {
  double max_diff = -1, diff;
  int i, j;

  //for each pixel find the difference and if it is max it stores
  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
      diff = fabs(old[i][j] - new[i][j]);
      if(diff > max_diff)  max_diff = diff;
    }
  }

  return max_diff;
}

int can_terminate(double **old, double **new, int m, int n, MPI_Comm comm) {
  double max_diff, global_max_diff;

  //calculates the max difference in the pixels
  max_diff = calculate_max_diff(old, new, m, n);

  //reduces to all of the processes the max difference in the pixels
  MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, comm);
  if(global_max_diff < MIN_DIFF && global_max_diff >= 0)
    return 1;

  return 0;
}

double calculate(double **edge, double **old, double **new, int m, int n, int mp, int np, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, char *filename) {
  int i, j, iter;
  MPI_Request request[8];
  MPI_Status status[8];
  double iter_time, average_iter_time = 0;

  init_mpi_datatypes_row_col(mpi_Datatypes, mp, np);

  MPI_Barrier(cart_info.comm);
  if(cart_info.id == MASTER)
    iter_time = MPI_Wtime();

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

    if(PRINTFREQ != 0 && (iter%PRINTFREQ==0 || iter == 1)) {
      //prints average pixel every PRINTFREQ
      print_average_pixel(old, iter, cart_info, m, n, mp, np, filename);

      //check if can terminate every PRINTFREQ
      if(can_terminate(old, new, mp, np, cart_info.comm))
        break;
    }

    //writes the new table to the old to send it
    for (i=1;i<mp+1;i++) {
      for (j=1;j<np+1;j++) {
        old[i][j]=new[i][j];
      }
    }

  }

  MPI_Barrier(cart_info.comm);
  if(cart_info.id == MASTER)
    iter_time = MPI_Wtime() - iter_time;

  //calculates average iteration time
  average_iter_time = (float)iter_time / (float)iter;
  return average_iter_time;
}
