#include "../include/finilization.h"

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

void finilization(int world_rank, int world_size, int argc, char **argv, char *filename, double average_iter_time, double **masterbuf, int m, int n) {
  FILE * fp;

  if(world_rank == MASTER) {
    fp = fopen ("./data/results.tsv", "a");
    fprintf(fp, "%s\t%d\t%f\n", filename, world_size, average_iter_time*1000.0);
    fclose(fp);

    if(argc == 2) {
      strcpy(filename, "./output/image");
      strcat(filename, &(argv[1][16]));
    }
    else {
      strcpy(filename, "./output/imagenew192x128.pgm");
    }
    pgmwrite(filename, &masterbuf[0][0], m, n);
  }

  MPI_Finalize();
}

void free_tables(int world_rank, double **masterbuf, double **edge, double **old, double **new) {
  if(world_rank == MASTER) {
    free(masterbuf);
  }
  free(edge);
  free(old);
  free(new);
}
