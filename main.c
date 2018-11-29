//B136013

#include "include/mpi_types.h"
#include "include/cart_info.h"
#include "include/initialization.h"
#include "include/calculation.h"
#include "include/finilization.h"

#define FILENAME_SIZE 128

int main (int argc, char** argv) {
  double **masterbuf = NULL, **edge, **old, **new;
  char filename[FILENAME_SIZE];
  int world_rank, world_size, m, n, mp, np, max_mp, max_np;
  double average_iter_time = 0;
  int dim[2] = {0, 0}, period[2] = {0, 1}, reorder = 0;
  MPI_Comm comm;
  Cart_info cart_info;
  Mpi_Datatypes mpi_Datatypes;

  //initialization processes for the program
  initialization(&world_rank, &world_size, filename, argc, argv, &m, &n, &masterbuf);

  //decomposes the problem to the given dimensions and returns the new subproblem sizes
  decomposition(world_size, m, n, dim, &mp, &np, &max_mp, &max_np);

  //create the new MPI virtual topology
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

  //discovers and create the Cart_Info structure based on the given topology
  cart_info = discoverCart(world_rank, comm, world_size, dim);

  //gets the mpi derived datatypes
  mpi_Datatypes = init_mpi_datatypes(n, mp, np, max_mp, max_np);

  //fixes the appropriate sizes for this processes' problem
  initialize_sizes(cart_info, &mp, &np, max_mp, max_np);

  //allocates the necessary 2D tables
  allocate_tables(&edge, &old, &new, mp, np);

  //initializes the allocated tables
  initialize_tables(old, mp, np, cart_info);

  //scatters the masterbuf to the processes
  scatter_masterbuf(masterbuf, edge, mp, np, cart_info, mpi_Datatypes);

  // performs the main calculation loop over the buffers and returns the average iteration time
  average_iter_time = calculate(edge, old, new, m, n, mp, np, cart_info, &mpi_Datatypes, &filename[12]);

  //gathers the processes' reconstructed parts to the masterbuf of the master process
  gather_masterbuf(masterbuf, old, mp, np, cart_info, mpi_Datatypes);

  //exports the output and finilizes the program
  finilization(world_rank, world_size, argc, argv, filename, average_iter_time, masterbuf, m, n);

  //frees the allocated tables
  free_tables(world_rank, masterbuf, edge, old, new);

  return 0;
}
