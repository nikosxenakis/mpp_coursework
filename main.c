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
  int dim[2] = {0, 0}, period[2] = {0, 1}, reorder = 1;
  MPI_Comm comm;
  Cart_info cart_info;
  Mpi_Datatypes mpi_Datatypes;

  initialization(&world_rank, &world_size, filename, argc, argv, &m, &n, &masterbuf);

  decomposition(world_size, m, n, dim, &mp, &np, &max_mp, &max_np);

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

  cart_info = discoverCart(world_rank, comm, world_size, dim);

  mpi_Datatypes = init_mpi_datatypes(n, mp, np, max_mp, max_np);

  initialize_sizes(cart_info, &mp, &np, max_mp, max_np);

  allocate_tables(&edge, &old, &new, mp, np);

  initialize_tables(old, mp, np, cart_info);

  scatter_masterbuf(masterbuf, edge, mp, np, cart_info, mpi_Datatypes);

  average_iter_time = calculate(edge, old, new, m, n, mp, np, cart_info, &mpi_Datatypes, &filename[12]);

  gather_masterbuf(masterbuf, old, mp, np, cart_info, mpi_Datatypes);

  finilization(world_rank, world_size, argc, argv, filename, average_iter_time, masterbuf, m, n);

  free_tables(world_rank, masterbuf, edge, old, new);

  return 0;
}
