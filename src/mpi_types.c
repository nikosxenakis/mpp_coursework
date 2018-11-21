#include "../include/mpi_types.h"

//initializes the mpi datatypes
Mpi_Datatypes init_mpi_datatypes(int n, int mp, int np, int max_mp, int max_np) {
  Mpi_Datatypes mpi_Datatypes;

  //continuous table
  MPI_Type_vector(mp, np, np+2, MPI_DOUBLE, &mpi_Datatypes.cont_table);
  MPI_Type_commit(&mpi_Datatypes.cont_table);

  //continuous table for the sized max problems
  MPI_Type_vector(max_mp, max_np, max_np+2, MPI_DOUBLE, &mpi_Datatypes.max_cont_table);
  MPI_Type_commit(&mpi_Datatypes.max_cont_table);

  //table for master
  MPI_Type_vector(mp, np, n, MPI_DOUBLE, &mpi_Datatypes.table);
  MPI_Type_commit(&mpi_Datatypes.table);

  //table for master with max sized problems
  MPI_Type_vector(max_mp, max_np, n, MPI_DOUBLE, &mpi_Datatypes.max_table);
  MPI_Type_commit(&mpi_Datatypes.max_table);

  return mpi_Datatypes;
}

//initializes the mpi datatypes about rows and columns
void init_mpi_datatypes_row_col(Mpi_Datatypes *mpi_Datatypes, int mp, int np) {

  MPI_Type_contiguous(np+2, MPI_DOUBLE, &(mpi_Datatypes->row));
  MPI_Type_commit(&(mpi_Datatypes->row));

  MPI_Type_vector(mp, 1, np+2, MPI_DOUBLE, &(mpi_Datatypes->column));
  MPI_Type_commit(&(mpi_Datatypes->column));

}
