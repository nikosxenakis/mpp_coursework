#include "../include/mpi_types.h"

Mpi_Datatypes init_mpi_datatypes(int n, int mp, int np, int max_mp, int max_np) {
  Mpi_Datatypes mpi_Datatypes;

  MPI_Type_contiguous(mp*np, MPI_DOUBLE, &mpi_Datatypes.cont_table);
  MPI_Type_commit(&mpi_Datatypes.cont_table);

  MPI_Type_contiguous(max_mp*max_np, MPI_DOUBLE, &mpi_Datatypes.max_cont_table);
  MPI_Type_commit(&mpi_Datatypes.max_cont_table);

  MPI_Type_vector(mp, np, n, MPI_DOUBLE, &mpi_Datatypes.table);
  MPI_Type_commit(&mpi_Datatypes.table);

  MPI_Type_vector(max_mp, max_np, n, MPI_DOUBLE, &mpi_Datatypes.max_table);
  MPI_Type_commit(&mpi_Datatypes.max_table);

  return mpi_Datatypes;
}

void init_mpi_datatypes_row_col(Mpi_Datatypes mpi_Datatypes, int mp, int np) {

  MPI_Type_contiguous(np+2, MPI_DOUBLE, &mpi_Datatypes.row);
  MPI_Type_commit(&mpi_Datatypes.row);

  MPI_Type_vector(mp, 1, np+2, MPI_DOUBLE, &mpi_Datatypes.column);
  MPI_Type_commit(&mpi_Datatypes.column);

}
