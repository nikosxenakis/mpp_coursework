#ifndef MPI_TYPES
#define MPI_TYPES

#include <mpi.h>

typedef struct Mpi_Datatypes {
  MPI_Datatype cont_table;
  MPI_Datatype table;
  MPI_Datatype max_cont_table;
  MPI_Datatype max_table;
  MPI_Datatype row;
  MPI_Datatype column;
} Mpi_Datatypes;

Mpi_Datatypes init_mpi_datatypes(int n, int mp, int np, int max_mp, int max_np);

void init_mpi_datatypes_row_col(Mpi_Datatypes *mpi_Datatypes, int mp, int np);

#endif
