#ifndef INITIALIZATION
#define INITIALIZATION

#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "arralloc.h"
#include "pgmio.h"
#include "cart_info.h"
#include "mpi_types.h"

#define MASTER 0

//calculates the boundary value
static double boundaryval(int i, int m);

//decomposes the problem to the given dimensions
extern void decomposition(int world_size, int m, int n, int dim[2], int *mp, int *np, int *max_mp, int *max_np);

//fixes the appropriate sizes
extern void initialize_sizes(Cart_info cart_info, int *mp, int *np, int max_mp, int max_np);

//allocates the necessary tables
extern void allocate_tables(double ***edge, double ***old, double ***new, int m, int n);

//initializes the allocated tables
extern void initialize_tables(double **old, int m, int n, Cart_info cart_info);

//include initialization processes for the program
extern void initialization(int *world_rank, int *world_size, char *filename, int argc, char **argv, int *m, int *n, double ***masterbuf);

//scatters the masterbuf to the workers
extern void scatter_masterbuf(double **masterbuf, double **edge, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes);

#endif
