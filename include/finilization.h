#ifndef FINILIZATION
#define FINILIZATION

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "pgmio.h"
#include "mpi_types.h"
#include "cart_info.h"

#define MASTER 0

void gather_masterbuf(double **masterbuf, double **old, int mp, int np, Cart_info cart_info, Mpi_Datatypes mpi_Datatypes);

void finilization(int world_rank, int world_size, int argc, char **argv, char *filename, double average_iter_time, double **masterbuf, int m, int n);

void free_tables(int world_rank, double **masterbuf, double **edge, double **old, double **new);

#endif
