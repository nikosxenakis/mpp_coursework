//B136013

#ifndef CALCULATION
#define CALCULATION

#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "mpi_types.h"
#include "cart_info.h"
#include "parameters.h"

#define MASTER 0
#define FILENAME_SIZE 128
#define MIN_DIFF 0.1

// performs the main calculation loop over the buffers
extern double calculate(double **edge, double **old, double **new, int m, int n, int mp, int np, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, char *filename);

#endif
