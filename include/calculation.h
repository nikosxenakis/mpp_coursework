#ifndef CALCULATION
#define CALCULATION

#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "mpi_types.h"
#include "cart_info.h"

// #define MAXITER   1000
#define MAXITER   1000000
#define MASTER 0
#define PRINTFREQ  100
#define FILENAME_SIZE 128
#define MIN_DIFF 0.1

extern double calculate(double **edge, double **old, double **new, int m, int n, int mp, int np, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, char *filename);

#endif