#ifndef CALCULATION
#define CALCULATION

#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "mpi_types.h"
#include "cart_info.h"

#define MASTER 0
#define FILENAME_SIZE 128
#define MIN_DIFF 0.1

// performs the main calculation loop over the buffers
extern double calculate(double **edge, double **old, double **new, int m, int n, int mp, int np, Cart_info cart_info, Mpi_Datatypes *mpi_Datatypes, char *filename);

#endif

//preprocessor definitions for testing various configurations
#define PROD
// #define TIMING_TEST
// #define AVERAGE_PIXEL_TEST
// #define TIMING_WITH_INTERVALS_TEST

#ifdef PROD
	#define MAXITER   20000
	#define PRINTFREQ  100
#endif

#ifdef AVERAGE_PIXEL_TEST
	#define MAXITER   20000
	#define PRINTFREQ  1
#endif

#ifdef TIMING_TEST
	#define MAXITER   1000
	#define PRINTFREQ  0
#endif

#ifdef TIMING_WITH_INTERVALS_TEST
	#define MAXITER   1000
	#define PRINTFREQ  1
#endif
