#ifndef CART_INFO
#define CART_INFO

#include <stdio.h>
#include <mpi.h>
#include <string.h>

#define MASTER 0

typedef struct Cart_info {
  int id;
  MPI_Comm comm;
  int world_size;
  int coord[2];
  int dim[2];
  int up;
  int down;
  int left;
  int right;
} Cart_info;

int has_left(Cart_info cart_info);

int has_right(Cart_info cart_info);

int is_top(Cart_info cart_info);

Cart_info discoverCart(int id, MPI_Comm comm, int world_size, int dim[2]);

#endif
