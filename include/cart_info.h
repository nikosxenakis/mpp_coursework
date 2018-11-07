#ifndef CART_INFO
#define CART_INFO

#include <stdio.h>
#include <mpi.h>

#define MASTER 0

typedef struct Cart_info {
  int id;
  MPI_Comm comm;
  int coord[2];
  int up;
  int down;
  int left;
  int right;
} Cart_info;

int has_left(Cart_info cart_info);

int has_right(Cart_info cart_info);

Cart_info discoverCart(int id, MPI_Comm comm, int dim[2]);

void print_cart_info(Cart_info cart_info);

#endif
