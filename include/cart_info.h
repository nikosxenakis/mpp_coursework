#ifndef CART_INFO
#define CART_INFO

#include <stdio.h>
#include <mpi.h>
#include <string.h>

#define MASTER 0

#define TOP_TO_BOTTOM 0
#define BOTTOM_TO_TOP 1
#define LEFT_TO_RIGHT 3
#define RIGHT_TO_LEFT 4

//contains information about the topology
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

// returns true if the process has a left process
int has_left(Cart_info cart_info);

// returns true if the process has a right process
int has_right(Cart_info cart_info);

// returns true if the process has a top process
int is_top(Cart_info cart_info);

// returns true if the process has a bottom process
int is_bottom(Cart_info cart_info);

// discovers and create the Cart_Info structure based on the given topology
Cart_info discoverCart(int id, MPI_Comm comm, int world_size, int dim[2]);

#endif
