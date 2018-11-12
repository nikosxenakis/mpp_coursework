#include "../include/cart_info.h"

#define MASTER 0

int has_left(Cart_info cart_info) {
  return cart_info.left != MPI_PROC_NULL;
}

int has_right(Cart_info cart_info) {
  return cart_info.right != MPI_PROC_NULL;
}

int is_top(Cart_info cart_info) {
  return cart_info.coord[1] + 1 == cart_info.dim[1];
}

Cart_info discoverCart(int id, MPI_Comm comm, int world_size, int dim[2]) {

  Cart_info cart_info;

  cart_info.id = id;
  cart_info.comm = comm;
  cart_info.world_size = world_size;
  memcpy(cart_info.dim, dim, 2*sizeof(int));

  MPI_Cart_coords(comm, id, 2, cart_info.coord);

  cart_info.coord[1] = cart_info.coord[1] + 1;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.up);
  cart_info.coord[1] = cart_info.coord[1] - 2;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.down);
  cart_info.coord[1] = cart_info.coord[1] + 1;

  cart_info.left = MPI_PROC_NULL;
  cart_info.right = MPI_PROC_NULL;

  if(cart_info.coord[0] != 0) {
    cart_info.coord[0] = cart_info.coord[0] - 1;
    MPI_Cart_rank(comm, cart_info.coord, &cart_info.left);
    cart_info.coord[0] = cart_info.coord[0] + 1;
  }

  if(cart_info.coord[0] + 1 != dim[0]) {
    cart_info.coord[0] = cart_info.coord[0] + 1;
    MPI_Cart_rank(comm, cart_info.coord, &cart_info.right);
    cart_info.coord[0] = cart_info.coord[0] - 1;
  }

  return cart_info;
}
