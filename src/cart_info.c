#include "../include/cart_info.h"

#define MASTER 0

int has_left(Cart_info cart_info) {
  return cart_info.left != -1;
}

int has_right(Cart_info cart_info) {
  return cart_info.right != -1;
}

Cart_info discoverCart(int id, MPI_Comm comm, int dim[2]) {

  // if(id == MASTER)
  //   printf("Topology: %d x %d\n", dim[0], dim[1]);

  Cart_info cart_info;

  cart_info.id = id;
  cart_info.comm = comm;

  MPI_Cart_coords(comm, id, 2, cart_info.coord);

  cart_info.coord[1] = cart_info.coord[1] + 1;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.up);
  cart_info.coord[1] = cart_info.coord[1] - 2;
  MPI_Cart_rank(comm, cart_info.coord, &cart_info.down);
  cart_info.coord[1] = cart_info.coord[1] + 1;

  cart_info.left = -1;
  cart_info.right = -1;

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

void print_cart_info(Cart_info cart_info) {
  printf("id = %d has up = %d, down = %d, left = %d, right = %d\n", cart_info.id, cart_info.up, cart_info.down, cart_info.left, cart_info.right);
}
