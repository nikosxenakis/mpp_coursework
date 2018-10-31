#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "include/pgmio.h"
#include "include/calculate.h"
#define ITERATIONS 4000
#define MASTER 0
#define P 4

void log_input(char *input_file, int *nx, int *ny) {
	pgmsize (input_file, nx, ny);
	printf("input file = %s\nPGM Size = %d x %d\n", input_file, *nx, *ny);
}

void read_buffer(double** buff, char* input_file, int nx, int ny) {
	size_t size_of_buff = ((ny*nx)*sizeof(double));
	*buff = (double*)(malloc(size_of_buff));
	pgmread (input_file, *buff, nx, ny);
}

void init_halo(double **t, int nx, int ny) {
	int i;

	for (i = 0; i < ny+2; ++i) {
		(*t)[i] = 255.0f;
		(*t)[(nx+1)*(ny+2)+i] = 255.0f;
	}

	for (i = 1; i < nx+1; ++i) {
		(*t)[(ny+2)*i] = 255.0f;
		(*t)[(ny+2)*i+(ny+1)] = 255.0f;
	}
}

void allocate_tables(double **edge, double **old, double **new, int nx, int ny) {
	size_t size_of_table = ((ny+2)*(nx+2)*sizeof(double));
	*edge = (double*)(malloc(size_of_table));
	*old = (double*)(malloc(size_of_table));
	*new = (double*)(malloc(size_of_table));
}

void send_subtable(double **table, int nx, int ny) {
	int i;
	int table_size = (nx*ny)/P;

	for(i=1; i<P; ++i) {
		double *table_start = &((*table)[i*table_size]);
		MPI_Send(table_start, table_size, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD);
	}
}

void recv_subtable(double **table, int nx, int ny) {
	int table_size = (nx*ny)/P;
	MPI_Recv(*table, table_size, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void init_tables(double **buff, double **edge, double **old, double **new, int nx, int ny) {
	int i,j;

	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			(*edge)[(i+1)*(ny+2)+(j+1)] = (*buff)[i*ny+j];
		}
	}

	init_halo(edge, nx, ny);
	init_halo(new, nx, ny);

	for (i = 0; i < nx+2; ++i) {
		for (j = 0; j < ny+2; ++j) {
			(*old)[i*(ny+2)+j] = 255.0f;
		}
	}
}

void calculate(double **edge, double **old, double **new, int nx, int ny) {
	int it,i,j;

	for(it=0; it<ITERATIONS; ++it) {
		for(i=1; i<nx+1; ++i) {
			for(j=1; j<ny+1; ++j) {
				(*new)[i*(ny+2)+j] = (
					(*old)[(i-1)*(ny+2)+j] +
					(*old)[(i+1)*(ny+2)+j] +
					(*old)[i*(ny+2)+(j-1)] +
					(*old)[i*(ny+2)+(j+1)] -
					(*edge)[i*(ny+2)+j]
				) / 4;
			}
		}

		for(i=1; i<nx+1; ++i) {
			for(j=1; j<ny+1; ++j) {
				(*old)[i*(ny+2)+j] = (*new)[i*(ny+2)+j];
			}
		}
	}
}

void export_buffer(char *output_file, double **buff, double **old, int nx, int ny) {
	int i,j;

	for(i=1; i<nx+1; ++i) {
		for(j=1; j<ny+1; ++j) {
			(*buff)[(i-1)*ny+(j-1)] = (*old)[i*(ny+2)+j];
		}
	}

	pgmwrite(output_file, *buff, nx, ny);
}

void free_master_tables(double **buff) {
	free (*buff);
}

void free_tables(double **edge, double **old, double **new) {
	free (*edge);
	free (*old);
	free (*new);
}

int main(int argc, char** argv) {

	char *input_file, *output_file = "output/edge192x128.pgm";
	int nx, ny;
	double *buff, *edge, *old, *new;
	int world_size, world_rank;

	int i;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if(argc == 2) {
		input_file = argv[1];
	}
	else {
		input_file = "resources/edge192x128.pgm";
	}

	log_input(input_file, &nx, &ny);

	if(world_rank == MASTER) {
		read_buffer(&buff, input_file, nx, ny);
	}

	allocate_tables(&edge, &old, &new, nx, ny);

	//send old and edges to others
	//master sends the subtables
	if(world_rank == MASTER) {
		send_subtable(&buff, nx, ny);
	}
	else {
		// not OLD
		recv_subtable(&old, nx, ny);
	}

	if(world_rank == MASTER) {
		init_tables(&buff, &edge, &old, &new, nx, ny);
		calculate(&edge, &old, &new, nx, ny);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//master receives the subtables

	if(world_rank == MASTER) {
		export_buffer(output_file, &buff, &old, nx, ny);
		free_master_tables(&buff);
	}

	free_tables(&edge, &old, &new);

	MPI_Finalize();

	return 0;
}
