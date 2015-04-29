#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define ARRAYSIZE 10
#define SUBARRAYSIZE 5

struct Array{
  double* value;
  int rows;
  int cols;
};

void initArray(struct Array* m, int rows, int cols) {
	m->value = (double*) malloc((size_t)(rows*cols*sizeof(double)));
	m->rows = rows;
	m->cols = cols;
	for (int i=0; i<rows*cols; i++) {
		m->value[i] = 0.0;
	}
}

int main(int argc, char **argv) {
	//Init MPI and get rank of node and size of MPI_COMM_WORLD
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Create and allocate matrix and then set values to 0.0
	struct Array *array = malloc(sizeof(struct Array));
	initArray(array, ARRAYSIZE, ARRAYSIZE);
	
	//Create MPI_Datatype
	int arraysizes[2] 		= {ARRAYSIZE,ARRAYSIZE};
	int subarraysizes[2] 	= {SUBARRAYSIZE, SUBARRAYSIZE};
	int subarraystarts[2] = {0, 0};
	MPI_Datatype type_subarray;
	MPI_Type_create_subarray(2, arraysizes, subarraysizes, subarraystarts, MPI_ORDER_C, MPI_DOUBLE, &type_subarray);
 	MPI_Type_commit(&type_subarray);
	
	if (rank == 0) {
		printf("Subarray before sending (start):\n");	
		for (int i=1;i<=array->cols/2;i++) {
			for (int j=2;j<=array->rows/2;j++) {
				array->value[i*array->rows+j] = 2;
				printf("%f ", array->value[i*array->rows+j]);
			}
			printf("\n");
		}
		printf("Subarray before sending (end):\n");	
	}
	
	//Send the subarray from node 0 to node 1
	MPI_Status status;
	if (rank == 0) {
		MPI_Send(array->value, 1, type_subarray, 1, 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(array->value, 1, type_subarray, 0, 0, MPI_COMM_WORLD, &status);
	}

	if (rank == 1) {
		printf("Entire array after sending (start):\n");	
		for (int i=0;i<array->cols;i++) {
			for (int j=0;j<array->rows;j++) {
				printf("%f ", array->value[i*array->rows+j]);
			}
			printf("\n");
		}
		printf("Entire array after sending (end):\n");	
	}
	
	//Deallocate and quit
	free(array->value);
	MPI_Finalize();
	return 0;
}
	
