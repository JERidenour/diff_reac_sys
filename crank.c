#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>


// numerical parameters
#define N 512
#define hx 1.0/(N-1)
#define hx2 hx2*hx2
#define MAXITER 50000


int main(int argc, char *argv[])
{

	//Psq is the number of processes, must be a square
	int P, Psq, rank, rc;
	int target_north, target_south, target_east, target_west;

	// MPI
	rc = MPI_Init(&argc, &argv);
	rc = MPI_Comm_size(MPI_COMM_WORLD, &Psq);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (N < Psq) {
		fprintf(stdout, "Too few discretization points...\n");
		exit(1);
	}

	//determine position in process mesh:
	P = sqrt(Psq);	//process mesh dimension
	int north = rank > P * P - P - 1;
	int south = rank < P;
	int west = rank % P == 0;
	int east = (rank+1) % P == 0;

	target_north = (1-north)*(rank + P) + north*(rank - (P-1)*P);
	target_south = (1-south)*(rank - P) + south*(rank + (P-1)*P);
	target_east = (1-east)*(rank + 1) + east*(rank - P + 1);
	target_west = (1-west)*(rank - 1) + west*(rank + P - 1);

	//printf("process rank: %d, tn: %d, ts: %d, te: %d, tw: %d \n", rank, target_north, target_south, target_east, target_west);
	

	MPI_Finalize();
	return 0;
}