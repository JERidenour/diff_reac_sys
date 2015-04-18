#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>


// numerical parameters
#define N 64
#define hx 1.0/(N-1)
#define hx2 hx2*hx2
#define MAXITER 1


int main(int argc, char *argv[])
{

	//Psq is the number of processes, must be a square
	int P, Psq, rank, rc;
	int target_north, target_south, target_east, target_west;

	// MPI
	rc = MPI_Init(&argc, &argv);
	rc = MPI_Comm_size(MPI_COMM_WORLD, &Psq);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (N % Psq != 0) {
		fprintf(stdout, "Psq must divide N...\n");
		exit(1);
	}

	//determine position in process mesh:
	P = sqrt(Psq);	//process mesh dimension
	int north = rank > (P-1)*P-1;
	int south = rank < P;
	int west = rank % P == 0;
	int east = (rank+1) % P == 0;

	//calculate rank of neighbors
	target_north = (1-north)*(rank + P) + north*(rank - (P-1)*P);
	target_south = (1-south)*(rank - P) + south*(rank + (P-1)*P);
	target_east = (1-east)*(rank + 1) + east*(rank - P + 1);
	target_west = (1-west)*(rank - 1) + west*(rank + P - 1);

	//printf("process rank: %d, tn: %d, ts: %d, te: %d, tw: %d \n", rank, target_north, target_south, target_east, target_west);

	//the long solution vector
	int n_local = N/Psq;
	double u[n_local];
	//example values (index)
	for(int i = 0; i < n_local; i++){
		u[i] = i;
	}

	int n = sqrt(n_local); //length of boundary data
	//vectors receiving boudary data
	double uin_north[n], uin_south[n], uin_east[n], uin_west[n];
	//vectors sending boundary data
	double uout_north[n], uout_south[n], uout_east[n], uout_west[n];
	for(int i=0;i<n;i++) {
		uout_north[i] = u[i+(n-1)*n];
		uout_south[i] = u[i];
		uout_east[i] = u[(i+1)*n-1];
		uout_west[i] = u[i*n];
	}


	// are we red? (if not, we are black)
	bool red = (rank % 2) ? false : true;
	if(red){
		MPI_Send(
			uout_north,
			n,
			MPI_DOUBLE,
			target_north,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_south,
			n,
			MPI_DOUBLE,
			target_south,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_east,
			n,
			MPI_DOUBLE,
			target_east,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_west,
			n,
			MPI_DOUBLE,
			target_west,
			0,
			MPI_COMM_WORLD);

		MPI_Recv( 
			uin_north,
			n,
			MPI_DOUBLE,
			target_north,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv( 
			uin_south,
			n,
			MPI_DOUBLE,
			target_south,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv( 
			uin_east,
			n,
			MPI_DOUBLE,
			target_east,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv( 
			uin_west,
			n,
			MPI_DOUBLE,
			target_west,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
	}else{
		MPI_Send(
			uout_north,
			n,
			MPI_DOUBLE,
			target_north,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_south,
			n,
			MPI_DOUBLE,
			target_south,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_east,
			n,
			MPI_DOUBLE,
			target_east,
			0,
			MPI_COMM_WORLD);
		MPI_Send(
			uout_west,
			n,
			MPI_DOUBLE,
			target_west,
			0,
			MPI_COMM_WORLD);

		MPI_Recv(
			uin_north,
			n,
			MPI_DOUBLE,
			target_north,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv(
			uin_south,
			n,
			MPI_DOUBLE,
			target_south,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);		
		MPI_Recv(
			uin_east,
			n,
			MPI_DOUBLE,
			target_east,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv(
			uin_west,
			n,
			MPI_DOUBLE,
			target_west,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);



	}

	// printf("rank: %d, uin_north: %f, %f, %f  \n", rank, uin_north[0], uin_north[1], uin_north[2]);
	// printf("rank: %d, uin_south: %f, %f, %f  \n", rank, uin_south[0], uin_south[1], uin_south[2]);
	// printf("rank: %d, uin_east: %f, %f, %f  \n", rank, uin_east[0], uin_east[1], uin_east[2]);
	// printf("rank: %d, uin_west: %f, %f, %f  \n", rank, uin_west[0], uin_west[1], uin_west[2]);

	// for(int i = 0; i < n; i++) printf("rank: %d, uout_north: %f  \n", rank, uout_north[i]);
	// for(int i = 0; i < n; i++) printf("rank: %d, uout_south: %f  \n", rank, uout_south[i]);
	// for(int i = 0; i < n; i++) printf("rank: %d, uout_east: %f  \n", rank, uout_east[i]);
	// for(int i = 0; i < n; i++) printf("rank: %d, uout_west: %f  \n", rank, uout_west[i]);

	MPI_Finalize();
	return 0;
}