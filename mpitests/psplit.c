#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ndirs 8	// how many neighbours
#define n_min 4	// minimum meaningful size of local subdomain

typedef enum { false, true } bool;

void printarr(int **data, int n, char *str);
double **alloc2DArray(int n);
int badGridParams(int rank, int P, int p, int N, int n);

int main(int argc, char **argv)
{

	int rank, P, p;
	int n, N;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

	// set size of global domain and grid
	if (argc > 1 && atoi(argv[1]) > P * n_min )
		N = atoi(argv[1]);	// use argument if big enough
	else
	{
		N = P * n_min;
		if (!rank)
		{
			fprintf(stdout, "Setting N to minimum %d x %d.\n", P, n_min );
			fflush(stdout);
		}
	}
	p = (int) sqrt(P);	// size of process grid
	n = (int) N / P;	// size of subdomain

	// check grid parameters
	if (badGridParams(rank, P, p, N, n))
	{
		MPI_Finalize();
		return 0;
	}
	else if (!rank)
	{
		fprintf(stdout, "P: %d, p: %d, N: %d, n: %d\n", P, p, N, n );
		fflush(stdout);
	}

	// if we get this far we are probably ok

	// side of local 2D array is N+2 with ghost rows and cols
	double **localData = alloc2DArray((n + 2) * (n + 2));

	// ghost borders and corners
	//MPI_Datatype north, east, south, west;
	//MPI_Datatype northEast, southEast, southWest, northWest;

	// array of the ghost submatrices numbered clockwise N E S W NE SE SW NW
	MPI_Datatype ghosts[ndirs];

	// list of coordinates for ghost parts
	int gCoords[] = {	0, 1,			// 	north starts
	                    1, n,			//	north size
	                    1, n + 1,		//	east starts
	                    n, 1,			// 	east size
	                    n + 1, 1,		// 	south starts
	                    1, n,			// 	south size
	                    1, 0,			//	west starts
	                    n, 1,			//  west size
	                    0, n + 1,		// 	NE starts
	                    1, 1,			//	NE size
	                    n + 1, n + 1, 	// 	SE starts
	                    1, 1,			// 	SE size
	                    n + 1, 1,		// 	SW starts
	                    1, 1,			//  SW size
	                    0, 0,			// 	NW starts
	                    1, 1			//  NW size
	                };

	// create the ghost subarrays
	for (int i = 0; i < ndirs; ++i)
	{
		int totalDim[2]  = {n + 2, n + 2};	// dim of local data array
		int starts[2] = {gCoords[i * 4], gCoords[(i * 4) + 1]};
		int subsizes[2]  = {gCoords[(i * 4) + 2], gCoords[(i * 4) + 3]};
		MPI_Type_create_subarray(2, totalDim, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ghosts[i]);
		MPI_Type_commit(&ghosts[i]);
	}

	// create mpi process mesh
	int meshCoords[2];	// my mesh coordinates
	int meshRank;		// my meshRank
	int dimensions[2] = {p, p}; // dimensions of mesh
	int wraparound[2] = {1, 1}; // wraparound

	MPI_Comm processMesh;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wraparound, 1, &processMesh);
	MPI_Comm_rank(processMesh, &meshRank);	// get my mesh rank
	MPI_Cart_coords(processMesh, meshRank, 2, meshCoords); 	// get my i,j coordinates in process grid

	//printf("Rank: %d, meshRank: %d, i,j: %d,%d\n", rank, meshRank, meshCoords[0], meshCoords[1]);
	

	// get neighbours rank within mesh communicator
	int hood[3][3];	// hood[i][j] = meshRank, hood[1][1] = meshRank
	//int nRank, eRank, sRank, wRank, neRank, seRank, swRank, nwRank;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			int xof = meshCoords[0] - 1;
			int yof = meshCoords[1] - 1;
			int crd[2] = {xof + i, yof + j};
			int rnk;
			MPI_Cart_rank(processMesh, crd, &rnk); 	// get rank
			hood[i][j] = rnk;
		}
	}

	printf("Rank: %d, meshRank: %d, i, j: %d,%d, hood: %d %d %d %d %d %d %d %d %d\n",
	       rank, meshRank, meshCoords[0], meshCoords[1],
	       hood[0][0], hood[0][1], hood[0][2],
	       hood[1][0], hood[1][1], hood[1][2],
	       hood[2][0], hood[2][1], hood[2][2]
	      );
	fflush(stdout);




	// free the ghost spirits
	for (int i = 0; i < ndirs; ++i)
		MPI_Type_free(&ghosts[i]);

	MPI_Finalize();

	// free local array
	free(localData[0]);
	free(localData);

	return 0;
}


/*	Returns nonzero if bad. Zero if OK.
*/
int badGridParams(int rank, int P, int p, int N, int n)
{
	// P must be a square
	if (p * p != P)
	{
		if (!rank)
		{
			fprintf(stderr, "P: %d, p: %d. P must be a square, exiting.\n", P, p );
		}
		return 1;
	}
	// n must be a square
	int nsq = (int) sqrt(n);
	if (nsq * nsq != n)
	{
		if (!rank)
		{
			fprintf(stderr, "n: %d. n must be a square, exiting.\n", n );
		}
		return 1;
	}
	// check that N can be evenly distributed over the processes
	if ( (n * P) != N)
	{
		if (!rank)
		{
			fprintf(stderr, "P: %d, n: %d, N: %d\n", P, n, N );
			fprintf(stderr, "Can't compute, n*P must equal N. Exiting.\n");
		}
		return 1;
	}

	// so far so good
	return 0;
}





void printarr(int **data, int n, char *str) {
	printf("-- %s --\n", str);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%3d ", data[i][j]);
		}
		printf("\n");
	}
}

double **alloc2DArray(int n) {
	double *data = malloc(n * n * sizeof(double));
	double **arr = malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
		arr[i] = &(data[i * n]);

	return arr;
}