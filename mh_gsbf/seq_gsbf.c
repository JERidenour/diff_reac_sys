#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define NDIRS 8	// how many neighbours
#define N_MIN 1	// minimum meaningful size of local subdomain

/*#define BORDER_NORTH [1][1]			// 	north starts
#define BORDER_EAST [1][n]			//	east starts
#define BORDER_SOUTH [n][1]			// 	south starts
#define BORDER_WEST [1][1]			//	west starts
#define BORDER_NORTHEAST [1][n]			// 	NE starts
#define BORDER_SOUTHEAST [n][n]		 	// 	SE starts
#define BORDER_SOUTHWEST [n][1]			// 	SW starts
#define BORDER_NORTHWEST [1][1]			// 	NW starts

#define GHOST_NORTH [0][1]
#define GHOST_EAST [1][n + 1]		//	east starts
#define GHOST_SOUTH [n + 1][1]		// 	south starts
#define GHOST_WEST [1][0]			//	west starts
#define GHOST_NORTHEAST [0][n + 1]		// 	NE starts
#define GHOST_SOUTHEAST [n + 1][n + 1] 	// 	SE starts
#define GHOST_SOUTHWEST [n + 1][0]		// 	SW starts
#define GHOST_NORTHWEST [0][0]			// 	NW starts

#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3
#define NORTHEAST 4
#define SOUTHEAST 5
#define SOUTHWEST 6
#define NORTHWEST 7
#define CENTER 8

#define VERTICAL 1
#define HORIZONTAL 2
#define DIAGONAL_UP 3
#define DIAGONAL_DOWN 4*/

#define I 0 	// matrix coordinate indexes
#define J 1

//=========================================================================
// your computation defines
//=========================================================================

// mpicc -o gsbf gsbf.c -std=c99 -Wall -O2

// run program like so:
// mpirun -n <P processes> gsbf <N size of global domain> <id>
// for example:
// mpirun -n 16 gsbf 16

#define MAXITER 50000

//#define BLOPSIZE 4

// my testcase

/*#define HT 0.25
#define DU 0.00002
#define DV 0.00001*/
/*#define DU 2.0*1e-5
#define DV 1.0*1e-5*/
/*#define F 0.026
#define K 0.0550*/

// JR's values
// #define HT 0.19
#define HT 0.019
#define DU 2.0*1e-6
#define DV 1.0*1e-6
#define F 0.0375
#define K 0.0634

#define U0 0.5		// non-zero starting value at blob
#define V0 0.25		// non-zero starting value at blob

//=========================================================================
// declarations
//=========================================================================

typedef enum { false, true } bool;

void printarr(double **data, int n, char *str);
double **alloc2DArray(int n);
int badGridParams(int rank, int P, int p, int N, int n);

int main(int argc, char **argv)
{
	// timing
	double tInit = 0;
	//double tComm = 0;
	double tCalc = 0;
	double startTimeInit, endTimeInit;
	//double startTimeComm, endTimeComm;
	double startTimeCalc, endTimeCalc;

	//=========================================================================
	// MPI and subdomain basic setup
	//=========================================================================

	//int rank, P, p;
	int n, N;

	MPI_Init(&argc, &argv);

	// ###################################
	// Start timing INITIALIZATION
	startTimeInit = MPI_Wtime();
	// ###################################

	/*
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	*/

	// set size of global domain and grid
	if (argc > 1 && atoi(argv[1]) > 16 ) {
		N = atoi(argv[1]);	// use argument if big enough
		fprintf(stdout, "Setting N and n to: %d x %d.\n", N, N );
	}

	else
	{
		N = 16;

		fprintf(stdout, "Setting N to minimum: %d x %d.\n", N, N );
		fflush(stdout);

	}
	/*p = (int) sqrt(P);	// size of process grid
	n = (int) N / p;	// size of subdomain*/

	n = N; // size of subdomain

	// check grid parameters
	/*if (badGridParams(rank, P, p, N, n))	// TODO: this function is just a hack
	{
		MPI_Finalize();
		return 0;
	}
	else if (!rank)
	{
		fprintf(stdout, "P: %d, p: %d, N: %d, n: %d\n", P, p, N, n );
		fflush(stdout);
	}*/

	// TODO: add an ID argument so we can identify the particular run
	// if no ID given, create an ID based on system time (timestamp)
	// suggestion concatenate this ID to all filenames so we get
	// filename_ID.txt


	//=========================================================================
	// MPI cartesian process mesh
	//=========================================================================
	/*
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

	// What color am I
	// even row AND even process = red, even row AND odd process = black
	// odd row AND even process = black, odd row AND odd prcocess = red
	int red;
	if ( (meshCoords[I] % 2) == 0 )	// even row
	{
		red = (meshRank % 2) ? 0 : 1;	// even process = red, else black
	}
	else	// odd row
	{
		red = (meshRank % 2) ? 1 : 0;	// even process = black, else red
	}
	// printf("MeshRank %d with coords [%d %d] is red? %d\n", meshRank, meshCoords[0], meshCoords[1], red);

	// gfigure out my neighbourgood
	int kernel[3][3];		// hood[i][j] = meshRank, hood[1][1] = my meshRank
	int neighbour[NDIRS];	// neighbour[DIRECTION] = rank of neighbour
	//int nRank, eRank, sRank, wRank, neRank, seRank, swRank, nwRank;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			int xof = meshCoords[I] - 1;
			int yof = meshCoords[J] - 1;
			int crd[2] = {xof + i, yof + j};
			int rnk;
			MPI_Cart_rank(processMesh, crd, &rnk); 	// get rank
			kernel[i][j] = rnk;
		}
	}
	neighbour[NORTH] = kernel[0][1];
	neighbour[EAST] = kernel[1][2];
	neighbour[SOUTH] = kernel[2][1];
	neighbour[WEST] = kernel[1][0];
	neighbour[NORTHEAST] = kernel[0][2];
	neighbour[SOUTHEAST] = kernel[2][2];
	neighbour[SOUTHWEST] = kernel[2][0];
	neighbour[NORTHWEST] = kernel[0][0];

	printf("MeshRank %d, red: %d with coords [%d %d] has hood N E S W NE SE SW NW [%d %d %d %d %d %d %d %d].\n",
	       meshRank,
	       red,
	       meshCoords[I],
	       meshCoords[J],
	       neighbour[NORTH],
	       neighbour[EAST],
	       neighbour[SOUTH],
	       neighbour[WEST],
	       neighbour[NORTHEAST],
	       neighbour[SOUTHEAST],
	       neighbour[SOUTHWEST],
	       neighbour[NORTHWEST]
	      );

	//=========================================================================
	// Submatrices for border(outgoing) and ghost(incoming)
	//=========================================================================

	// ghost border and corners
	//MPI_Datatype north, east, south, west;
	//MPI_Datatype northEast, southEast, southWest, northWest;

	// array of the ghost(incoming) submatrices numbered clockwise N E S W NE SE SW NW
	MPI_Datatype ghost[NDIRS];

	// list of coordinates for ghost parts
	int ghostSizes[] = {	1, n,			//	north size
	                        n, 1,			// 	east size
	                        1, n,			// 	south size
	                        n, 1,			//  west size
	                        1, 1,			//	NE size
	                        1, 1,			// 	SE size
	                        1, 1,			//  SW size
	                        1, 1			//  NW size
	                   };
	// create the border subarrays
	for (int i = 0; i < NDIRS; ++i)
	{
		int totalDim[2]  = {n + 2, n + 2};	// dim of local data array
		int starts[2] = {0, 0};
		int subsizes[2]  = {ghostSizes[(i * 2)], ghostSizes[(i * 2) + 1]};
		MPI_Type_create_subarray(2, totalDim, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ghost[i]);
		MPI_Type_commit(&ghost[i]);
	}

	// array of the border(outgoing) submatrices numbered clockwise N E S W NE SE SW NW
	MPI_Datatype border[NDIRS];

	// list of coordinates for border parts
	int borderSizes[] = {	1, n,			//	north size
	                        n, 1,			// 	east size
	                        1, n,			// 	south size
	                        n, 1,			//  west size
	                        1, 1,			//	NE size
	                        1, 1,			// 	SE size
	                        1, 1,			//  SW size
	                        1, 1			//  NW size
	                    };

	// create the border subarrays
	for (int i = 0; i < NDIRS; ++i)
	{
		int sizes[2]  = {n + 2, n + 2};	// dim of local data array
		int starts[2] = {0, 0};
		int subsizes[2]  = {borderSizes[(i * 2)], borderSizes[(i * 2) + 1]};
		MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &border[i]);
		MPI_Type_commit(&border[i]);
	}
	*/
	//=========================================================================
	// Computation
	//=========================================================================

	double hx = 1.0 / (N - 1);
	double ht = HT;
	double ru = (ht * DU) / (hx * hx);
	double rv = (ht * DV) / (hx * hx);

	/*if (meshRank == 0)
	{
		printf("hx: %.15f, ht: %.15f, ru: %.15f, rv: %.15f\n", hx, ht, ru, rv );
	}*/


	// allocate data vectors
	double **u = alloc2DArray((n) * (n));
	double **unew = alloc2DArray((n) * (n));

	double **v = alloc2DArray((n) * (n));
	double **vnew = alloc2DArray((n) * (n));

	// init
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			u[i][j] = 1.0;
			unew[i][j] = 1.0;
			v[i][j] = 0.0;
			vnew[i][j] = 0.0;
		}
	}

	// For one process:
	int c = floor(n / 2);
	int cc = floor(n / 4);

	for (int i = cc ; i < (c + cc); i++) {
		for (int j = cc ; j < (c + cc); j++) {
			u[i][j] = U0;
			v[i][j] = V0;
		}
	}



// For four processes:
	/*if (Psq == 4) {

		int r = 50;
		if (rank == 0) {
			for (int i = n_inner - r; i < n_inner; i++) {
				for (int j = n_inner - r; j < n_inner; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 1) {
			for (int i = 0; i < r; i++) {
				for (int j = n_inner - r; j < n_inner; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 2) {
			for (int i = n_inner - r; i < n_inner; i++) {
				for (int j = 0; j < r; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 3) {
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < r; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
	}*/

// For sixteen processes:
	/*if (P == 16) {

		if(meshRank == 0) {
			printf("Creating initial condition for 16 processes.\n");
		}

		// int r = 50;
		if (meshRank == 5 || meshRank == 6 || meshRank == 9 || meshRank == 10) {
			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= n; j++) {
					u[i][j] = 0.5;
					v[i][j] = 0.25;
				}
			}
		}
	}*/

// For sixty-four processes
	/*if (Psq == 64) {

		int r = 50;
		if (rank == 27) {
			for (int i = n_inner - r; i < n_inner; i++) {
				for (int j = n_inner - r; j < n_inner; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 28) {
			for (int i = 0; i < r; i++) {
				for (int j = n_inner - r; j < n_inner; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 35) {
			for (int i = n_inner - r; i < n_inner; i++) {
				for (int j = 0; j < r; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
		if (rank == 36) {
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < r; j++) {
					u[i + j * n_inner] = 0.5;
					v[i + j * n_inner] = 0.25;
				}
			}
		}
	}*/

	/*// inital blobs
	if (meshRank == 0) {
		for (int i = 1; i <= BLOPSIZE; ++i)
		{
			for (int j = 1; j <= BLOPSIZE; ++j)
			{
				u[i][j] = 0.5;
				v[i][j] = 0.25;
				unew[i][j] = 0.5;
				vnew[i][j] = 0.25;
			}
		}
	}
	if (meshRank == P - 1) {
		for (int i = 1; i <= BLOPSIZE; ++i)
		{
			for (int j = 1; j <= BLOPSIZE; ++j)
			{
				u[i][j] = 0.5;
				v[i][j] = 0.25;
				unew[i][j] = 0.5;
				vnew[i][j] = 0.25;
			}
		}
	}*/

// ###################################
// END timing INITIALIZATION
	endTimeInit = MPI_Wtime();
	tInit = endTimeInit - startTimeInit;
// ###################################

	for (int iter = 0; iter < MAXITER; ++iter)	// begin computation iterations
	{
		/*
		// ###################################
		// START timing COMMUNICATION
		startTimeComm = MPI_Wtime();
		// ###################################

		// exchange data
		int err = 0;
		MPI_Status status;

		if (red)
		{
			// talk to north
			err = MPI_Sendrecv(	&(u BORDER_NORTH),
			                    1,
			                    border[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    &(u GHOST_NORTH),
			                    1,
			                    ghost[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_NORTH),
			                    1,
			                    border[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    &(v GHOST_NORTH),
			                    1,
			                    ghost[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );

			// talk to south
			err = MPI_Sendrecv(	&(u BORDER_SOUTH),
			                    1,
			                    border[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    &(u GHOST_SOUTH),
			                    1,
			                    ghost[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_SOUTH),
			                    1,
			                    border[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    &(v GHOST_SOUTH),
			                    1,
			                    ghost[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );

			// talk to east
			err = MPI_Sendrecv(	&(u BORDER_EAST),
			                    n,
			                    border[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    &(u GHOST_EAST),
			                    n,
			                    ghost[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_EAST),
			                    n,
			                    border[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    &(v GHOST_EAST),
			                    n,
			                    ghost[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );

			// talk to west
			err = MPI_Sendrecv(	&(u BORDER_WEST),
			                    n,
			                    border[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    &(u GHOST_WEST),
			                    n,
			                    ghost[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_WEST),
			                    n,
			                    border[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    &(v GHOST_WEST),
			                    n,
			                    ghost[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );

		}
		else	// black
		{
			// talk to south
			err = MPI_Sendrecv(	&(u BORDER_SOUTH),
			                    1,
			                    border[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    &(u GHOST_SOUTH),
			                    1,
			                    ghost[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_SOUTH),
			                    1,
			                    border[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    &(v GHOST_SOUTH),
			                    1,
			                    ghost[SOUTH],
			                    neighbour[SOUTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );

			// talk to north
			err = MPI_Sendrecv(	&(u BORDER_NORTH),
			                    1,
			                    border[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    &(u GHOST_NORTH),
			                    1,
			                    ghost[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_NORTH),
			                    1,
			                    border[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    &(v GHOST_NORTH),
			                    1,
			                    ghost[NORTH],
			                    neighbour[NORTH],
			                    VERTICAL,
			                    processMesh,
			                    &status
			                  );

			// talk to west
			err = MPI_Sendrecv(	&(u BORDER_WEST),
			                    n,
			                    border[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    &(u GHOST_WEST),
			                    n,
			                    ghost[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_WEST),
			                    n,
			                    border[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    &(v GHOST_WEST),
			                    n,
			                    ghost[WEST],
			                    neighbour[WEST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );

			// talk to east
			err = MPI_Sendrecv(	&(u BORDER_EAST),
			                    n,
			                    border[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    &(u GHOST_EAST),
			                    n,
			                    ghost[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );
			err = MPI_Sendrecv(	&(v BORDER_EAST),
			                    n,
			                    border[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    &(v GHOST_EAST),
			                    n,
			                    ghost[EAST],
			                    neighbour[EAST],
			                    HORIZONTAL,
			                    processMesh,
			                    &status
			                  );

		} // end data exchange

		// ###################################
		// END timing COMMUNICATION
		endTimeComm = MPI_Wtime();
		tComm += (endTimeComm - startTimeComm);	// add to total
		// ###################################
		*/
		// ###################################
		// START timing CALCULATION
		startTimeCalc = MPI_Wtime();
		// ###################################

		/*int firstelem = 0;
		int lastelem = n - 1;*/
		// compute inner
		for (int i = 1; i <= n - 2; ++i)
		{
			for (int j = 1; j <= n - 2; ++j)
			{

				unew[i][j] = 	u[i][j] +
				                ru * ( u[i - 1][j] - (4 * u[i][j]) + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] ) +
				                ht * ( -u[i][j] * v[i][j] * v[i][j] + F * (1.0 - u[i][j])	 );

				vnew[i][j] = 	v[i][j] +
				                rv * ( v[i - 1][j] - (4 * v[i][j]) + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] ) +
				                ht * ( u[i][j] * v[i][j] * v[i][j] - (F + K) * v[i][j] );
			}
		}

		// compute non-corner top row: i = 0, i-1 = n-1
		int topRow = 0;
		int topRowWrap = n - 1;
		for (int j = 1; j <= n - 2; ++j)
		{

			unew[topRow][j] = 	u[topRow][j] +
			                    ru * ( u[topRowWrap][j] - (4 * u[topRow][j]) + u[topRow + 1][j] + u[topRow][j - 1] + u[topRow][j + 1] ) +
			                    ht * ( -u[topRow][j] * v[topRow][j] * v[topRow][j] + F * (1.0 - u[topRow][j])	 );

			vnew[topRow][j] = 	v[topRow][j] +
			                    rv * ( v[topRowWrap][j] - (4 * v[topRow][j]) + v[topRow + 1][j] + v[topRow][j - 1] + v[topRow][j + 1] ) +
			                    ht * ( u[topRow][j] * v[topRow][j] * v[topRow][j] - (F + K) * v[topRow][j] );
		}

		// compute non-corner bottom row: i = n-1, i+1 = 0
		int botRow = n - 1;
		int botRowWrap = 0;
		for (int j = 1; j <= n - 2; ++j)
		{
			unew[botRow][j] = 	u[botRow][j] +
			                    ru * ( u[botRow - 1][j] - (4 * u[botRow][j]) + u[botRowWrap][j] + u[botRow][j - 1] + u[botRow][j + 1] ) +
			                    ht * ( -u[botRow][j] * v[botRow][j] * v[botRow][j] + F * (1.0 - u[botRow][j])	 );

			vnew[botRow][j] = 	v[botRow][j] +
			                    rv * ( v[botRow - 1][j] - (4 * v[botRow][j]) + v[botRowWrap][j] + v[botRow][j - 1] + v[botRow][j + 1] ) +
			                    ht * ( u[botRow][j] * v[botRow][j] * v[botRow][j] - (F + K) * v[botRow][j] );
		}

		// compute non-corner left column: j = 0, j-1 = n-1
		int leftCol = 0;
		int leftColWrap = n - 1;
		for (int i = 1; i <= n - 2; ++i)
		{

			unew[i][leftCol] = 	u[i][leftCol] +
			                    ru * ( u[i - 1][leftCol] - (4 * u[i][leftCol]) + u[i + 1][leftCol] + u[i][leftColWrap] + u[i][leftCol + 1] ) +
			                    ht * ( -u[i][leftCol] * v[i][leftCol] * v[i][leftCol] + F * (1.0 - u[i][leftCol])	 );

			vnew[i][leftCol] = 	v[i][leftCol] +
			                    rv * ( v[i - 1][leftCol] - (4 * v[i][leftCol]) + v[i + 1][leftCol] + v[i][leftColWrap] + v[i][leftCol + 1] ) +
			                    ht * ( u[i][leftCol] * v[i][leftCol] * v[i][leftCol] - (F + K) * v[i][leftCol] );
		}

		// compute non-corner right column: j = n-1, j+1 = 0
		int rightCol = n - 1;
		int rightColWrap = 0;
		for (int i = 1; i <= n - 2; ++i)
		{

			unew[i][rightCol] = 	u[i][rightCol] +
			                        ru * ( u[i - 1][rightCol] - (4 * u[i][rightCol]) + u[i + 1][rightCol] + u[i][rightCol - 1] + u[i][rightColWrap] ) +
			                        ht * ( -u[i][rightCol] * v[i][rightCol] * v[i][rightCol] + F * (1.0 - u[i][rightCol])	 );

			vnew[i][rightCol] = 	v[i][rightCol] +
			                        rv * ( v[i - 1][rightCol] - (4 * v[i][rightCol]) + v[i + 1][rightCol] + v[i][rightCol - 1] + v[i][rightColWrap] ) +
			                        ht * ( u[i][rightCol] * v[i][rightCol] * v[i][rightCol] - (F + K) * v[i][rightCol] );
		}

		// compute upper left corner
		// i = 0, j = 0, i-1 = n-1, j-1 = n-1
		unew[0][0] = 	u[0][0] +
		                ru * ( u[n - 1][0] - (4 * u[0][0]) + u[0 + 1][0] + u[0][n - 1] + u[0][0 + 1] ) +
		                ht * ( -u[0][0] * v[0][0] * v[0][0] + F * (1.0 - u[0][0])	 );

		vnew[0][0] = 	v[0][0] +
		                rv * ( v[n - 1][0] - (4 * v[0][0]) + v[0 + 1][0] + v[0][n - 1] + v[0][0 + 1] ) +
		                ht * ( u[0][0] * v[0][0] * v[0][0] - (F + K) * v[0][0] );

		// compute upper right corner
		// i = 0, j = n-1, i-1 = n-1, j+1 = 0
		unew[0][rightCol] = 	u[0][rightCol] +
		                        ru * ( u[n - 1][rightCol] - (4 * u[0][rightCol]) + u[0 + 1][rightCol] + u[0][rightCol - 1] + u[0][0] ) +
		                        ht * ( -u[0][rightCol] * v[0][rightCol] * v[0][rightCol] + F * (1.0 - u[0][rightCol])	 );

		vnew[0][rightCol] = 	v[0][rightCol] +
		                        rv * ( v[n - 1][rightCol] - (4 * v[0][rightCol]) + v[0 + 1][rightCol] + v[0][rightCol - 1] + v[0][0] ) +
		                        ht * ( u[0][rightCol] * v[0][rightCol] * v[0][rightCol] - (F + K) * v[0][rightCol] );

		// compute lower right corner
		// i = n-1, i+1 = 0, j = n-1, j+1 = 0
		unew[botRow][rightCol] = 	u[botRow][rightCol] +
		                            ru * ( u[botRow - 1][rightCol] - (4 * u[botRow][rightCol]) + u[0][rightCol] + u[botRow][rightCol - 1] + u[botRow][0] ) +
		                            ht * ( -u[botRow][rightCol] * v[botRow][rightCol] * v[botRow][rightCol] + F * (1.0 - u[botRow][rightCol])	 );

		vnew[botRow][rightCol] = 	v[botRow][rightCol] +
		                            rv * ( v[botRow - 1][rightCol] - (4 * v[botRow][rightCol]) + v[0][rightCol] + v[botRow][rightCol - 1] + v[botRow][0] ) +
		                            ht * ( u[botRow][rightCol] * v[botRow][rightCol] * v[botRow][rightCol] - (F + K) * v[botRow][rightCol] );

		// compute lower left corner
		// i = n-1, i+1 = 0, j = 0, j-1 = n-1
		unew[botRow][0] = 	u[botRow][0] +
		                    ru * ( u[botRow - 1][0] - (4 * u[botRow][0]) + u[0][0] + u[botRow][0 - 1] + u[botRow][n - 1] ) +
		                    ht * ( -u[botRow][0] * v[botRow][0] * v[botRow][0] + F * (1.0 - u[botRow][0])	 );

		vnew[botRow][0] = 	v[botRow][0] +
		                    rv * ( v[botRow - 1][0] - (4 * v[botRow][0]) + v[0][0] + v[botRow][0 - 1] + v[botRow][n - 1] ) +
		                    ht * ( u[botRow][0] * v[botRow][0] * v[botRow][0] - (F + K) * v[botRow][0] );

		// swap pointers
		double **tmpnew;

		tmpnew = unew;
		unew = u;
		u = tmpnew;	// u = unew

		tmpnew = vnew;
		vnew = v;
		v = tmpnew;	// u = unew

		// ###################################
		// END timing CALCULATION
		endTimeCalc = MPI_Wtime();
		tCalc += (endTimeCalc - startTimeCalc);
		// ###################################



	} // end computation iteration


	printf("Sequential t_init: %f, t_comm: %f, t_calc: %f.\n", tInit, 0.0, tCalc);

//=========================================================================
// Output to file
//=========================================================================

	FILE *fp;	// TODO: concatenate all filenames with a unique run ID

	size_t firstelem = 0;	// we use all domain, no ghosts
	size_t lastelem = n - 1;

	// write u
	fp = fopen("u_seq.txt", "w");	// open new file to write
	for (int i = firstelem; i <= lastelem; ++i) {
		for (int j = firstelem; j <= lastelem; ++j)
		{
			fprintf(fp, "%.15f ", u[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	// write v
	fp = fopen("v_seq.txt", "w");	// open new file to write
	for (int i = firstelem; i <= lastelem; ++i) {
		for (int j = firstelem; j <= lastelem; ++j)
		{
			fprintf(fp, "%.15f ", v[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("Seq wrote [%d %d] elements.\n", n, n);

	// TODO: write timing info

//=========================================================================
// Cleanup
//=========================================================================

// TODO: right now there is no proper cleanup, fix this

// free the ghost spirits
	/*for (int i = 0; i < NDIRS; ++i)
		MPI_Type_free(&ghost[i]);*/

	//MPI_Finalize();

// free local array
	/*free(u[0]);
	free(u);*/

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
	// N must be a square
	/*int nsq = (int) sqrt(N);
	if (nsq * nsq != N)
	{
		if (!rank)
		{
			fprintf(stderr, "n: %d. n must be a square, exiting.\n", n );
		}
		return 1;
	}*/
	// check that N can be evenly distributed over the processes
	if ( (n * p) != N)
	{
		if (!rank)
		{
			fprintf(stderr, "P: %d, n: %d, N: %d\n", P, n, N );
			fprintf(stderr, "Can't compute, n*p must equal N. Exiting.\n");
		}
		return 1;
	}

	// so far so good
	return 0;
}

void printarr(double **data, int n, char *str) {
	printf("-- %s --\n", str);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// printf("%.3f ", data[i][j]);
			printf("%4d ", (int) data[i][j]);
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