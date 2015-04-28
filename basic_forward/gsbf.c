#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define NDIRS 8	// how many neighbours
#define N_MIN 1	// minimum meaningful size of local subdomain

#define BORDER_NORTH [1][1]			// 	north starts
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
#define DIAGONAL_DOWN 4

#define I 0
#define J 1

//=========================================================================
// your computation defines
//=========================================================================

// run program like so:
// mpirun -n <P processes> psplit <N size of global domain>
// for example:
// mpirun -n 16 psplit 16

#define MAXITER 10

#define BLOPSIZE 4
#define HT 0.25

#define DU 0.00002
#define DV 0.00001
#define F 0.026
#define K 0.0550

//=========================================================================
// declarations
//=========================================================================

typedef enum { false, true } bool;

void printarr(double **data, int n, char *str);
double **alloc2DArray(int n);
int badGridParams(int rank, int P, int p, int N, int n);

int main(int argc, char **argv)
{

	//=========================================================================
	// MPI and subdomain basic setup
	//=========================================================================

	int rank, P, p;
	int n, N;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

	// set size of global domain and grid
	if (argc > 1 && atoi(argv[1]) > P * N_MIN )
		N = atoi(argv[1]);	// use argument if big enough
	else
	{
		N = P * N_MIN;
		if (!rank)
		{
			fprintf(stdout, "Setting N to minimum: %d x %d.\n", P, N_MIN );
			fflush(stdout);
		}
	}
	p = (int) sqrt(P);	// size of process grid
	n = (int) N / p;	// size of subdomain

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


	//=========================================================================
	// MPI cartesian process mesh
	//=========================================================================

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

	//=========================================================================
	// Computation
	//=========================================================================

	double hx = 1.0 / (N - 1);
	double ht = HT;
	double ru = (ht * DU) / (hx * hx);
	double rv = (ht * DV) / (hx * hx);

	if (meshRank == 0)
	{
		printf("hx: %.15f, ht: %.15f, ru: %.15f, rv: %.15f\n", hx, ht, ru, rv );
	}


	// allocate data vectors
	double **u = alloc2DArray((n + 2) * (n + 2));
	double **unew = alloc2DArray((n + 2) * (n + 2));

	double **v = alloc2DArray((n + 2) * (n + 2));
	double **vnew = alloc2DArray((n + 2) * (n + 2));

	// init
	for (int i = 0; i < n + 2; ++i)
	{
		for (int j = 0; j < n + 2; ++j)
		{
			u[i][j] = 1.0;
			unew[i][j] = 1.0;
			v[i][j] = 0.0;
			vnew[i][j] = .0;
		}
	}

	// inital blobs
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
	}



	for (int iter = 0; iter < MAXITER; ++iter)	// begin computation iterations
	{


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
			// MPI_Get_count(&status, ghost[NORTH], &rcount);
			// printf("%d received: %d from NORTH.\n", meshRank, rcount);

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
			/*MPI_Get_count(&status, ghost[SOUTH], &rcount);
			printf("%d received: %d from SOUTH.\n", meshRank, rcount);*/

		} // end data exchange



		// compute
		for (int i = 1; i <= n; ++i)
		{
			for (int j = 0; j <= n; ++j)
			{
				unew[i][j] = 	u[i][j] +
				                ru * ( u[i - 1][j] - (4 * u[i][j]) + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] ) +
				                ht * ( -u[i][j] * v[i][j] * v[i][j] + F * (1.0 - u[i][j])	 );

				vnew[i][j] = 	v[i][j] +
				                rv * ( v[i - 1][j] - (4 * v[i][j]) + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] ) +
				                ht * ( u[i][j] * v[i][j] * v[i][j] - (F + K) * v[i][j] );
			}
		}

		/*for (int i = 1; i <= n; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				u[i][j] = unew[i][j];
				v[i][j] = vnew[i][j];
			}
		}*/

		// swap pointers
		double **tmpnew;

		tmpnew = unew;
		unew = u;
		u = tmpnew;	// u = unew

		tmpnew = vnew;
		vnew = v;
		v = tmpnew;	// u = unew



	} // end computation iteration

	//=========================================================================
	// Output to file
	//=========================================================================

	FILE *fp;

	if (meshRank == 0) {

		size_t firstelem = 1;
		size_t lastelem = n;

		// write u
		fp = fopen("u.txt", "w");	// open new file to write
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%.15f ", u[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		// write v
		fp = fopen("v.txt", "w");	// open new file to write
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%.15f ", v[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("Master wrote [%d %d] elements.\n", n, n);

		// write info
		fp = fopen("pPnN.txt", "w");	// open new file to write
		fprintf(fp, "%d %d %d %d\n", p, P, n, N);
		fclose(fp);
		printf("Master wrote pPnN.\n");

	} else {

		int go;
		MPI_Recv(	&go,
		            1,
		            MPI_INT,
		            meshRank - 1,
		            666,	// tag
		            MPI_COMM_WORLD,
		            MPI_STATUS_IGNORE);

		size_t firstelem = 1;
		size_t lastelem = n;

		// write u
		fp = fopen("u.txt", "a");	// open to append
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%.15f ", u[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		// write v
		fp = fopen("v.txt", "a");	// open to append
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%.15f ", v[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("Process %d appended [%d %d] elements.\n", meshRank, n, n);
	}

	if (meshRank < P - 1) {
		MPI_Send(	&meshRank,	// send my meshRank to next process
		            1,
		            MPI_INT,
		            meshRank + 1,
		            666,
		            MPI_COMM_WORLD);
	}

	//=========================================================================
	// Cleanup
	//=========================================================================

	// free the ghost spirits
	/*for (int i = 0; i < NDIRS; ++i)
		MPI_Type_free(&ghost[i]);*/

	MPI_Finalize();

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
	int nsq = (int) sqrt(N);
	if (nsq * nsq != N)
	{
		if (!rank)
		{
			fprintf(stderr, "n: %d. n must be a square, exiting.\n", n );
		}
		return 1;
	}
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