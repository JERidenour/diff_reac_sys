#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ndirs 8	// how many neighbours
#define n_min 1	// minimum meaningful size of local subdomain

/*#define NBR_NORTH hood[0][1]
#define NBR_EAST hood[1][2]
#define NBR_SOUTH hood[2][1]
#define NBR_WEST hood[1][0]
#define NBR_NORTHEAST hood[0][2]
#define NBR_SOUTHEAST hood[2][2]
#define NBR_SOUTHWEST hood[2][0]
#define NBR_NORTHWEST hood[0][0]
#define NBR_CENTER hood[1][1]*/

/*#define BORDER_EAST_START {1,n}
#define BORDER_EAST_SIZE {n,1}
#define GHOST_EAST_START {1,n + 1}
#define GHOST_EAST_SIZE {n,1}*/

#define DATA_BORDER_NORTH localData[1][1]			// 	north starts
#define DATA_BORDER_EAST localData[1][n]			//	east starts
#define DATA_BORDER_SOUTH localData[n][1]			// 	south starts
#define DATA_BORDER_WEST localData[1][1]			//	west starts
#define DATA_BORDER_NORTHEAST localData[1][n]			// 	NE starts
#define DATA_BORDER_SOUTHEAST localData[n][n]		 	// 	SE starts
#define DATA_BORDER_SOUTHWEST localData[n][1]			// 	SW starts
#define DATA_BORDER_NORTHWEST localData[1][1]			// 	NW starts

#define DATA_GHOST_NORTH localData[0][1]
#define DATA_GHOST_EAST localData[1][n + 1]		//	east starts
#define DATA_GHOST_SOUTH localData[n + 1][1]		// 	south starts
#define DATA_GHOST_WEST localData[1][0]			//	west starts
#define DATA_GHOST_NORTHEAST localData[0][n + 1]		// 	NE starts
#define DATA_GHOST_SOUTHEAST localData[n + 1][n + 1] 	// 	SE starts
#define DATA_GHOST_SOUTHWEST localData[n + 1][0]		// 	SW starts
#define DATA_GHOST_NORTHWEST localData[0][0]			// 	NW starts

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
	if (argc > 1 && atoi(argv[1]) > P * n_min )
		N = atoi(argv[1]);	// use argument if big enough
	else
	{
		N = P * n_min;
		if (!rank)
		{
			fprintf(stdout, "Setting N to minimum: %d x %d.\n", P, n_min );
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
	int neighbour[ndirs];	// neighbour[DIRECTION] = rank of neighbour
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
	MPI_Datatype ghost[ndirs];

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
	for (int i = 0; i < ndirs; ++i)
	{
		int totalDim[2]  = {n + 2, n + 2};	// dim of local data array
		int starts[2] = {0, 0};
		int subsizes[2]  = {ghostSizes[(i * 2)], ghostSizes[(i * 2) + 1]};
		MPI_Type_create_subarray(2, totalDim, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ghost[i]);
		MPI_Type_commit(&ghost[i]);
	}

	// array of the border(outgoing) submatrices numbered clockwise N E S W NE SE SW NW
	MPI_Datatype border[ndirs];

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
	for (int i = 0; i < ndirs; ++i)
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

	// TODO: computation

	// make test matrix
	double **localData = alloc2DArray((n + 2) * (n + 2));
	for (int i = 0; i < n + 2; ++i)
	{
		for (int j = 0; j < n + 2; ++j)
		{
			localData[i][j] = (i * (n + 2)) + j + 1 + rank * 100;
		}
	}

	// send stuff to neighbours
	int err = 0;
	/*int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	            int dest, int sendtag,
	            void *recvbuf, int recvcount, MPI_Datatype recvtype,
	            int source, int recvtag,
	            MPI_Comm comm, MPI_Status *status)*/

	MPI_Status status;
	int rcount;
	//int sizes[2]  = {n + 2, n + 2};	// dim of local data array
	//MPI_Datatype border[NORTH];
	//MPI_Datatype ghost[SOUTH];

	if (red)
	{
		// talk to north
		err = MPI_Sendrecv(	&(DATA_BORDER_NORTH),
		                    1,
		                    border[NORTH],
		                    neighbour[NORTH],
		                    VERTICAL,
		                    &(DATA_GHOST_NORTH),
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
		err = MPI_Sendrecv(	&(DATA_BORDER_SOUTH),
		                    1,
		                    border[SOUTH],
		                    neighbour[SOUTH],
		                    VERTICAL,
		                    &(DATA_GHOST_SOUTH),
		                    1,
		                    ghost[SOUTH],
		                    neighbour[SOUTH],
		                    VERTICAL,
		                    processMesh,
		                    &status
		                  );

		// talk to east
		err = MPI_Sendrecv(	&(DATA_BORDER_EAST),
		                    n,
		                    border[EAST],
		                    neighbour[EAST],
		                    HORIZONTAL,
		                    &(DATA_GHOST_EAST),
		                    n,
		                    ghost[EAST],
		                    neighbour[EAST],
		                    HORIZONTAL,
		                    processMesh,
		                    &status
		                  );

		// talk to west
		err = MPI_Sendrecv(	&(DATA_BORDER_WEST),
		                    n,
		                    border[WEST],
		                    neighbour[WEST],
		                    HORIZONTAL,
		                    &(DATA_GHOST_WEST),
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
		err = MPI_Sendrecv(	&(DATA_BORDER_SOUTH),
		                    1,
		                    border[SOUTH],
		                    neighbour[SOUTH],
		                    VERTICAL,
		                    &(DATA_GHOST_SOUTH),
		                    1,
		                    ghost[SOUTH],
		                    neighbour[SOUTH],
		                    VERTICAL,
		                    processMesh,
		                    &status
		                  );

		// talk to north
		err = MPI_Sendrecv(	&(DATA_BORDER_NORTH),
		                    1,
		                    border[NORTH],
		                    neighbour[NORTH],
		                    VERTICAL,
		                    &(DATA_GHOST_NORTH),
		                    1,
		                    ghost[NORTH],
		                    neighbour[NORTH],
		                    VERTICAL,
		                    processMesh,
		                    &status
		                  );

		// talk to west
		err = MPI_Sendrecv(	&(DATA_BORDER_WEST),
		                    n,
		                    border[WEST],
		                    neighbour[WEST],
		                    HORIZONTAL,
		                    &(DATA_GHOST_WEST),
		                    n,
		                    ghost[WEST],
		                    neighbour[WEST],
		                    HORIZONTAL,
		                    processMesh,
		                    &status
		                  );

		// talk to east
		err = MPI_Sendrecv(	&(DATA_BORDER_EAST),
		                    n,
		                    border[EAST],
		                    neighbour[EAST],
		                    HORIZONTAL,
		                    &(DATA_GHOST_EAST),
		                    n,
		                    ghost[EAST],
		                    neighbour[EAST],
		                    HORIZONTAL,
		                    processMesh,
		                    &status
		                  );
		/*MPI_Get_count(&status, ghost[SOUTH], &rcount);
		printf("%d received: %d from SOUTH.\n", meshRank, rcount);*/

	}


	//=========================================================================
	// Output to terminal
	//=========================================================================

	if (meshRank == 0) {

		printf("MeshRank %d reporting.\n", meshRank );
		printarr(localData, n + 2, "LocalData after comm.");

	} else {

		int go;
		MPI_Recv(	&go,
		            1,
		            MPI_INT,
		            rank - 1,
		            666,	// tag
		            MPI_COMM_WORLD,
		            MPI_STATUS_IGNORE);

		printf("MeshRank %d reporting.\n", meshRank );
		printarr(localData, n + 2, "LocalData after comm.");
	}

	if (rank < P - 1) {
		MPI_Send(	&rank,	// send my rank to next process
		            1,
		            MPI_INT,
		            rank + 1,
		            666,
		            MPI_COMM_WORLD);
	}

	//=========================================================================
	// Output to file
	//=========================================================================

	/*
	// int globalsizes[2] = {p*(n + 2), p*(n + 2)};	// dims of global array including ghosts
	int globalsizes[2] = {p*n, p*n};	// dims of global array including ghosts
	int localsizes[2] = {n, n};				// dims of local array ex ghosts
	int globalstarts[2] = {meshCoords[I]*n, meshCoords[J]*n};

	MPI_Datatype savedomain;
	err = MPI_Type_create_subarray(2, globalsizes, localsizes,
	                                 globalstarts, MPI_ORDER_C, MPI_DOUBLE, &savedomain);
	MPI_Type_commit(&savedomain);

	MPI_File mpifile;
	err = MPI_File_open(processMesh, "out.dat",
	                      MPI_MODE_WRONLY | MPI_MODE_CREATE,
	                      MPI_INFO_NULL, &mpifile);
	err = MPI_File_set_view(mpifile, 0, MPI_DOUBLE,
	                          savedomain, "native", MPI_INFO_NULL);
	err = MPI_File_write_all(mpifile, &(DATA_BORDER_NORTH),
	                           6 * 8, MPI_DOUBLE, &status);
	err = MPI_File_close(&mpifile);
	*/


	FILE *fp;
	//MPI_Status rstat;

	if (meshRank == 0) {
		fp = fopen("out.txt", "w");	// open new file to write
		size_t firstelem = 1;
		size_t lastelem = n;
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%f ", localData[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("Master wrote [%d %d] elements.\n", n,n);
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

		fp = fopen("out.txt", "a");	// open to append
		size_t firstelem = 1;
		size_t lastelem = n;
		for (int i = firstelem; i <= lastelem; ++i) {
			for (int j = 1; j <= lastelem; ++j)
			{
				fprintf(fp, "%f ", localData[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("Process %d appended [%d %d] elements.\n", meshRank, n,n);
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
	/*for (int i = 0; i < ndirs; ++i)
		MPI_Type_free(&ghost[i]);*/

	MPI_Finalize();

	// free local array
	/*free(localData[0]);
	free(localData);*/

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