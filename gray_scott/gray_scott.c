#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "gray_scott.h"
#include <time.h>

// numerical parameters
#define Nglobal 256*256
#define ht 0.19
#define F 0.0375
#define K 0.0634
#define Du 2.0*1e-5
#define Dv 1*1e-5
#define u0 1.0
#define v0 0.0
#define MAXITER 50000

int main(int argc, char *argv[])
{

	/* Initialize cpu timing */
	clock_t begin, end;
	double time_spent;
	begin = clock();

	//Psq is the number of processes, must be a square
	int P, Psq, rank, rc;
	int target_north, target_south, target_east, target_west;

	// MPI
	rc = MPI_Init(&argc, &argv);
	rc = MPI_Comm_size(MPI_COMM_WORLD, &Psq);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (Nglobal % Psq != 0) {
		fprintf(stdout, "Psq must divide Nglobal...\n");
		exit(1);
	}

	//determine position in process mesh:
	P = sqrt(Psq);	//process mesh dimension
	int north = rank > (P - 1) * P - 1;
	int south = rank < P;
	int east = (rank + 1) % P == 0;
	int west = rank % P == 0;

	//calculate rank of neighbors
	target_north = (1 - north) * (rank + P) + north * (rank - (P - 1) * P);
	target_south = (1 - south) * (rank - P) + south * (rank + (P - 1) * P);
	target_east = (1 - east) * (rank + 1) + east * (rank - P + 1);
	target_west = (1 - west) * (rank - 1) + west * (rank + P - 1);

	//dimension sizes
	csi N_inner = Nglobal / Psq;
	csi n_inner = sqrt(N_inner);
	csi n = n_inner;
	csi N = N_inner;

	double *u, *v, *unew, *vnew ;
	// allocate and initialize u and v
	u = (double*) calloc( N, sizeof(double)) ;
	v = (double*) calloc( N, sizeof(double)) ;
	for (int i = 0; i < N; ++i) {
		u[i] = u0;
		v[i] = v0;
	}
	// allocate work vectors unew and vnew
	unew = (double*) calloc( N , sizeof(double)) ;
	vnew = (double*) calloc( N , sizeof(double)) ;

	// initial values
	if (rank == 0 || rank == 3) {
		int c = floor(n_inner / 2);
		int r = floor(n_inner / 20);
		for (int i = c - r ; i < c + r; i++) {
			for (int j = c - r ; j < c + r; j++) {
				u[i + j * n_inner] = 0.5;
				v[i + j * n_inner] = 0.25;
			}
		}
	}
	if (rank == 2) {
		int c = floor(n_inner / 5);
		int r = floor(n_inner / 20);
		for (int i = c - r ; i < c + r; i++) {
			for (int j = c - r ; j < c + r; j++) {
				u[i + j * n_inner] = 0.5;
				v[i + j * n_inner] = 0.25;
			}
		}
	}
	if (rank == 1) {
		int c = floor(2*n_inner / 3);
		int r = floor(n_inner / 20);
		for (int i = c - r ; i < c + r; i++) {
			for (int j = c - r ; j < c + r; j++) {
				u[i + j * n_inner] = 0.5;
				v[i + j * n_inner] = 0.25;
			}
		}
	}

	//vectors receiving boudary data
	double bin_north[2 * n_inner], bin_south[2 * n_inner], bin_east[2 * n_inner], bin_west[2 * n_inner];
	//vectors sending boundary data
	double bout_north[2 * n_inner], bout_south[2 * n_inner], bout_east[2 * n_inner], bout_west[2 * n_inner];

	// CSparse data
	cs *T ;
	cs *IpTu, *IpTv ;

	//==========================
	// Initialize matrices
	//==========================

	csi nz = N + 2 * (N - 1) + 2 * (N - n) ;	// main + 2*(sub) + 2*(subsub) diagonals
	double hx = 1.0 / (sqrt(Nglobal) - 1);
	double hx2 = hx * hx;
	double su = Du * ht / (hx2);
	double sv = Dv * ht / (hx2);

	int err = 0;	// should be 0 for no error

	// create IpT

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Du, hx, ht);
	IpTu = cs_compress(T) ;
	cs_spfree (T) ;

	// double ru = (Du / (hx * hx)) * ht;
	// printf("ru = %f \n", ru);
	// printf("Du = %f \n", Du);
	// printf("hx = %f \n", hx);
	// printf("ht = %f \n", ht);

	// double rv = (Dv / (hx * hx)) * ht;
	// printf("rv = %f \n", rv);
	// printf("Dv = %f \n", Dv);

	// cs_print(IpTu, 0);

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Dv, hx, ht);
	IpTv = cs_compress(T) ;
	cs_spfree (T) ;

	// cs_print(IpTv, 0);


	//==========================
	// Begin iteration
	//==========================

	for (int i = 0; i < MAXITER; i++) {

		//==========================
		// communication step
		//==========================

		// build boundary output data
		for (int i = 0; i < n_inner; i++) {
			bout_north[i] = u[i + (n_inner - 1) * n_inner];
			bout_north[i + n_inner] = v[i + (n_inner - 1) * n_inner];

			bout_south[i] = u[i];
			bout_south[i + n_inner] = v[i];

			bout_east[i] = u[(i + 1) * n_inner - 1];
			bout_east[i + n_inner] = v[(i + 1) * n_inner - 1];

			bout_west[i] = u[i * n_inner];
			bout_west[i + n_inner] = v[i * n_inner];
		}

		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bout_north: %f \n", rank, bout_north[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bout_south: %f \n", rank, bout_south[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bout_east: %f \n", rank, bout_east[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bout_west: %f \n", rank, bout_west[i]);
		// }

		// *****************************
		// SEND/RECV FROM 0

		// send/recv 0 --> 2 (northward)
		if (rank == 0) {
			MPI_Send(
			    bout_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 2) {
			MPI_Recv(
			    bin_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 0 --> 2 (southward)
		if (rank == 0) {
			MPI_Send(
			    bout_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 2) {
			MPI_Recv(
			    bin_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 0 --> 1 (eastward)
		if (rank == 0) {
			MPI_Send(
			    bout_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 1) {
			MPI_Recv(
			    bin_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 0 --> 1 (westward)
		if (rank == 0) {
			MPI_Send(
			    bout_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 1) {
			MPI_Recv(
			    bin_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// *****************************
		// SEND/RECV FROM 1

		// send/recv 1 --> 3 (northward)
		if (rank == 1) {
			MPI_Send(
			    bout_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 3) {
			MPI_Recv(
			    bin_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 1 --> 3 (southward)
		if (rank == 1) {
			MPI_Send(
			    bout_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 3) {
			MPI_Recv(
			    bin_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 1 --> 0 (eastward)
		if (rank == 1) {
			MPI_Send(
			    bout_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 0) {
			MPI_Recv(
			    bin_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 1 --> 0 (westward)
		if (rank == 1) {
			MPI_Send(
			    bout_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 0) {
			MPI_Recv(
			    bin_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// *****************************
		// SEND/RECV FROM 2

		// send/recv 2 --> 0 (northward)
		if (rank == 2) {
			MPI_Send(
			    bout_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 0) {
			MPI_Recv(
			    bin_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 2 --> 0 (southward)
		if (rank == 2) {
			MPI_Send(
			    bout_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 0) {
			MPI_Recv(
			    bin_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 2 --> 3 (eastward)
		if (rank == 2) {
			MPI_Send(
			    bout_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 3) {
			MPI_Recv(
			    bin_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 2 --> 3 (westward)
		if (rank == 2) {
			MPI_Send(
			    bout_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 3) {
			MPI_Recv(
			    bin_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// *****************************
		// SEND/RECV FROM 3

		// send/recv 3 --> 1 (northward)
		if (rank == 3) {
			MPI_Send(
			    bout_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 1) {
			MPI_Recv(
			    bin_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 3 --> 1 (southward)
		if (rank == 3) {
			MPI_Send(
			    bout_south,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 1) {
			MPI_Recv(
			    bin_north,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 3 --> 2 (eastward)
		if (rank == 3) {
			MPI_Send(
			    bout_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 2) {
			MPI_Recv(
			    bin_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// send/recv 3 --> 2 (westward)
		if (rank == 3) {
			MPI_Send(
			    bout_west,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
		}
		if (rank == 2) {
			MPI_Recv(
			    bin_east,
			    2 * n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		//==========================
		// Update step
		//==========================

		//  update step one
		for (int i = 0; i < N; ++i) {
			unew[i] = u[i] + ht * ( (-u[i] * v[i] * v[i]) + (F * (1.0 - u[i])) ) ;
		}

		// update step two
		err = cs_gaxpy(IpTu, u, unew) ;

		/// update step three
		for (int i = 0; i < n_inner; i++) {
			unew[i + (n_inner - 1)*n_inner] += su * bin_north[i];
			unew[i] +=  su * bin_south[i];
			unew[(i + 1)*n_inner - 1] += su * bin_east[i];
			unew[i * n_inner] += su * bin_west[i];
		}

		//  update step one
		for (int i = 0; i < N; ++i) {
			vnew[i] = v[i] + ht * ( (u[i] * v[i] * v[i]) - ((F + K) * v[i]) ) ;
		}

		// update step two
		err = cs_gaxpy(IpTv, v, vnew) ;

		// update step three
		for (int i = 0; i < n_inner; i++) {
			vnew[i + (n_inner - 1)*n_inner] += sv * bin_north[i + n_inner];
			vnew[i] += sv * bin_south[i + n_inner];
			vnew[(i + 1)*n_inner - 1] += sv * bin_east[i + n_inner];
			vnew[i * n_inner] += sv * bin_west[i + n_inner];
		}

		// reset values
		for (int i = 0; i < N; ++i) {
			u[i] = unew[i];
			v[i] = vnew[i];
		}

		for (int i = 0; i < 2 * n_inner; i++) {
			bout_north[i] = 0;
			bout_south[i] = 0;
			bout_east[i] = 0;
			bout_west[i] = 0;
			bin_north[i] = 0;
			bin_south[i] = 0;
			bin_east[i] = 0;
			bin_west[i] = 0;
		}

		//==========================
		// end iteration
		//==========================
	}
	//==========================
	// print results
	//==========================

	FILE *fp;
	if (rank == 0) {


		char txt[50], rankstr[50], pro[50];
		strcpy(txt, ".txt");
		strcpy(pro, "pro_");
		sprintf(rankstr, "%d", rank);
		strcat(pro, rankstr);
		strcat(pro, txt);
		//printf("rank: %d, dest: %s\n",rank, dest);

		fp = fopen(pro, "w");	// open new file to write

		//TODO: make these print as matrices
		for (int i = 0; i < n_inner; ++i) {
			for (int j = 0; j < n_inner; j++) {
				fprintf(fp, "%f     ", u[j + i * n_inner]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

	} else {
		int go;
		MPI_Recv(	&go,
		            1,
		            MPI_INT,
		            rank - 1,
		            666,	// tag
		            MPI_COMM_WORLD,
		            MPI_STATUS_IGNORE);

		char txt[50], rankstr[50], pro[50];
		strcpy(txt, ".txt");
		strcpy(pro, "pro_");
		sprintf(rankstr, "%d", rank);
		strcat(pro, rankstr);
		strcat(pro, txt);
		//printf("rank: %d, dest: %s\n",rank, dest);

		fp = fopen(pro, "w");	// open new file to write

		for (int i = 0; i < n_inner; ++i) {
			for (int j = 0; j < n_inner; j++) {
				fprintf(fp, "%f     ", u[j + i * n_inner]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	if (rank < Psq - 1) {
		MPI_Send(	&rank,	// send my rank to next process
		            1,
		            MPI_INT,
		            rank + 1,
		            666,
		            MPI_COMM_WORLD);
	}


	// FILE *fp;
	// if (rank == 0) {
	// 	fp = fopen("u_vals_center.txt", "w");	// open new file to write
	// 	for (int i = 0; i < N_inner; ++i) {
	// 		fprintf(fp, "%f\n", u[i]);
	// 	}
	// 	fclose(fp);
	// }

	/* end cpu timing  */
	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
	printf("Process %d, time spent: %f \n", rank, time_spent);

	//==========================
	// free resources
	//==========================

	cs_spfree (IpTu) ;
	cs_spfree (IpTv) ;

	free(u) ;
	free(v) ;
	free(unew) ;
	free(vnew) ;

	MPI_Finalize();
	return 0;
}