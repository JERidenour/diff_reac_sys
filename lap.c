#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "lap.h"

// numerical parameters
#define Nglobal 4*64
#define ht 0.025
#define F 0.026
#define K 0.057
#define Du 2.0*1e-5
#define Dv 1*1e-5
#define u0 1.0
#define v0 0.0
#define MAXITER 1000


int main(int argc, char *argv[])
{

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
	if (Psq % 2 != 0) {
		fprintf(stdout, "Number of processes must be even...\n");
		exit(1);
	}

	//determine position in process mesh:
	P = sqrt(Psq);	//process mesh dimension
	int north = rank > (P - 1) * P - 1;
	int south = rank < P;
	int west = rank % P == 0;
	int east = (rank + 1) % P == 0;

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

	// initial values
	if (rank == 0 || rank == 3 || rank == 2) {
		int c = floor(N_inner / 2);
		int r = floor(n_inner / 2);
		u[c+r] = 0.5;
		v[c+r] = 1;
	}

	// for (int i = 0; i < N_inner; ++i) {
	// 	u[i] = i;
	// }

	// allocate work vectors unew and vnew
	unew = (double*) calloc( N , sizeof(double)) ;
	vnew = (double*) calloc( N , sizeof(double)) ;

	//vectors receiving boudary data
	double bin_north[2 * n_inner], bin_south[2 * n_inner], bin_east[2 * n_inner], bin_west[2 * n_inner];
	//vectors sending boundary data
	double bout_north[2 * n_inner], bout_south[2 * n_inner], bout_east[2 * n_inner], bout_west[2 * n_inner];

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

	// CSparse data
	cs *T ;
	cs *IpTu, *IpTv ;

	//==========================
	// Initialize matrices
	//==========================

	csi nz = N + 2 * (N - 1) + 2 * (N - n) ;	// main + 2*(sub) + 2*(subsub) diagonals
	double hx = 1.0 / (Nglobal - 1);
	double hx2 = hx * hx;
	double su = Du * ht / (hx2);
	double sv = Dv * ht / (hx2);

	int err = 0;	// should be 0 for no error

	// create IpT

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Du, hx, ht);
	IpTu = cs_compress(T) ;
	cs_spfree (T) ;

	//cs_print(IpTu, 0);

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Dv, hx, ht);
	IpTv = cs_compress(T) ;
	cs_spfree (T) ;

	//cs_print(IpTv, 0);


	//==========================
	// Begin iteration
	//==========================

	for (int i = 0; i < MAXITER; i++) {

		//==========================
		// Integration step
		//==========================

		//  update step one
		for (int i = 0; i < N; ++i) {
			unew[i] = ht * ( (-u[i] * v[i] * v[i]) + (F * (1.0 - u[i])) ) ;
		}
		/// update step two
		for (int i = 0; i < n_inner; i++) {
			unew[i + (n_inner - 1)*n_inner] += su * bin_east[i];
			unew[i] +=  su * bin_west[i];
			unew[(i + 1)*n_inner - 1] += su * bin_north[i];
			unew[i * n_inner] += su * bin_south[i];
		}
		// update step three
		err = cs_gaxpy(IpTu, u, unew) ;
		//printf("iteration: %d, rank: %d, cs_gaxpy(IpTu, u, unew):  %d \n", i, rank, err);

		//  update step one
		for (int i = 0; i < N; ++i) {
			vnew[i] = ht * ( (u[i] * v[i] * v[i]) - ((F + K) * v[i]) ) ;
		}
		// update step two
		for (int i = 0; i < n_inner; i++) {
			vnew[i + (n_inner - 1)*n_inner] += sv * bin_east[i + n_inner];
			vnew[i] += sv * bin_west[i + n_inner];
			vnew[(i + 1)*n_inner - 1] += sv * bin_north[i + n_inner];
			vnew[i * n_inner] += sv * bin_south[i + n_inner];
		}
		// update step three
		err = cs_gaxpy(IpTv, v, vnew) ;
		//printf("iteration: %d, rank: %d, cs_gaxpy(IpTv, v, vnew):  %d \n", i, rank, err);

		// reset values
		for (int i = 0; i < N; ++i) {
			u[i] = unew[i];
			v[i] = vnew[i];
		}

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

		// are we red? (if not, we are black)
		bool red = (rank % 2) ? false : true;
		if (red) {
			MPI_Send(
			    bout_north,
			    n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_north,
			    n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		} else {
			MPI_Send(
			    bout_north,
			    n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_north,
			    n_inner,
			    MPI_DOUBLE,
			    target_north,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

			MPI_Send(
			    bout_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD);
			MPI_Recv(
			    bin_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
		}

		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bin_north: %f \n", rank, bin_north[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bin_south: %f \n", rank, bin_south[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bin_east: %f \n", rank, bin_east[i]);
		// }
		// for (int i = 0; i < 2 * n_inner; i++) {
		// 	printf("rank: %d, bin_west: %f \n", rank, bin_west[i]);
		// }

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
				fprintf(fp, "%f     ", u[j+i*n_inner]);
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
				fprintf(fp, "%f     ", u[j+i*n_inner]);
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