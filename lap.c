#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "lap.h"

// numerical parameters
#define Nglobal 81
/*#define hx 1.0/(Nglobal-1)
#define hx2 hx2*hx2*/
#define ht 0.0001
#define F 0.034
#define K 0.063
#define Du 2.0*1e-5
#define Dv 1*1e-5
#define u0 1.0
#define v0 0.0
#define MAXITER 2500


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
	//csi N = n * n;
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
	if (rank == 4) {
		// int r = N / 2;
		int m = n_inner;
		// printf("rank = %d \n", rank);
		// printf("r = %d, m = %d \n", r, m );
		for (int i = 1; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				u[1 + j * m] = 0.5;
				v[1 + j * m] = 1.25;
				u[i + 0 * m] = 0.5;
				v[i + 0 * m] = 1.25;
			}
		}
	}

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
	//cs *ImTu, *ImTv, *IpTu, *IpTv ;
	cs *IpTu, *IpTv ;

	//==========================
	// Initialize matrices
	//==========================

	csi nz = N + 2 * (N - 1) + 2 * (N - n) ;	// main + 2*(sub) + 2*(subsub) diagonals
	double hx = 1.0 / (Nglobal - 1);
	double hx2 = hx * hx;
	double su = Du * ht / (hx2);
	double sv = Dv * ht / (hx2);

	    // create ImT

	    int err = 0;	// should be 0 for no error

	// T = cs_spalloc (N, N, nz, 1, 1) ;
	// err = create_ImT(T, n, N, nz, Du, hx, ht);
	// ImTu = cs_compress(T) ;
	// cs_spfree (T) ;

	// T = cs_spalloc (N, N, nz, 1, 1) ;
	// err = create_ImT(T, n, N, nz, Dv, hx, ht);
	// ImTv = cs_compress(T) ;
	// cs_spfree (T) ;

	// create IpT

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Du, hx, ht);
	IpTu = cs_compress(T) ;
	cs_spfree (T) ;

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_IpT(T, n, N, nz, Dv, hx, ht);
	IpTv = cs_compress(T) ;
	cs_spfree (T) ;



	//==========================
	// Begin iteration
	//==========================

	for (int i = 0; i < MAXITER; i++) {

		//==========================
		// Integration step
		//==========================

		// RHSu = sparse(IpTu*u + ht*(-u.*(v.^2) + F*(1-u)));
		for (int i = 0; i < N; ++i) {
			unew[i] = ht * ( (-u[i] * v[i] * v[i]) + (F * (1.0 - u[i])) ) ;	// RHSu part 1
		}
		//add boundary values to rhs
		for (int i = 0; i < n_inner; i++) {
			unew[i + (n_inner - 1)*n_inner] += bin_north[i];
			unew[i] +=  bin_south[i];
			unew[(i + 1)*n_inner - 1] += bin_east[i];
			unew[i * n_inner] += bin_west[i];
		}
		err = cs_gaxpy(IpTu, u, unew) ; // RHSu part 2

		// RHSv = sparse(IpTv*v + ht*(u.*(v.^2) - (F+k)*v));
		for (int i = 0; i < N; ++i) {
			vnew[i] = ht * ( (u[i] * v[i] * v[i]) - ((F + K) * v[i]) ) ;	// RHSv part 1
		}
		//add boundary values to rhs
		for (int i = 0; i < n_inner; i++) {
			vnew[i + (n_inner - 1)*n_inner] += bin_north[i + n_inner];
			vnew[i] += bin_south[i + n_inner];
			vnew[(i + 1)*n_inner - 1] += bin_east[i + n_inner];
			vnew[i * n_inner] += bin_west[i + n_inner];
		}
		err = cs_gaxpy(IpTv, v, vnew) ;	// RHSv part 2

		//solve linear system (backslash)
		// err = cs_lusol(0, ImTu, unew, 0) ;
		// err = cs_lusol(0, ImTv, vnew, 0) ;

		for (int i = 0; i < N; ++i) {
			u[i] = unew[i];
			v[i] = vnew[i];
		}

		//==========================
		// communication step
		//==========================

		for (int i = 0; i < n_inner; i++) {
			bout_north[i] = -su * u[i + (n_inner - 1) * n_inner];
			bout_north[i + n_inner] = -sv * v[i + (n_inner - 1) * n_inner];

			bout_south[i] = -su * u[i];
			bout_south[i + n_inner] = -sv * v[i];

			bout_east[i] = -su * u[(i + 1) * n_inner - 1];
			bout_east[i + n_inner] = -sv * v[(i + 1) * n_inner - 1];

			bout_west[i] = -su * u[i * n_inner];
			bout_west[i + n_inner] = -sv * v[i * n_inner];
		}

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
			MPI_Send(
			    bout_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
			MPI_Send(
			    bout_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
			MPI_Send(
			    bout_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
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
			MPI_Recv(
			    bin_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
			MPI_Recv(
			    bin_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
			MPI_Recv(
			    bin_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
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
			MPI_Send(
			    bout_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD);
			MPI_Send(
			    bout_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD);
			MPI_Send(
			    bout_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
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
			MPI_Recv(
			    bin_south,
			    n_inner,
			    MPI_DOUBLE,
			    target_south,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
			MPI_Recv(
			    bin_east,
			    n_inner,
			    MPI_DOUBLE,
			    target_east,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
			MPI_Recv(
			    bin_west,
			    n_inner,
			    MPI_DOUBLE,
			    target_west,
			    0,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);
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
		fp = fopen("u_vals.txt", "w");	// open new file to write
		for (int i = 0; i < N_inner; ++i) {
			fprintf(fp, "%f\n", u[i]);
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

		fp = fopen("u_vals.txt", "a");	// open to append

		for (int i = 0; i < N_inner; ++i) {
			fprintf(fp, "%f\n", u[i]);
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

	//==========================
	// free resources
	//==========================

	// cs_spfree (ImTu) ;
	// cs_spfree (ImTv) ;
	cs_spfree (IpTu) ;
	cs_spfree (IpTv) ;

	free(u) ;
	free(v) ;
	free(unew) ;
	free(vnew) ;

	MPI_Finalize();
	return 0;
}