#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "cn.h"

// numerical parameters
#define Nglobal 64
/*#define hx 1.0/(Nglobal-1)
#define hx2 hx2*hx2*/
#define ht 1.0
#define F 0.034
#define K 0.065
#define Du 1.0
#define Dv 0.5
#define u0 1.0 ;
#define v0 0.0 ;
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
	if (Nglobal % Psq != 0) {
		fprintf(stdout, "Psq must divide Nglobal...\n");
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

	// length of solution vector, 
	// as well as one side of NxN laplacian matrix
	csi N = Nglobal/Psq;

	//printf("process rank: %d, tn: %d, ts: %d, te: %d, tw: %d \n", rank, target_north, target_south, target_east, target_west);
	double *u, *v, *unew, *vnew ;
	// allocate and initialize u and v
	u = (double*) malloc( N * sizeof(double)) ;
	v = (double*) malloc( N * sizeof(double)) ;
	for (int i = 0; i < N; ++i) {
		u[i] = u0;
		v[i] = v0;
	}

	// TODO: create inital values

	// allocate work vectors unew and vnew
	unew = (double*) malloc( N * sizeof(double)) ;
	vnew = (double*) malloc( N * sizeof(double)) ;

	//the long solution vector
	

	/*double u[N];
	double v[N];*/
	//example values (index)
	/*for(int i = 0; i < N; i++){
		u[i] = i;
	}*/

	csi n = sqrt(N); //length of boundary data

	// TODO: define boundary vectors for v

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

	// CSparse data
	cs *T ;
	cs *ImTu, *ImTv, *IpTu, *IpTv ;

	//==========================
	// Initialize matrices
	//==========================

	csi nz = N + 2 * (N - 1) + 2 * (N - n) ;	// main + 2*(sub) + 2*(subsub) diagonals
	double hx = 1.0/(Nglobal-1);
	/*double hx2 = hx2*hx2;*/

	// create ImT

	int err = 0;	// should be 0 for no error

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_ImT(T, n, N, nz, Du, hx, ht);
	ImTu = cs_compress(T) ;
	cs_spfree (T) ;

	T = cs_spalloc (N, N, nz, 1, 1) ;
	err = create_ImT(T, n, N, nz, Dv, hx, ht);
	ImTv = cs_compress(T) ;
	cs_spfree (T) ;

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

	//==========================
	// Integration step
	//==========================

	// RHSu = sparse(IpTu*u + ht*(-u.*(v.^2) + F*(1-u)));
	for (int i = 0; i < N; ++i) {
		unew[i] = ht * ( (-u[i] * v[i] * v[i]) + (F * (1.0 - u[i])) ) ;	// RHSu part 1
	}
	err = cs_gaxpy(IpTu, u, unew) ; // RHSu part 2

	// RHSv = sparse(IpTv*v + ht*(u.*(v.^2) - (F+k)*v));
	for (int i = 0; i < N; ++i) {
		vnew[i] = ht * ( (u[i] * v[i] * v[i]) - ((F + K) * v[i]) ) ;	// RHSv part 1
	}
	err = cs_gaxpy(IpTv, v, vnew) ;	// RHSv part 2

	// u_new = ImTu\RHSu;
	// v_new = ImTv\RHSv;
	err = cs_lusol(0, ImTu, unew, 0) ;
	err = cs_lusol(0, ImTv, vnew, 0) ;


	//==========================
	// communication step
	//==========================

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

	//==========================
	// end iteration
	//==========================

	cs_spfree (ImTu) ;
	cs_spfree (ImTv) ;
	cs_spfree (IpTu) ;
	cs_spfree (IpTv) ;

	free(u) ;
	free(v) ;
	free(unew) ;
	free(vnew) ;

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