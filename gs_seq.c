#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

#include "gs.h"

/*#define N 4			// NxN dim of domain
#define dt 0.5
#define Du 2e-5
#define Dv 1e-5
#define F 0.034
#define K 0.063*/


int main(int argc, char const *argv[])
{
	size_t nsteps = 1;
	double N = 4;
	double dt = 0.5;
	double Du = 0.00002;
	double Dv = 0.00001;
	/*double Du = 1.0;
	double Dv = 0.5;*/
	double F = 0.034;
	double K = 0.063;

	double N2 = N * N;	// N2xN2 dim of matrix

	double h = 1 / (N - 1);	// hx
	double h2 = h * h;

	// u and v work vectors
	double *u = (double*) calloc(N2, sizeof(double));		// zero initialized
	for (int i = 0; i < N2; ++i) u[i] = 1;					// set initial u to 1
	double *v = (double*) calloc(N2, sizeof(double));		// zero initialized
	//for (int i = 0; i < N2; ++i) v[i] = 1;
	double *unew = (double*) calloc(N2, sizeof(double));	// zero initialized
	double *vnew = (double*) calloc(N2, sizeof(double));	// zero initialized
	double *uold = (double*) calloc(N2, sizeof(double));	// zero initialized
	// for (int i = 0; i < N2; ++i) uold[i] = 1;
	double *vold = (double*) calloc(N2, sizeof(double));	// zero initialized
	//for (int i = 0; i < N2; ++i) vold[i] = 1;

	double *guess = (double*) calloc(N2, sizeof(double));
	for (int i = 0; i < N2; ++i) guess[i] = 0.5;
	//double ImTu[N2][N2];
	//double ImTv[N2][N2];

	double *ImTu = (double*) calloc(N2 * N2, sizeof(double));	// the matrix (massive waste of space)
	double *ImTv = (double*) calloc(N2 * N2, sizeof(double));	// the matrix (massive waste of space)

	double *IpTu = (double*) calloc(N2 * N2, sizeof(double));	// the matrix (massive waste of space)
	double *IpTv = (double*) calloc(N2 * N2, sizeof(double));	// the matrix (massive waste of space)

	double *rhsu = (double*) calloc(N2, sizeof(double));
	double *rhsv = (double*) calloc(N2, sizeof(double));

	printf("Creating ImT.\n");
	create_ImT(ImTu, N, N2, Du, h, dt);
	create_ImT(ImTv, N, N2, Dv, h, dt);

	/*printf("Creating IpT.\n");
	create_IpT(IpTu, N, N2, Du, h, dt);
	create_IpT(IpTv, N, N2, Dv, h, dt);*/

	/*printf("Creating RHSu.\n");
	create_RHSu(rhsu, u, v, uold, vold, N, N2, Du, F, h, dt);
	printf("Creating RHSv.\n");
	create_RHSv(rhsv, u, v, uold, vold, N, N2, Dv, F, K, h, dt);*/

	printf("ImTu:\n");
	print_matrix(ImTu,N2);
	printf("ImTv:\n");
	print_matrix(ImTv,N2);

	/*printf("IpTu:\n");
	print_matrix(IpTu,N2);
	printf("IpTv:\n");
	print_matrix(IpTv,N2);*/

	

	printf("u:\n");
	print_vector(u,N2);
	printf("v:\n");
	print_vector(v,N2);

	printf("uold:\n");
	print_vector(uold,N2);
	printf("vold:\n");
	print_vector(vold,N2);

	for (int step = 0; step < nsteps; ++step)
	{
		/*for (int i = 0; i < N2; ++i)
		{
			rhsu[i] = 0;
			rhsv[i] = 0;
		}*/
		create_RHSu(rhsu, u, v, uold, vold, N, N2, Du, F, h, dt);
		create_RHSv(rhsv, u, v, uold, vold, N, N2, Dv, F, K, h, dt);

		printf("RHSu:\n");
		print_vector(rhsu,N2);
		printf("RHSv:\n");
		print_vector(rhsv,N2);

		gauss_seidel(ImTu, rhsu, unew, guess, N2, 50);
		gauss_seidel(ImTv, rhsv, vnew, guess, N2, 50);

		for (int i = 0; i < N2; ++i)
		{
			uold[i] = u[i];
			vold[i] = v[i];
			u[i] = unew[i];
			v[i] = vnew[i];
		}

	}

	printf("unew:\n");
	print_vector(u, N2);
	printf("vnew:\n");
	print_vector(v, N2);

	return 0;
}
