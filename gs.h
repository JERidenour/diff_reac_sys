#ifndef GS_H
#define GS_H

void print_matrix(double *T, N) {
	printf("[\n");
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			printf("%.4f ", T[N * i + j] );
		}
		printf("\n");
	}
	printf("]\n");
}

void print_vector(double *u, N) {
	printf("[\n");
	for (int i = 0; i < N; ++i)
	{
		printf("%.4f ", u[i] );
	}
	printf("\n]\n");
}


/*	N is length of matrix.
	n is length of block.
	*T is pointer to zero-allocated array.
	This is a constant matrix.
*/
void create_ImT(double *T, size_t n, size_t N, double D, double hx, double ht) {

	double ru = (D / (hx * hx)) * (ht * 0.5);

	double s = 1 - (ru * (-4)) ;
	double s2 = -ru ;

	size_t skip = n - 1 ;

	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		//err += cs_entry(T, k, k, s) ;
		// T[k][k] = s ;
		// T + (N*k) + k = s ;
		T [(N * k) + k] = s ;
	}
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
		// if (k % n != skip) T[k][k+1] = s2 ;
		// if (k % n != skip) T + (N*k) + k + 1 = s2 ;
		if (k % n != skip) T[(N * k) + k + 1] = s2 ;

	}
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
		// if (k % n != skip) T[k+1][k] = s2 ;
		// if (k % n != skip) T + N*(k+1) + k = s2 ;
		if (k % n != skip) T [N * (k + 1) + k] = s2 ;
	}
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k + n, k, s2) ;
		// T[k+n][k] = s2 ;
		// T + N*(k+n) + k = s2;
		T[N * (k + n) + k] = s2;
	}
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k, k + n , s2) ;
		// T[k][k+n] = s2;
		// T + (N*k) + k + n = s2 ;
		T [(N * k) + k + n] = s2 ;
	}
}



/* Needs *rhs array of length N */
void create_RHSu(	double *rhs,
                    /*const double *T,*/
                    const double *u,
                    const double *v,
                    const double *uold,
                    const double *vold,
                    size_t n,
                    size_t N,
                    double D,
                    double F,
                    double hx,
                    double ht) {
	// IpTu*u
	// + 1.5*ht*(-u.*(v.^2) + F*(1-u)) - 0.5*ht*(-uold.*(vold.^2) + F*(1-uold))

	double ru = (D / (hx * hx)) * (ht * 0.5) ;
	double s = 1 + (ru * (-4)) ;
	double s2 = ru ;

	double alfa = 1.5 * ht ;
	double beta = 0.5 * ht ;

	size_t skip = n - 1 ;

	for (int j = 0; j < N; ++j)
	{
		rhs[j] = 0;
	}
	// printf("Main diagonal in RHSu.\n");
	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		// err += cs_entry(T, k, k, s) ;
		//T[k][k] = s * u[k] ;
		rhs[k] += s * u[k] ;
	}
	// printf("Sup diagonal in RHSu.\n");
	// printf("N: %d\n", N);
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
		// if (k % n != skip) T[k][k+1] = s2 * u[k+1] ;
		//printf("k: %d\n",k);

		if (k % n != skip) rhs[k] += s2 * u[k + 1] ;
	}
	// printf("Sub diagonal in RHSu.\n");
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
		// if (k % n != skip) T[k+1][k] = s2 * u[k] ;
		if (k % n != skip) rhs[k + 1] += s2 * u[k] ;
	}
	// printf("Subsub diagonal in RHSu.\n");
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k + n, k, s2) ;
		//T[k+n][k] = s2 * u[k];
		rhs[k + n] += s2 * u[k];
	}
	// printf("Supsup diagonal in RHSu.\n");
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k, k + n , s2) ;
		// T[k][k+n] = s2 * u[k+n];
		rhs[k] += s2 * u[k + n];
	}

	/*printf("IpTuu:\n");
	print_vector(rhs,N);*/

	for (int j = 0; j < N; ++j)
	{
		rhs[j] += alfa * ( -u[j] * v[j] * v[j] + F * (1 - u[j]) ) ;
	}

	for (int j = 0; j < N; ++j)
	{
		rhs[j] -= beta * ( -uold[j] * vold[j] * vold[j] + F * (1 - uold[j]) ) ;
	}
}

void create_RHSv(	double *rhs,
                    const double *u,
                    const double *v,
                    const double *uold,
                    const double *vold,
                    size_t n,
                    size_t N,
                    double D,
                    double F,
                    double K,
                    double hx,
                    double ht) {
	// RHSv = sparse(IpTv*v + 1.5*ht*(u.*(v.^2) - (F+k)*v) ...
	// - 0.5*ht*(uold.*(vold.^2) - (F+k)*vold));

	double ru = (D / (hx * hx)) * (ht * 0.5) ;
	double s = 1 + (ru * (-4)) ;
	double s2 = ru ;

	double alfa = 1.5 * ht ;
	double beta = 0.5 * ht ;

	size_t skip = n - 1 ;

	// set the AB term
	// for (int j = 0; j < N; ++j)
	// {
	// 	rhs[j] = alfa * ( u[j] * v[j] * v[j] - (F + K) * v[j] ) - beta * ( uold[j] * vold[j] * vold[j] - (F + K) * vold[j] ) ;
	// }
	// print_vector(rhs,N);
	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		// err += cs_entry(T, k, k, s) ;
		//T[k][k] = s * u[k] ;
		rhs[k] += s * v[k] ;
	}
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
		// if (k % n != skip) T[k][k+1] = s2 * u[k+1] ;
		if (k % n != skip) rhs[k] += s2 * v[k + 1] ;
	}
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
		// if (k % n != skip) T[k+1][k] = s2 * u[k] ;
		if (k % n != skip) rhs[k + 1] += s2 * v[k] ;
	}
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k + n, k, s2) ;
		//T[k+n][k] = s2 * u[k];
		rhs[k + n] += s2 * v[k];
	}
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k, k + n , s2) ;
		// T[k][k+n] = s2 * u[k+n];
		rhs[k] += s2 * v[k + n];
	}

	for (int j = 0; j < N; ++j)
	{
		rhs[j] += alfa * ( u[j] * v[j] * v[j] - (F + K) * v[j] ) ;
	}

	for (int j = 0; j < N; ++j)
	{
		rhs[j] -= beta * ( uold[j] * vold[j] * vold[j] - (F + K) * vold[j] ) ;
	}

	/*printf("IpTvv:\n");
	print_vector(rhs,N);*/
}

void gauss_seidel(	const double *A,
                    const double *b,
                    double *x,
                    const double *x0,
                    size_t N,
                    size_t maxiter) {

	double omega = 1;	// 1 => GS

	for (int i = 0; i < N; ++i)
	{
		x[i] = x0[i];
	}

	for (int iter = 0; iter < maxiter; ++iter)
	{
		double sigma;
		double left, right;
		for (int i = 0; i < N; ++i)
		{
			sigma = 0;
			for (int j = 0; j < N; ++j)
			{
				 if (j != i) sigma += A[N * i + j] * x[j];
				//sigma += A[N * i + j] * x[j];
			}
			// x[i] = (1.0 / A[N * i + i]) * (b[i] - sigma);	// pure GS
			
			x[i] = (1-omega)*x[i] + (omega/A[N * i + i])*(b[i]-sigma);
		}
	}
}



#endif