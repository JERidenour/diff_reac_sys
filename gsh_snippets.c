/* Same as ImT but (I + T) */
/*void create_IpT(double *T, size_t n, size_t N, double D, double hx, double ht) {

	double ru = (D / (hx * hx)) * (ht * 0.5);

	double s = 1 + (ru * (-4)) ;
	double s2 = ru ;


	size_t skip = n - 1 ;

	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		// err += cs_entry(T, k, k, s) ;
		// T[k][k] = s ;
		T[N * k + k] = s ;
	}
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
		// if (k % n != skip) T[k][k+1] = s2 ;
		if (k % n != skip) T[(N * k) + k + 1] = s2 ;
	}
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		// if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
		// if (k % n != skip) T[k+1][k] = s2 ;
		if (k % n != skip) T[N * (k + 1) + k] = s2 ;
	}
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k + n, k, s2) ;
		// T[k + n][k] = s2 ;
		T[N * (k + n) + k] = s2 ;
	}
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		// err += cs_entry(T, k, k + n , s2) ;
		// T[k][k + n] = s2 ;
		T[(N * k) + k + n] = s2 ;
	}
}*/