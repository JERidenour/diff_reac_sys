#ifndef CN_H
#define CN_H

#include "CSparse/Include/cs.h"

/* Pass correctly allocated cs *T 
which will point to a new n^2xn^2 matrix when successful. 
Lord knows what abomination it points to if failed.
Returns 0 if all went well, or a positive integer signifying
the number of times that cs_entry() failed.*/
int create_ImT(cs *T, csi n, csi N, csi nz, double D, double hx, double ht) {

	double ru = (D / (hx * hx)) * (ht * 0.5);

	double s = 1 - (ru * (-4)) ;
	double s2 = -ru ;

	csi skip = n - 1 ;

	int err = 0;

	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		err += cs_entry(T, k, k, s) ;
	}
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
	}
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
	}
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		err += cs_entry(T, k + n, k, s2) ;
	}
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		err += cs_entry(T, k, k + n , s2) ;
	}

	return err ;
}

/* Same as ImT but (I + T) */
int create_IpT(cs *T, csi n, csi N, csi nz, double D, double hx, double ht) {

	double ru = (D / (hx * hx)) * (ht * 0.5);

	double s = 1 + (ru * (-4)) ;
	double s2 = ru ;

	csi skip = n - 1 ;

	int err = 0;

	// add main diagonal, i = k, j=k
	for (int k = 0; k < N; ++k)
	{
		err += cs_entry(T, k, k, s) ;
	}
	// add sup, i = k, j = k+1
	for (int k = 0; k < N - 1; ++k)
	{
		if (k % n != skip) err += cs_entry(T, k, k + 1, s2) ;
	}
	// add sub, i = k+1, j = k
	for (int k = 0; k < N - 1; ++k)
	{
		if (k % n != skip) err += cs_entry(T, k + 1, k, s2) ;
	}
	// add subsub, i = k+n+1, j = k
	for (int k = 0; k < N - n; ++k)
	{
		err += cs_entry(T, k + n, k, s2) ;
	}
	// add supsup, i = k, j = k+n+1
	for (int k = 0; k < N - n; ++k)
	{
		err += cs_entry(T, k, k + n , s2) ;
	}

	return err ;
}

#endif