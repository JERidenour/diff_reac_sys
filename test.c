#include <stdio.h>

int main(int argc, char const *argv[])
{
	int n_inner = 4;
	double u[n_inner * n_inner], unew[n_inner * n_inner];
	double bin_north[n_inner], bin_south[n_inner], bin_east[n_inner], bin_west[n_inner];

	for (int i = 0; i < n_inner; ++i)
	{
		bin_north[i] = 1;
		bin_south[i] = 1;
		bin_east[i] = 1;
		bin_west[i] = 1;

	}

	for (int i = 0; i < n_inner * n_inner; ++i)
	{
		unew[i] = 0;
	}

	for (int i = 0; i < n_inner; i++) {
		unew[i + (n_inner - 1)*n_inner] -= bin_north[i];
		unew[i] -= bin_south[i];
		unew[(i + 1)*n_inner - 1] -= bin_east[i];
		unew[i * n_inner] -= bin_west[i];
	}

	printf("unew: \n");
	for (int i = 0; i < n_inner * n_inner; ++i)
	{
		printf("%.4f \n", unew[i]);
	}
	return 0;
}