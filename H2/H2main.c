#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

void getNewConfig(double *r, double d, gsl_rng * r) {
	int coodinate = floor(gsl_rng_uniform(r)*6);
	r[coordinate] += 2*d*gsl_rng_uniform(r) - d;
}

double metropolis(double delta, gsl_rng * r, int N, double fraction_use){
	double x_m = 0; double x_t = x_m + delta * (gsl_rng_uniform(r) - 0.5);
	double y_m = 0; double y_t = y_m + delta * (gsl_rng_uniform(r) - 0.5); 
	double z_m = 0; double z_t = z_m + delta * (gsl_rng_uniform(r) - 0.5);
	double p_m = evalWeightFunction3(x_m, y_m, z_m); double p_t;
	double I = 0; double normfactor = 1/(N * fraction_use);
	int counter = 0;
	for (int i = 0; i < N; i++) {
		x_t = x_m + delta * (gsl_rng_uniform(r) - 0.5);
		y_t = y_m + delta * (gsl_rng_uniform(r) - 0.5); 
		z_t = z_m + delta * (gsl_rng_uniform(r) - 0.5);
		
		p_t = evalWeightFunction3(x_t, y_t, z_t);
		if (gsl_rng_uniform(r) < p_t/p_m) {
			x_m = x_t; 
			y_m = y_t;
			z_m = z_t;
			p_m = p_t;
			counter++;

		}
		if ( i > N * (1 - fraction_use) ) {
			I += evalFunction3(x_m, y_m, z_m) * normfactor;
		}
	}
	printf("accaptancerate = %.4f\n", (double)counter / N);
	return I;
}