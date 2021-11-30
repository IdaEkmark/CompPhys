#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

void getNewConfig(double *R_m, double *R_t, double d, gsl_rng * r) {
	int coordinate = floor(gsl_rng_uniform(r)*6);
	for (int i=0; i<6; i++){
		if (i == coordinate) {
			R_t[i] = R_m[i] + 2*d*gsl_rng_uniform(r) - d;
		} else {
			R_t[i] = R_m[i];
		}
	}
	
}

double metropolis(double d, gsl_rng * r, double alpha, int N){
	double R_m[6] = {0, 0, 0, 0, 0, 0}; double R_t[6];
	getNewConfig(R_m, R_t, d, r);
	/*
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
	printf("accaptancerate = %.4f\n", (double)counter / N);*/
	return 1.0;//I;
}

int main(){
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r, 12); double d = 1.0; double alpha = 0.1; int N = 10000;
	double k = metropolis(d, r, alpha, N);
	
	return 0;
}