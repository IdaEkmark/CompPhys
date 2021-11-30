#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#define PI 3.141592653589

void getNewConfig(double *R_m, double *R_t, double d, gsl_rng * r) {
	int coordinate = floor(gsl_rng_uniform(r)*6);
	for (int j=0; j<6; j++){
		if (j == coordinate) {
			R_t[j] = R_m[j] + 2*d*gsl_rng_uniform(r) - d;
		} else {
			R_t[j] = R_m[j];
		}
	}
	
}

double evalWaveFun(double *R, double alpha) {
    double r1 = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
    double r2 = sqrt(R[3]*R[3] + R[4]*R[4] + R[5]*R[5]);
    double r12 = sqrt( (R[0]-R[3])*(R[0]-R[3]) + (R[1]-R[4])*(R[1]-R[4]) + (R[2]-R[5])*(R[2]-R[5]) );

    return exp(-2.0*(r1+r2) + 1.0 / (2.0*( 1.0/r12 + alpha )) );
}

double evalLocalEnergy(double *R, double alpha) {
    double r1 = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
    double r2 = sqrt(R[3]*R[3] + R[4]*R[4] + R[5]*R[5]);
    double r12 = sqrt( (R[0]-R[3])*(R[0]-R[3]) + (R[1]-R[4])*(R[1]-R[4]) + (R[2]-R[5])*(R[2]-R[5]) );
    double denom = 1 + alpha * r12;
    double dotProd = (R[0]/r1 - R[3]/r2) * (R[0]-R[3]) + (R[1]/r1 - R[4]/r2) * (R[1]-R[4]) + (R[2]/r1 - R[5]/r2) * (R[2]-R[5]);

    return -4.0 + dotProd / (r12 * denom*denom) - 1.0/(r12 * denom*denom*denom) - 1.0/(4.0 * denom*denom*denom*denom) + 1/r12;
}

double evalDensityNormalization(double alpha, int N_norm, double d, gsl_rng * r) {
	double I = 0; double g_i; double psi_i; double R_i[6]; double xi_i; 
	for (int i=0; i<N_norm; i++) {
		for (int j=0; j<6; j++) {
			xi_i = 2*gsl_rng_uniform(r);
			R_i[j] = 2*10*d/PI * acos(1 - xi_i) - 10*d;
		}
		psi_i = evalWaveFun(R_i, alpha);
		g_i = psi_i * psi_i;
		for (int j=0; j<6; j++) {
			g_i *= 2/PI / sin(PI * R_i[j]); 
		}
		I += g_i;
	}
	I /= N_norm;
	return I;
}

double evalDensityFun(double *R, double alpha, double densityNorm) {
	double psi = evalWaveFun(R, alpha);
	double rho = psi * psi / densityNorm;
	return rho;
}

double metropolis(double d, gsl_rng * r, double alpha, int N, int N_norm) {
	double R_m[6] = {0, 0, 0, 0, 0, 0}; double R_t[6];
	double densityNorm = evalDensityNormalization(alpha, N_norm, d, r);
	
	double p_m = evalDensityFun(R_m, alpha, densityNorm); double p_t;
	/*
	double I = 0; double normfactor = 1/(N * fraction_use);
	int counter = 0;
	for (int i = 0; i < N; i++) {
		getNewConfig(R_m, R_t, d, r);
		
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
	gsl_rng_set(r, 12); double d = 1.0; double alpha = 0.1; int N = 10000; int N_norm = 10000;
	double k = metropolis(d, r, alpha, N, N_norm);
	
	return 0;
}