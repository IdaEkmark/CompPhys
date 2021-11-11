/******************************************************************************
 * E1code4
 * 
 * Compile me as:
 * clang -c fft.c -o fft.o -lgsl -lgslcblas
 * clang E2code1.c fft.o -o ./Executable_files/<execuutable name> -lgsl -lgslcblas
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>
#include "fft.h"  // interface to fft routine

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm 
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double **V, double **Q,
            double dt, double *m, double kappa)
{
    double q[n_particles];
    double v[n_particles];
    double a[n_particles];
    q[0] = Q[0][0];
    q[1] = Q[1][0];
    q[2] = Q[2][0];
    v[0] = V[0][0];
    v[1] = V[1][0];
    v[2] = V[2][0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
		
        /* Save the displacement of the three atoms */
        Q[0][i] = q[0];
        Q[1][i] = q[1];
        Q[2][i] = q[2];
        V[0][i] = v[0];
		V[1][i] = v[1];
		V[2][i] = v[2];
    }
}

/*
 * constructs time array
 * @array - array to be filled with time values
 * @start - start value
 * @len_t - number of times stamps in array
 * @dt - time step between two consecutive times
*/
void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
    }
}

/*
 * Saves t, E, K and P to a file (Usually .csv).
 * @fname - File name 
 * @time_array - array of time values
 * @Q - n_particles x n_timesteps array with displacement values
 * @V - n_particles x n_timesteps array with velocity values
 * @n_timesteps - number of timesteps
 * @n_particles - number of particles
*/
void saveEKPtoFile(char *fname, double *time_array,
		   double *E, double *K, double *P, int n_timesteps) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "time, E, K, P\n");

    for(int i = 0; i < n_timesteps; ++i){
	    fprintf(fp, "%f", time_array[i]);
		fprintf(fp, ", %f", E[i]);
		fprintf(fp, ", %f", K[i]);
		fprintf(fp, ", %f", P[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/*
 * Saves t, q_i and v_i to a file (Usually .csv).
 * @fname - File name 
 * @time_array - array of time values
 * @Q - n_particles x n_timesteps array with displacement values
 * @V - n_particles x n_timesteps array with velocity values
 * @n_timesteps - number of timesteps
 * @n_particles - number of particles
*/
void saveQVtoFile(char *fname, double *time_array,
		   double **Q, double **V, int n_timesteps, int n_particles) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "time");
    for (int j = 0; j < n_particles; ++j) {
        fprintf(fp, ", Q%i, V%i", j, j);
    }
    fprintf(fp, "\n");

    for(int i = 0; i < n_timesteps; ++i){
	    fprintf(fp, "%f", time_array[i]);

        for(int j = 0; j < n_particles; ++j){
            fprintf(fp, ", %f, %f", Q[j][i], V[j][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/*
 * writes to file
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
 * @n_points - number of points
*/
void write_to_file(char *fname, double *frequencies,
		   double *spectrum, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "signal, frequency\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", spectrum[i], frequencies[i]);
    }
    fclose(fp);
}

void calculateEnergy(int n_timesteps, int n_particles, double **V, double **Q,
	     double dt, double *m, double kappa, double *K, double *P, double *E) 
{
	for (int t = 0; t < n_timesteps; t++) { 
		K[t] = 0;
		P[t] = 0;

		// Handle boundaries
		P[t] += kappa * Q[0][t] * Q[0][t] / 2.0;
		P[t] += kappa * Q[n_particles-1][t] * Q[n_particles-1][t] / 2.0;

		// Handle everything else
		for (int p = 0; p < n_particles; p++) {
			K[t] += m[p] * V[p][t] * V[p][t] / 2.0;
			if (p < n_particles-1) {
				P[t] += kappa * ( Q[p+1][t] - Q[p][t] ) * ( Q[p+1][t] - Q[p][t] ) / 2.0;
			}
		}
		E[t] = K[t] + P[t];
	}
}

void extractColumnVector(double **M, double *V, int index, int Mrows) {
	for (int i = 0; i < Mrows; i++) {
		V[i] = M[index][i];
	}
}

int main()
{
    double kn_per_m = 6.24151e1; // 1e3 J/m^2 * (1e-10 m/Å)^2 / (1.602e-19 J/eV)
    double tMax = 0.25;
    int n_t = 2.5e4; double dt = tMax / (double) n_t; int n_p = 3; double kappa = kn_per_m;
    double m[n_p];
	for (int i=0; i<n_p; i++) {
		m[i] = 12.0 / 9649.0;
	}
    double **Q; double **V;
    
	Q = malloc(n_p * sizeof *Q);
	for (int i=0; i < n_p; i++){
		Q[i] = malloc(n_t * sizeof *Q[i]);
	}
	V = malloc(n_p * sizeof *V);
	for (int i=0; i < n_p; i++){
		V[i] = malloc(n_t * sizeof *V[i]);
	}
	
	Q[0][0] = 0.01;
	Q[1][0] = 0.005;
	Q[2][0] = -0.005;
	V[0][0] = 0;
	V[1][0] = 0;
	V[2][0] = 0;
	
	velocity_verlet(n_t, n_p, V, Q, dt, m, kappa);
	
	double *time_array;
	time_array = malloc(n_t * sizeof (double));
    arange(time_array, 0, n_t, dt);
    saveQVtoFile("5/E1u5_dt2.5e4.csv", time_array, Q, V, n_t, n_p);


	double *frequencies;
	frequencies = malloc(n_t * sizeof (double));
	for (int i = 0; i < n_t; i++) {
		frequencies[i] = i / (dt * n_t) - 1.0 / (2.0 * dt);
	}

	/*
	 * Do the fft
	 */
	double *fftd_data1; double *fftd_data2; double *fftd_data3;
	fftd_data1 = malloc(n_t * sizeof (double));
	fftd_data2 = malloc(n_t * sizeof (double));
	fftd_data3 = malloc(n_t * sizeof (double));

	
	powerspectrum(Q[0], fftd_data1, n_t);
	powerspectrum_shift(fftd_data1, n_t);
	powerspectrum(Q[1], fftd_data2, n_t);
	powerspectrum_shift(fftd_data2, n_t);
	powerspectrum(Q[2], fftd_data3, n_t);
	powerspectrum_shift(fftd_data3, n_t);

	/*
     * Dump fft and frequencies to file
     */
    write_to_file("5/E1u5_PS_q1_dt2.5e4.csv", fftd_data1, frequencies, n_t);
	write_to_file("5/E1u5_PS_q2_dt2.5e4.csv", fftd_data2, frequencies, n_t);
	write_to_file("5/E1u5_PS_q3_dt2.5e4.csv", fftd_data3, frequencies, n_t);

	free(frequencies);

	double *E; double *K; double *P;
	E = malloc(n_t * sizeof (double));
	K = malloc(n_t * sizeof (double));
	P = malloc(n_t * sizeof (double));

	calculateEnergy(n_t, n_p, V, Q, dt, m, kappa, K, P, E);
	saveEKPtoFile("5/E1u5_Energies_dt2.5e4.csv", time_array, E, K, P, n_t);
	
	free(E);
	free(K);
	free(P);

	free(time_array);

	free(fftd_data1);
	free(fftd_data2);
	free(fftd_data3);

	for (int i=0; i < n_p; i++){
		free(Q[i]);
	}
	free(Q);
	
	for (int i=0; i< n_p; i++){
		free(V[i]);
	}
	free(V);
	return 0;
}