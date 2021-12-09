/* E4 main c file */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"  // interface to fft routine

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define K_B 1.380649e-23
#define PI 3.141592653589

void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
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
void brownian_verlet(int n_timesteps, double eta, double temperature, double mass, double omega0, double (*pos), double (*vel),
		double dt, gsl_rng * r)
{
    double q;
    double v;
    double a;

    double c0 = exp(-eta*dt);
    double v_th = sqrt(K_B * temperature / mass);

    double gaussian_r;
    
    q = pos[0];
    v = vel[0];

    // a(t)
    a = -omega0*omega0 * q;
    //printf("n_t: %i\n", n_timesteps);

    for (int t = 1; t < n_timesteps + 1; t++) {
        //printf("Time: %i\n",t);
        //printf("Pos: %.2f\n",q);

        // v(t+dt/2)
        gaussian_r = gsl_ran_gaussian(r, 1.0);
        v = dt * 0.5 * a + sqrt(c0) * v + v_th * sqrt(1-c0) * gaussian_r;
        
        // q(t+dt)
        q += dt * v;
        
        // a(t+dt)
        a = -omega0*omega0 * q;
        gaussian_r = gsl_ran_gaussian(r, 1.0);
        
        // v(t+dt)
        v = (dt * 0.5 * a + v) * sqrt(c0) + v_th * sqrt(1-c0) * gaussian_r;
		
        // Save the displacement of the three atoms
        pos[t] = q;
        vel[t] = v;
    }
}

/*
 * writes to file
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
 * @n_points - number of points
*/
void saveDataToFile(char *fname, double *yvals, double *tvals, int n_points, int n_skip)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "xData, yData\n");
    for(int i = 0; i < n_points; i = i + n_skip){
	    fprintf(fp, "%e, %e\n", tvals[i], yvals[i]);
    }
    fclose(fp);
}

/*
 * reads time_array and signal data from file
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
*/
void read_data(char *fname, double *time_array, double *signal, int startind, int skiplen)
{
    FILE *fp = fopen(fname, "r");

    /* if file no found
     * error out and exit code 1
     */
    if(fp == NULL){
	perror("error:");
	exit(1);
    }

    /* skip header */
    fseek(fp, strlen("xData, yData\n"), SEEK_SET);
    char line[128] = {0};
    char *token;
    int i = 0; int j = 0;
    while(fgets(line, sizeof(line), fp) != NULL){
        token = strtok(line, ",");
        if ((i % skiplen == 0) && (i >= startind)) {
            time_array[j] = strtod(token, NULL);
            //printf("%.2e\n", time_array[j]);
        }
        token = strtok(NULL, ",");
        if ((i % skiplen == 0) && (i >= startind)) {
            signal[j] = strtod(token, NULL);
            j++;
        }
        i++;
        memset(line, 0, sizeof(line));
        token = NULL;
    }
    fclose(fp);
}

void evalMean(double **inputmatrix, double *outputvector, int N, int nt) {
	for(int t=0; t<nt+1; t++) {
		outputvector[t] = 0;
		for(int i=0; i<N; i++) {
			outputvector[t] += inputmatrix[i][t];
		}
		outputvector[t] /= N;
	}
}

void evalVariance(double **inputmatrix, double *outputvector, double *mean, int N, int nt) {
	for(int t=0; t<nt+1; t++) {
		outputvector[t] = 0;
		for(int i=0; i<N; i++) {
			outputvector[t] += pow((inputmatrix[i][t] - mean[t]),2);
		}
		outputvector[t] /= N;
	}
}

void evalFokkerPlanck(double **X, double **V, double **f, int t, double dx, double dv, int nx, int nv,
		double x_max, double x_min, double v_max, double v_min, int N) 
{
	double x; double v; int ix; int iv;
	
	for (int ix=0; ix<nx; ix++) {
		for (int iv=0; iv<nv; iv++) {
			f[ix][iv] = 0;
		}
	}

	for (int i=0; i<N; i++) {
		x = X[i][t];
		v = V[i][t];
		
		ix = floor((x-x_min)/(x_max-x_min)*nx);
		iv = floor((v-v_min)/(v_max-v_min)*nv);
		f[ix][iv] += 1/((double)N);
	}
}

void evalDensityFunctions(double **X,  double **V, double *rhoX, double *rhoV, int t, double dx, double dv, int nx, int nv,
		double x_max, double x_min, double v_max, double v_min, int N)
{
	int ix; int iv;
	double **f;
	f = malloc(nx * sizeof *f);
	for (int i=0; i < nx; i++){
		f[i] = malloc((nv) * sizeof *f);
	}
	
	evalFokkerPlanck(X, V, f, t, dx, dv, nx, nv, x_max, x_min, v_max, v_min, N);
	
	for (ix=0; ix<nx; ix++) {
		rhoX[ix] = 0;
		for (iv=0; iv<nv; iv++) {
			rhoX[ix] += f[ix][iv]*dv;
		}
	}
	for (iv=0; iv<nv; iv++) {
		rhoV[iv] = 0;
		for (ix=0; ix<nx; ix++) {
			rhoV[iv] += f[ix][iv]*dx;
		}
	}
	
	for (int i=0; i < nx; i++){
		free(f[i]);
	}
	free(f);
}

void runtask2trajectory() {
    double dt = 1e-6; double t_max = 4e-3; int nt = (int) (t_max/dt); 
    double tau_low = 147.3e-6; double tau_high = 48.5e-6; double eta_low = 1/tau_low; double eta_high = 1/tau_high; double omega0 = 3.1e3*2*PI;
    double r = 2.79e-6/2; double rho = 2.65e3; double m = (4.0/3.0 * PI * pow(r,3)) * rho; double T = 297.0; 
    double x0 = 0.1e-6; double v0 = 2e-3; 
    double *X_low; double *V_low; double *X_high; double *V_high; double *time;

    X_low = malloc((nt+1) * sizeof(double));
    V_low = malloc((nt+1) * sizeof(double));
    X_high = malloc((nt+1) * sizeof(double));
    V_high = malloc((nt+1) * sizeof(double));

    X_low[0] = x0;
    V_low[0] = v0;
    X_high[0] = x0;
    V_high[0] = v0;
    
    gsl_rng * rg = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rg, 12*5);

    brownian_verlet(nt, eta_low, T, m, omega0, X_low, V_low, dt, rg);
    brownian_verlet(nt, eta_high, T, m, omega0, X_high, V_high, dt, rg);

    time = malloc((nt+1) * sizeof(double));
	arange(time, 0, nt+1, dt);

    saveDataToFile("2/position_tau147.3e-6_5.csv", X_low, time, nt+1, 1);
    saveDataToFile("2/velocity_tau147.3e-6_5.csv", V_low, time, nt+1, 1);
    saveDataToFile("2/position_tau48.5e-6_5.csv", X_high, time, nt+1, 1);
    saveDataToFile("2/velocity_tau48.5e-6_5.csv", V_high, time, nt+1, 1);

    free(time);
    free(X_low);
    free(V_low);
    free(X_high);
    free(V_high);
}

void runtask2mean() {
    double dt = 1e-6; double t_max = 4e-3; int nt = (int) (t_max/dt); 
    double tau_low = 147.3e-6; double tau_high = 48.5e-6; double eta_low = 1/tau_low; double eta_high = 1/tau_high; double omega0 = 3.1e3*2*PI;
    double r = 2.79e-6/2; double rho = 2.65e3; double m = (4.0/3.0 * PI * pow(r,3)) * rho; double T = 297.0; 
    double x0 = 0.1e-6; double v0 = 2e-3; 
    double **X_low; double **V_low; double **X_high; double **V_high; double *time; int N = 10000; int i;
    
    X_low = malloc(N * sizeof *X_low);
	for (i=0; i < N; i++){
		X_low[i] = malloc((nt+1) * sizeof *X_low[i]);
	}
	V_low = malloc(N * sizeof *V_low);
	for (int i=0; i < N; i++){
		V_low[i] = malloc((nt+1) * sizeof *V_low[i]);
	}
    X_high = malloc(N * sizeof *X_high);
	for (int i=0; i < N; i++){
		X_high[i] = malloc((nt+1) * sizeof *X_high[i]);
	}
	V_high = malloc(N * sizeof *V_high);
	for (int i=0; i < N; i++){
		V_high[i] = malloc((nt+1) * sizeof *V_high[i]);
	}
	
	for (i=0; i<N; i++) {
	    X_low[i][0] = x0;
	    V_low[i][0] = v0;
	    X_high[i][0] = x0;
	    V_high[i][0] = v0;
	}
    
    gsl_rng * rg = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rg, 12*12);
    
    for (int i=0; i<N; i++) {
    	brownian_verlet(nt, eta_low, T, m, omega0, X_low[i], V_low[i], dt, rg);
		brownian_verlet(nt, eta_high, T, m, omega0, X_high[i], V_high[i], dt, rg);
    }
    
    printf("Brownian done\n");
    
    double *X_mean_low; double *V_mean_low; double *X_mean_high; double *V_mean_high;
    double *X_var_low; double *V_var_low; double *X_var_high; double *V_var_high;
    
    X_mean_low = malloc((nt+1) * sizeof(double));
    V_mean_low = malloc((nt+1) * sizeof(double));
    X_mean_high = malloc((nt+1) * sizeof(double));
    V_mean_high = malloc((nt+1) * sizeof(double));
    X_var_low = malloc((nt+1) * sizeof(double));
    V_var_low = malloc((nt+1) * sizeof(double));
    X_var_high = malloc((nt+1) * sizeof(double));
    V_var_high = malloc((nt+1) * sizeof(double));
    
    evalMean(X_low, X_mean_low, N, nt);
    evalVariance(X_low, X_var_low, X_mean_low, N, nt);
    printf("X_low done\n");
    evalMean(V_low, V_mean_low, N, nt);
    evalVariance(V_low, V_var_low, V_mean_low, N, nt);
    printf("V_low done\n");
    
    evalMean(X_high, X_mean_high, N, nt);
    evalVariance(X_high, X_var_high, X_mean_high, N, nt);
    printf("X_high done\n");
    evalMean(V_high, V_mean_high, N, nt);
    evalVariance(V_high, V_var_high, V_mean_high, N, nt);
    printf("V_high done\n");
    
    time = malloc((nt+1) * sizeof(double));
	arange(time, 0, nt+1, dt);

    saveDataToFile("2/meanposition_tau147.3e-6.csv", X_mean_low, time, nt+1, 1);
    saveDataToFile("2/varposition_tau147.3e-6.csv", X_var_low, time, nt+1, 1);
    
    saveDataToFile("2/meanvelocity_tau147.3e-6.csv", V_mean_low, time, nt+1, 1);
    saveDataToFile("2/varvelocity_tau147.3e-6.csv", V_var_low, time, nt+1, 1);
    
    saveDataToFile("2/meanposition_tau48.5e-6.csv", X_mean_high, time, nt+1, 1);
    saveDataToFile("2/varposition_tau48.5e-6.csv", X_var_high, time, nt+1, 1);
    
    saveDataToFile("2/meanvelocity_tau48.5e-6.csv", V_mean_high, time, nt+1, 1);
    saveDataToFile("2/varvelocity_tau48.5e-6.csv", V_var_high, time, nt+1, 1);
    
    free(time);
    free(X_mean_low);
    free(V_mean_low);
    free(X_mean_high);
    free(V_mean_high);
    free(X_var_low);
    free(V_var_low);
    free(X_var_high);
    free(V_var_high);
    for (i=0; i < N; i++){
		free(X_low[i]);
	}
	free(X_low);
    for (i=0; i < N; i++){
		free(V_low[i]);
	}
	free(V_low);
    for (i=0; i < N; i++){
		free(X_high[i]);
	}
	free(X_high);
    for (i=0; i < N; i++){
		free(V_high[i]);
	}
	free(V_high);
}

void runtask3() {
    double dt = 1e-6; double t_max = 4e-3; int nt = (int) (t_max/dt); 
    double tau_low = 147.3e-6; double tau_high = 48.5e-6; double eta_low = 1/tau_low; double eta_high = 1/tau_high; double omega0 = 3.1e3*2*PI;
    double r = 2.79e-6/2; double rho = 2.65e3; double m = (4.0/3.0 * PI * pow(r,3)) * rho; double T = 297.0; 
    double x0 = 0.1e-6; double v0 = 2e-3; 
    double **X_low; double **V_low; double **X_high; double **V_high; int N = 10000; int i;
    
    X_low = malloc(N * sizeof *X_low);
	for (i=0; i < N; i++){
		X_low[i] = malloc((nt+1) * sizeof *X_low[i]);
	}
	V_low = malloc(N * sizeof *V_low);
	for (int i=0; i < N; i++){
		V_low[i] = malloc((nt+1) * sizeof *V_low[i]);
	}
    X_high = malloc(N * sizeof *X_high);
	for (int i=0; i < N; i++){
		X_high[i] = malloc((nt+1) * sizeof *X_high[i]);
	}
	V_high = malloc(N * sizeof *V_high);
	for (int i=0; i < N; i++){
		V_high[i] = malloc((nt+1) * sizeof *V_high[i]);
	}
	
	for (i=0; i<N; i++) {
	    X_low[i][0] = x0;
	    V_low[i][0] = v0;
	    X_high[i][0] = x0;
	    V_high[i][0] = v0;
	}
    
    gsl_rng * rg = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rg, 12*12);
    
    for (int i=0; i<N; i++) {
    	brownian_verlet(nt, eta_low, T, m, omega0, X_low[i], V_low[i], dt, rg);
		brownian_verlet(nt, eta_high, T, m, omega0, X_high[i], V_high[i], dt, rg);
    }
    
    printf("Brownian done\n");
    
    double *X_mean_low; double *V_mean_low; double *X_mean_high; double *V_mean_high;
    double *X_var_low; double *V_var_low; double *X_var_high; double *V_var_high;
    
    X_mean_low = malloc((nt+1) * sizeof(double));
    V_mean_low = malloc((nt+1) * sizeof(double));
    X_mean_high = malloc((nt+1) * sizeof(double));
    V_mean_high = malloc((nt+1) * sizeof(double));
    X_var_low = malloc((nt+1) * sizeof(double));
    V_var_low = malloc((nt+1) * sizeof(double));
    X_var_high = malloc((nt+1) * sizeof(double));
    V_var_high = malloc((nt+1) * sizeof(double));
    
    evalMean(X_low, X_mean_low, N, nt);
    evalVariance(X_low, X_var_low, X_mean_low, N, nt);
    printf("X_low done\n");
    evalMean(V_low, V_mean_low, N, nt);
    evalVariance(V_low, V_var_low, V_mean_low, N, nt);
    printf("V_low done\n");
    
    evalMean(X_high, X_mean_high, N, nt);
    evalVariance(X_high, X_var_high, X_mean_high, N, nt);
    printf("X_high done\n");
    evalMean(V_high, V_mean_high, N, nt);
    evalVariance(V_high, V_var_high, V_mean_high, N, nt);
    printf("V_high done\n");
    
    double t = 2e-3; 
    int it = (int) (t/t_max* ((double)nt)); double dx = 1e-9; double dv = 5e-5;
    printf("%i\n", it);
    double x_max_low = X_mean_low[it] + 5*sqrt(X_var_low[it]); double x_min_low = X_mean_low[it] - 5*sqrt(X_var_low[it]);
    double v_max_low = V_mean_low[it] + 5*sqrt(V_var_low[it]); double v_min_low = V_mean_low[it] - 5*sqrt(V_var_low[it]);
    double x_max_high = X_mean_high[it] + 5*sqrt(X_var_high[it]); double x_min_high = X_mean_high[it] - 5*sqrt(X_var_high[it]);
    double v_max_high = V_mean_high[it] + 5*sqrt(V_var_high[it]); double v_min_high = V_mean_high[it] - 5*sqrt(V_var_high[it]);
	
    int nx_low = (int) ((x_max_low-x_min_low)/dx); int nv_low = (int) ((v_max_low-v_min_low)/dv);
    int nx_high = (int) ((x_max_high-x_min_high)/dx); int nv_high = (int) ((v_max_high-v_min_high)/dv);
	
    double *rhoX_low; double *rhoV_low;
	rhoX_low = malloc(nx_low * sizeof(double));
	rhoV_low = malloc(nv_low * sizeof(double));
	
	double *rhoX_high; double *rhoV_high;
	rhoX_high = malloc(nx_high * sizeof(double));
	rhoV_high = malloc(nv_high * sizeof(double));
	
	evalDensityFunctions(X_low, V_low, rhoX_low, rhoV_low, it, dx, dv, nx_low, nv_low, x_max_low, x_min_low, v_max_low, v_min_low, N);
	printf("rho_low done\n");
	evalDensityFunctions(X_high, V_high, rhoX_high, rhoV_high, it, dx, dv, nx_high, nv_high, x_max_high, x_min_high, v_max_high, v_min_high, N);
	printf("rho_high done\n");
	
	double *x_low; double *v_low; 
	x_low = malloc(nx_low * sizeof(double));
	v_low = malloc(nv_low * sizeof(double));
	arange(x_low, x_min_low, nx_low, dx);
	arange(v_low, v_min_low, nv_low, dv);
    saveDataToFile("3/rhoX_tau147.3e-6_t2.csv", rhoX_low, x_low, nx_low, 1);
    saveDataToFile("3/rhoV_tau147.3e-6_t2.csv", rhoV_low, v_low, nv_low, 1);
    free(x_low);
    free(v_low);

	double *x_high; double *v_high; 
	x_high = malloc(nx_high * sizeof(double));
	v_high = malloc(nv_high * sizeof(double));
	arange(x_high, x_min_high, nx_high, dx);
	arange(v_high, v_min_high, nv_high, dv);
    saveDataToFile("3/rhoX_tau48.5e-6_t2.csv", rhoX_high, x_high, nx_high, 1);
    saveDataToFile("3/rhoV_tau48.5e-6_t2.csv", rhoV_high, v_high, nv_high, 1);
    free(x_high);
    free(v_high);
    
    
    
    free(rhoX_low);
    free(rhoV_low);
    free(rhoX_high);
    free(rhoV_high);
    free(X_mean_low);
    free(V_mean_low);
    free(X_mean_high);
    free(V_mean_high);
    free(X_var_low);
    free(V_var_low);
    free(X_var_high);
    free(V_var_high);
    for (i=0; i < N; i++){
		free(X_low[i]);
	}
	free(X_low);
    for (i=0; i < N; i++){
		free(V_low[i]);
	}
	free(V_low);
    for (i=0; i < N; i++){
		free(X_high[i]);
	}
	free(X_high);
    for (i=0; i < N; i++){
		free(V_high[i]);
	}
	free(V_high);
}


int main() {
	//runtask2trajectory();
	//runtask2mean();
	runtask3();
    return 0;
}