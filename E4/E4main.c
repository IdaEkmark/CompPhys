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

void evalVelocityCorrelationStandard(int n_timesteps, double *vel, double *vcf, int vcf_intervall) {
	double vcf_time;
	int time = 0;
	while (time + vcf_intervall < n_timesteps + 1) {
		vcf_time = 0;
		for (int time_sum = time + 1; time_sum <= time + vcf_intervall; time_sum++ ) {
			vcf_time += vel[time_sum] * vel[time_sum - time]/vcf_intervall;
		}
		vcf[time] = vcf_time;
		time++;
	}
}

void runtask1() {
    double dt1 = 5e-6; double dt2 = 1e-6; double t_max = 4e-3; int nt1 = (int) (t_max/dt1); int nt2 = (int) (t_max/dt2); 
    double tau_low = 147.3e-6; double tau_high = 48.5e-6; double eta_low = 1/tau_low; double eta_high = 1/tau_high; double omega0 = 3.1e3*2*PI;
    double r = 2.79e-6/2; double rho = 2.65e3; double m = (4.0/3.0 * PI * pow(r,3)) * rho; double T = 297.0; double x0 = 0; double v0 = 0; 
    double *X1; double *V1; double *X2; double *V2; double *X3; double *V3; double *X4; double *V4; double *time;

    X1 = malloc(nt1 * sizeof(double));
    V1 = malloc(nt1 * sizeof(double));
    X2 = malloc(nt1 * sizeof(double));
    V2 = malloc(nt1 * sizeof(double));
    X3 = malloc(nt2 * sizeof(double));
    V3 = malloc(nt2 * sizeof(double));
    X4 = malloc(nt2 * sizeof(double));
    V4 = malloc(nt2 * sizeof(double));

    X1[0] = x0;
    V1[0] = v0;
    X2[0] = x0;
    V2[0] = v0;
    X3[0] = x0;
    V3[0] = v0;
    X4[0] = x0;
    V4[0] = v0;
    
    gsl_rng * rg = gsl_rng_alloc(gsl_rng_mt19937);

    brownian_verlet(nt1, eta_low, T, m, omega0, X1, V1, dt1, rg);
    brownian_verlet(nt1, eta_high, T, m, omega0, X2, V2, dt1, rg);

    brownian_verlet(nt2, eta_low, T, m, omega0, X3, V3, dt2, rg);
    brownian_verlet(nt2, eta_high, T, m, omega0, X4, V4, dt2, rg);

    time = malloc((nt1+1) * sizeof(double));
	arange(time, -2e-3, nt1+1, dt1);

    saveDataToFile("1/position_tau147.3e-6_dt5e-6.csv", X1, time, nt1+1, 1);
    saveDataToFile("1/velocity_tau147.3e-6_dt5e-6.csv", V1, time, nt1+1, 1);
    saveDataToFile("1/position_tau48.5e-6_dt5e-6.csv", X2, time, nt1+1, 1);
    saveDataToFile("1/velocity_tau48.5e-6_dt5e-6.csv", V2, time, nt1+1, 1);

    free(time);
    time = malloc((nt2+1) * sizeof(double));
	arange(time, -2e-3, nt2+1, dt2);

    saveDataToFile("1/position_tau147.3e-6_dt1e-6.csv", X3, time, nt2+1, 1);
    saveDataToFile("1/velocity_tau147.3e-6_dt1e-6.csv", V3, time, nt2+1, 1);
    saveDataToFile("1/position_tau48.5e-6_dt1e-6.csv", X4, time, nt2+1, 1);
    saveDataToFile("1/velocity_tau48.5e-6_dt1e-6.csv", V4, time, nt2+1, 1);

    free(time);

    free(X1);
    free(V1);
    free(X2);
    free(V2);
    free(X3);
    free(V3);
    free(X4);
    free(V4);
}

void runtask2a() {
    int ntNew = 41; double dtau = 50e-6;
    double time_array[ntNew];
    double vDataLow50[ntNew];
    double vDataHigh50[ntNew];
    read_data("1/velocity_tau147.3e-6_dt1e-6.csv", time_array, vDataLow50, 2000, 50);
    read_data("1/velocity_tau48.5e-6_dt1e-6.csv", time_array, vDataHigh50, 2000, 50);

    /*
     * Construct array with frequencies
     */ 
    double frequencies[ntNew];
    for(int i = 0; i < ntNew; i++){
	    frequencies[i] = i / (dtau * (ntNew-1)) - 1.0 / (2.0 * dtau);
    }

    /*
     * Do the fft
     */
    double fftd_dataLow50[ntNew]; double fftd_dataHigh50[ntNew];
    powerspectrum(vDataLow50, fftd_dataLow50, ntNew);
    powerspectrum_shift(fftd_dataLow50, ntNew);
    powerspectrum(vDataHigh50, fftd_dataHigh50, ntNew);
    powerspectrum_shift(fftd_dataHigh50, ntNew);

    saveDataToFile("2/power_147.3e-6_dtau50dt.csv", fftd_dataLow50, frequencies, ntNew, 1);
    saveDataToFile("2/power_48.5e-6_dtau50dt.csv", fftd_dataHigh50, frequencies, ntNew, 1);

    ntNew = 81; dtau = 25e-6;
    double time_array2[ntNew];
    double vDataLow25[ntNew];
    double vDataHigh25[ntNew];
    read_data("1/velocity_tau147.3e-6_dt1e-6.csv", time_array2, vDataLow25, 2000, 25);
    read_data("1/velocity_tau48.5e-6_dt1e-6.csv", time_array2, vDataHigh25, 2000, 25);

    /*
     * Construct array with frequencies
     */ 
    double frequencies2[ntNew];
    for(int i = 0; i < ntNew; i++){
	    frequencies2[i] = i / (dtau * (ntNew-1)) - 1.0 / (2.0 * dtau);
    }

    /*
     * Do the fft
     */
    double fftd_dataLow25[ntNew]; double fftd_dataHigh25[ntNew];
    powerspectrum(vDataLow25, fftd_dataLow25, ntNew);
    powerspectrum_shift(fftd_dataLow25, ntNew);
    powerspectrum(vDataHigh25, fftd_dataHigh25, ntNew);
    powerspectrum_shift(fftd_dataHigh25, ntNew);

    saveDataToFile("2/power_147.3e-6_dtau25dt.csv", fftd_dataLow25, frequencies2, ntNew, 1);
    saveDataToFile("2/power_48.5e-6_dtau25dt.csv", fftd_dataHigh25, frequencies2, ntNew, 1);
}

void runtask2b() {
    double dt = 1e-6; double t_max = 10e-3; int nt = (int) (t_max/dt);
    double tau_low = 147.3e-6; double tau_high = 48.5e-6; double eta_low = 1/tau_low; double eta_high = 1/tau_high; double omega0 = 3.1e3*2*PI;
    double r = 2.79e-6/2; double rho = 2.65e3; double m = (4.0/3.0 * PI * pow(r,3)) * rho; double T = 297.0; double x0 = 0; double v0 = 0; 
    double *X1; double *V1; double *X2; double *V2; double *time;

    X1 = malloc(nt * sizeof(double));
    V1 = malloc(nt * sizeof(double));
    X2 = malloc(nt * sizeof(double));
    V2 = malloc(nt * sizeof(double));

    X1[0] = x0;
    V1[0] = v0;
    X2[0] = x0;
    V2[0] = v0;
    
    gsl_rng * rg = gsl_rng_alloc(gsl_rng_mt19937);
    
    brownian_verlet(nt, eta_low, T, m, omega0, X1, V1, dt, rg);
    brownian_verlet(nt, eta_high, T, m, omega0, X2, V2, dt, rg);

    time = malloc((nt+1) * sizeof(double));
	arange(time, -3e-3, nt+1, dt);

    saveDataToFile("2/position_tau147.3e-6_dt1e-6_Long.csv", X1, time, nt+1, 1);
    saveDataToFile("2/velocity_tau147.3e-6_dt1e-6_Long.csv", V1, time, nt+1, 1);
    saveDataToFile("2/position_tau48.5e-6_dt1e-6_Long.csv", X2, time, nt+1, 1);
    saveDataToFile("2/velocity_tau48.5e-6_dt1e-6_Long.csv", V2, time, nt+1, 1);

    // Compute 4 different power spectra
    int ntNew = 321; double dtau = 25e-6; int startind = 2000;
    double vDataLow25[ntNew];
    double vDataHigh25[ntNew];

    for (int i=0; i<ntNew; i++) {
        vDataLow25[i] = V1[startind + 25*i];
        vDataLow25[i] = V1[startind + 25*i];
    }

    /*
     * Construct array with frequencies
     */
    double frequencies[ntNew];
    for(int i = 0; i < ntNew; i++){
	    frequencies[i] = i / (dtau * (ntNew-1)) - 1.0 / (2.0 * dtau);
    }

    /*
     * Do the fft
     */
    double fftd_dataLow25[ntNew]; double fftd_dataHigh25[ntNew];
    powerspectrum(vDataLow25, fftd_dataLow25, ntNew);
    powerspectrum_shift(fftd_dataLow25, ntNew);
    powerspectrum(vDataHigh25, fftd_dataHigh25, ntNew);
    powerspectrum_shift(fftd_dataHigh25, ntNew);

    free(time);

    free(X1);
    free(V1);
    free(X2);
    free(V2);
}

void runtask3(){
	double dt = 1e-6; double t_max = 8e-3; int nt = (int) (t_max/dt);
    double vDataLow[nt];
    double vDataHigh[nt];
    double time_array[nt];
    read_data("1/velocity_tau147.3e-6_dt1e-6.csv", time_array, vDataLow, 2000, 1);
    read_data("1/velocity_tau48.5e-6_dt1e-6.csv", time_array, vDataHigh, 2000, 1);
    
    double *VCF_high; double *VCF_low; int vcf_intervall = 5000; double *time;
    VCF_high = malloc((nt+1-vcf_intervall) * sizeof(double));
    VCF_low = malloc((nt+1-vcf_intervall) * sizeof(double));
    time = malloc((nt+1-vcf_intervall) * sizeof(double));
	arange(time, 0.0, nt+1-vcf_intervall, dt);
	
    evalVelocityCorrelationStandard(nt, vDataHigh, VCF_high, vcf_intervall);
    evalVelocityCorrelationStandard(nt, vDataLow, VCF_low, vcf_intervall);
    
    saveDataToFile("3/velcorr_147.3e-6.csv", VCF_low, time, nt+1-vcf_intervall, 1);
    saveDataToFile("3/velcorr_48.5e-6.csv", VCF_high, time, nt+1-vcf_intervall, 1);
    
	free(VCF_high);
	free(VCF_low);
	free(time);
    
}

int main() {
    //runtask1();
    runtask2a();
    //runtask2b();
	//runtask3();
    return 0;
}