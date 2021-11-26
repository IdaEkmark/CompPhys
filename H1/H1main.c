/*
 H1main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "H1lattice.h"
#include "H1potential.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>

#define N_ATOMS 256
#define K_B 8.6e-5
#define KAPPA_T 2.10
#define PI 3.141592653589
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
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
	    fprintf(fp, "%f, %f\n", tvals[i], yvals[i]);
    }
    fclose(fp);
}

void saveQPtoFile(char *fname, double *time_array,
		double (*pos)[N_ATOMS][3], double (*mom)[N_ATOMS][3], int n_timesteps, int n_particles, int n_p_skip, int n_t_skip) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "time");
    for (int i = 0; i < n_particles; i = i+n_p_skip) {
    	for (int j = 0; j < 3; j++) {
    		fprintf(fp, ", Q%i%i", i, j);
    	}
    	for (int j = 0; j < 3; j++) {
			fprintf(fp, ", P%i%i", i, j);
		}
    }
    fprintf(fp, "\n");

    for(int t = 0; t < n_timesteps; t = t + n_t_skip){
	    fprintf(fp, "%f", time_array[t]);
	    
	    for (int i = 0; i < n_particles; i = i+n_p_skip) {
			for (int j = 0; j < 3; j++) {
				fprintf(fp, ", %f", pos[t][i][j]);
			}
			for (int j = 0; j < 3; j++) {
				fprintf(fp, ", %f", mom[t][i][j]);
			}
		}
        fprintf(fp, "\n");
    }
    fclose(fp);
}

double get_E_kin(double momenta[][3], double n_particles, double mass) 
{
	double E_kin = 0;
	double mass_inv = 1/mass;
	
	for (int p = 0; p < n_particles; p++) {
		for (int i = 0; i < 3; i++) {
			E_kin += 0.5 * momenta[p][i] * momenta[p][i] * mass_inv; 
		}
	}
	return E_kin;
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
void velocity_verlet(int n_timesteps, int n_particles, double a0, double (*pos)[N_ATOMS][3], double (*mom)[N_ATOMS][3],
		double dt, double mass, int N)
{
    double q[n_particles][3];
    double p[n_particles][3];
    double f[n_particles][3];

    int i; int j;
    double mass_inv = 1/mass;

    
    for (i=0; i< n_particles; i++) {
        for (j=0; j<3; j++) {
            q[i][j] = pos[0][i][j];
            p[i][j] = mom[0][i][j];
        }
    }
    get_forces_AL(f,q,N*a0,n_particles);
    
    for (int t = 1; t < n_timesteps + 1; t++) {
        // v(t+dt/2)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
        
        // q(t+dt)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                q[i][j] += dt * p[i][j]*mass_inv;
            }
        }
        
        // a(t+dt)
        get_forces_AL(f,q,N*a0,n_particles);
        
        // v(t+dt)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
		
        // Save the displacement of the three atoms
        for (i = 0; i<n_particles; i++) {
            for (j = 0; j<3; j++) {
                pos[t][i][j] = q[i][j];
                mom[t][i][j] = p[i][j];
            }
        }
    }
}

double getInstantaneousTemperature(double momentum[][3], int n_particles, double mass) {
	double E_kin = get_E_kin(momentum, n_particles, mass);
	double T_inst = 2 * E_kin / (3 * n_particles * K_B);
	return T_inst;
}

double getInstantaneousPressure(double positions[][3], double momentum[][3], int n_particles, double cell_length, double mass){
	double E_kin = get_E_kin(momentum, n_particles, mass);
	double virial = get_virial_AL(positions, cell_length, n_particles);
	double P_inst = (2 * E_kin + virial)/(3.0 * cell_length * cell_length * cell_length);
	return P_inst;
}

double velocity_verlet_equi(int n_timesteps, int n_particles, double a0, double dt, 
		double mass, int N, double T_eq, double P_eq, double tau_T, double tau_P, 
		double (*pos)[N_ATOMS][3], double (*mom)[N_ATOMS][3],  double *a0_vec)
{
    double q[n_particles][3];
    double p[n_particles][3];
    double f[n_particles][3];
    		
    int i; int j;
    double mass_inv = 1/mass; double tau_T_inv = 1/tau_T; double tau_P_inv = 1/tau_P; 
    double T_inst; double alpha_T; double P_inst; double alpha_P;
    
    for (i=0; i< n_particles; i++) {
        for (j=0; j<3; j++) {
            q[i][j] = pos[0][i][j];
            p[i][j] = mom[0][i][j];
        }
    }
    
    
    get_forces_AL(f,q,N*a0,n_particles);
    
    for (int t = 1; t < n_timesteps + 1; t++) {
        // v(t+dt/2)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
        
        // q(t+dt)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                q[i][j] += dt * p[i][j]*mass_inv;
            }
        }
        
        // a(t+dt)
        get_forces_AL(f,q,N*a0,n_particles);
        
        // v(t+dt)
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
		
        T_inst = getInstantaneousTemperature(p, n_particles, mass);
        alpha_T = 1 + 2 * dt * tau_T_inv * ( T_eq - T_inst ) / T_inst; 
        
        P_inst = getInstantaneousPressure(q, p, n_particles, N*a0, mass);
        alpha_P = 1.0 - KAPPA_T * dt * tau_P_inv * ( P_eq - P_inst );
        
        for (i = 0; i < n_particles; i++) {
        	for (j = 0; j < 3; j++) {
        		p[i][j] *= sqrt(alpha_T);
        		q[i][j] *= cbrt(alpha_P);
        	}
        }
        a0 *= cbrt(alpha_P);
        
        // Save the displacement of the three atoms
        for (i = 0; i<n_particles; i++) {
            for (j = 0; j<3; j++) {
                pos[t][i][j] = q[i][j];
                mom[t][i][j] = p[i][j];
            }
        }
        a0_vec[t] = a0;
    }
    return a0;
}

double evalSquaredDisplacementAtTime(int n_particles, double (*pos)[N_ATOMS][3], int time, int time_sum) {
	double sd_time = 0;
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < 3; j++) {
			sd_time += (pos[time_sum][i][j] - pos[time_sum - time][i][j]) * (pos[time_sum][i][j] - pos[time_sum - time][i][j]) / n_particles; 
		}
	}
	return sd_time;
}

void evalMeanSquaredDisplacement(int n_timesteps, int n_particles, double (*pos)[N_ATOMS][3], double *msd, int msd_intervall) {
	double msd_time;
	int time = 0;
	while (time + msd_intervall < n_timesteps + 1) {
		msd_time = 0;
		for (int time_sum = time + 1; time_sum <= time + msd_intervall; time_sum++ ) {
			msd_time += evalSquaredDisplacementAtTime(n_particles, pos, time, time_sum)/msd_intervall;
		}
		msd[time] = msd_time;
		time++;
	}
	
}

double evalVelocityCorrelationAtTime(int n_particles, double (*mom)[N_ATOMS][3], double mass, int time, int time_sum) {
	double vc_time = 0;
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < 3; j++) {
			vc_time += mom[time_sum][i][j] * mom[time_sum - time][i][j] / (n_particles * mass * mass); 
		}
	}
	return vc_time;
}

void evalVelocityCorrelationStandard(int n_timesteps, int n_particles, double (*mom)[N_ATOMS][3], double mass,
		double *vcf, int vcf_intervall) {
	double vcf_time;
	int time = 0;
	while (time + vcf_intervall < n_timesteps + 1) {
		vcf_time = 0;
		for (int time_sum = time + 1; time_sum <= time + vcf_intervall; time_sum++ ) {
			vcf_time += evalVelocityCorrelationAtTime(n_particles, mom, mass, time, time_sum)/vcf_intervall;
		}
		vcf[time] = vcf_time;
		time++;
	}
	
}

void evalVelocityCorrelationFast(int n_timesteps, int n_particles, double (*mom)[N_ATOMS][3], double mass,
		double *vcf, int vcf_intervall) {
	int t; int i;
	double velExtended_x[2*vcf_intervall]; 
	double velExtended_y[2*vcf_intervall];
	double velExtended_z[2*vcf_intervall];
	double invFourierVelExtended_x[2*(2*vcf_intervall)]; 
	double invFourierVelExtended_y[2*(2*vcf_intervall)];
	double invFourierVelExtended_z[2*(2*vcf_intervall)];
	double fourier_x[2*vcf_intervall]; 
	double fourier_y[2*vcf_intervall];
	double fourier_z[2*vcf_intervall];
	double normfactor = (double) vcf_intervall/((double)n_timesteps * (double)n_particles) ;
	
	/*Declare wavetable and workspace for fft*/
	gsl_fft_complex_wavetable *comp;
	gsl_fft_complex_workspace *work_c;
	
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work_r;

	/*Allocate space for wavetable and workspace for fft*/
	work_c = gsl_fft_complex_workspace_alloc(vcf_intervall);
	comp = gsl_fft_complex_wavetable_alloc(vcf_intervall);
		
	work_r = gsl_fft_real_workspace_alloc(vcf_intervall);
	real = gsl_fft_real_wavetable_alloc(vcf_intervall);
	
	for(t = 0; t < vcf_intervall; t++) {
		vcf[t] = 0;
	}
	int n_t_start;
	
	for(i = 0; i < n_particles; i++) {
		n_t_start = 0;
		while(n_t_start + vcf_intervall <= n_timesteps) {
			for(t = 0; t < vcf_intervall; t++) {
				velExtended_x[t] = mom[t + n_t_start][i][0]/mass;
				velExtended_y[t] = mom[t + n_t_start][i][1]/mass;
				velExtended_z[t] = mom[t + n_t_start][i][2]/mass;
			}
			for(t = vcf_intervall; t < 2*vcf_intervall; t++) {
				velExtended_x[t] = 0;
				velExtended_y[t] = 0;
				velExtended_z[t] = 0;
			}
			
			/*make copy of data to avoid messing with data in the transform*/
			for (t = 0; t < 2*vcf_intervall; t++)	{
				REAL(invFourierVelExtended_x, t) = velExtended_x[t];
				IMAG(invFourierVelExtended_x, t) = 0.0;
				REAL(invFourierVelExtended_y, t) = velExtended_y[t];
				IMAG(invFourierVelExtended_y, t) = 0.0;
				REAL(invFourierVelExtended_z, t) = velExtended_z[t];
				IMAG(invFourierVelExtended_z, t) = 0.0;
			}
	
			/*Do the fft*/
			gsl_fft_complex_inverse(invFourierVelExtended_x, 1, vcf_intervall, comp, work_c);
			gsl_fft_complex_inverse(invFourierVelExtended_y, 1, vcf_intervall, comp, work_c);
			gsl_fft_complex_inverse(invFourierVelExtended_z, 1, vcf_intervall, comp, work_c);
			 
			for (t = 0; t < 2*vcf_intervall; t++)	{
				fourier_x[t] = REAL(invFourierVelExtended_x, t)*REAL(invFourierVelExtended_x, t) + 
								IMAG(invFourierVelExtended_x, t)*IMAG(invFourierVelExtended_x, t);
				fourier_y[t] = REAL(invFourierVelExtended_y, t)*REAL(invFourierVelExtended_y, t) + 
								IMAG(invFourierVelExtended_y, t)*IMAG(invFourierVelExtended_y, t);

				fourier_z[t] = REAL(invFourierVelExtended_z, t)*REAL(invFourierVelExtended_z, t) + 
								IMAG(invFourierVelExtended_z, t)*IMAG(invFourierVelExtended_z, t);
			}
			gsl_fft_real_transform(fourier_x, 1, vcf_intervall, real, work_r);
			gsl_fft_real_transform(fourier_y, 1, vcf_intervall, real, work_r);
			gsl_fft_real_transform(fourier_z, 1, vcf_intervall, real, work_r);
			
			for(t = 0; t < vcf_intervall; t++) {
				vcf[t] += (fourier_x[2*t+1] + fourier_y[2*t+1] + fourier_z[2*t+1])*normfactor;
			}
			
			n_t_start += vcf_intervall;
			
		}
	}
}

void powerspectrum(double *data, double *powspec, int n, double dt, double *omegas, int n_omegas, double domega) /* input data, output powspec_data, number of timesteps */
{
	double omega = 0;
	for(int i = 0; i < n_omegas; i++){
		omegas[i] = omega;
		powspec[i] = 0;
		for(int j = 0; j < n; j++) {
			powspec[i] += data[j] * cos(omega * data[j]); 
		}
		powspec[i] *= 2*dt;
		omega += domega;
	}
}

void runTask1() {
	double *a0s; int N = 4; int natoms = N_ATOMS; double X[natoms][3];
	double *E_pots; int na = 11; double da = 0.1/((double)na-1);

	a0s = malloc(na * sizeof(double));
	arange(a0s, 4, na, da);

	E_pots = malloc(na * sizeof(double));

	for (int i=0; i < na; i++) {
		init_fcc(X, N, a0s[i]);
		E_pots[i] = get_energy_AL(X, N * a0s[i] , natoms);
	}

	saveDataToFile("1/Epots_a0.csv", E_pots, a0s, na, 1);

	free(E_pots);
	free(a0s);
}

void runTask2(double dt) {
	double a0 = 4.03; int N = 4; int n_t = 10000; int natoms = N_ATOMS; double mass = 27.0 / 9649.0;
	double (*positions)[natoms][3]; double (*momenta)[natoms][3]; double standardpositions[natoms][3];
	double *E_pot; double *E_kin; double *time;
	int i; int j; int t;

    positions = malloc((n_t+1) * sizeof *positions);
    momenta = malloc((n_t+1) * sizeof *momenta);
    E_pot = malloc((n_t + 1) * sizeof(double));
    E_kin = malloc((n_t + 1) * sizeof(double));
    
	init_fcc(standardpositions, N, a0);
	
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions[0][i][j] = standardpositions[i][j] + a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
            momenta[0][i][j] = 0.0;
		}
	}
	
	velocity_verlet(n_t, natoms, a0, positions, momenta, dt, mass, N);
	
	for (t = 0; t < n_t + 1; t++) {
		E_pot[t] = get_energy_AL(positions[t], N*a0, natoms);
		E_kin[t] = get_E_kin(momenta[t], natoms, mass);
	}
	
	time = malloc((n_t+1) * sizeof(double));
	arange(time, 0.0, n_t+1, dt);
	
	if (dt == 0.001) {
		saveDataToFile("2/Epot_dt0.001.csv", E_pot, time, n_t, 1);
		saveDataToFile("2/Ekin_dt0.001.csv", E_kin, time, n_t, 1);
	}
	else if (dt == 0.01) {
		saveDataToFile("2/Epot_dt0.01.csv", E_pot, time, n_t, 1);
		saveDataToFile("2/Ekin_dt0.01.csv", E_kin, time, n_t, 1);
	}
	else if (dt == 0.02) {
		saveDataToFile("2/Epot_dt0.02.csv", E_pot, time, n_t, 1);
		saveDataToFile("2/Ekin_dt0.02.csv", E_kin, time, n_t, 1);
	}
	
	
    free(positions);
    free(momenta);
    free(E_pot);
    free(E_kin);
}

void runTask3() {
	double a0 = 4.03; double mass = 27.0 / 9649.0;
	int N = 4; int n_t_equi = 10000; int n_t = 30000; double dt = 1e-3; int natoms = N_ATOMS; int n_p_skip = 85; int n_t_skip = 1;
	double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
	double *temperature; double *pressure; double *time;  double *timeTP; double *timeA0;
	double T_eq = 500 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
	int i; int j; int t;

	positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
	momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
	a0_equi = malloc((n_t_equi) * sizeof(double));
	temperature = malloc((n_t_equi + n_t) * sizeof(double));
	pressure = malloc((n_t_equi + n_t) * sizeof(double));
	
	init_fcc(positions_equi[0], N, a0);
	
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
			momenta_equi[0][i][j] = 0.0;
		}
	}
		
	velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
			positions_equi, momenta_equi, a0_equi);
	
	for ( t = 0; t < n_t_equi; t++ ) {
		temperature[t] = getInstantaneousTemperature(momenta_equi[t+1], natoms, mass);
		pressure[t] = getInstantaneousPressure(positions_equi[t+1], momenta_equi[t+1], natoms, N*a0_equi[t], mass);
	}
	
	double (*positions)[natoms][3]; double (*momenta)[natoms][3];
	
	positions = malloc((n_t+1) * sizeof *positions);
	momenta = malloc((n_t+1) * sizeof *momenta);
			
	
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions[0][i][j] = positions_equi[n_t_equi][i][j];
			momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
		}
	}
	
	free(positions_equi);
	free(momenta_equi);
	
	velocity_verlet(n_t, natoms, a0_equi[n_t_equi-1], positions, momenta, dt, mass, N);
	
	for ( t = 0; t < n_t; t++ ) {
		temperature[t + n_t_equi] = getInstantaneousTemperature(momenta[t+1], natoms, mass);
		pressure[t + n_t_equi] = getInstantaneousPressure(positions[t+1], momenta[t+1], natoms, N*a0_equi[n_t_equi-1], mass);
	}
	
	time = malloc((n_t+1) * sizeof(double));
	arange(time, 0.0, n_t+1, dt);
	
	timeTP = malloc((n_t_equi+n_t) * sizeof(double));
	arange(timeTP, 0.0, (n_t_equi+n_t), dt);
	
	timeA0 = malloc((n_t_equi) * sizeof(double));
	arange(timeA0, 0.0, (n_t_equi), dt);
	
	saveQPtoFile("3/QP_dt0.001.csv", time, positions, momenta, n_t, natoms, n_p_skip, n_t_skip);
	saveDataToFile("3/T_dt0.001.csv", temperature, timeTP, n_t_equi+n_t, n_t_skip);
	saveDataToFile("3/P_dt0.001.csv", pressure, timeTP, n_t_equi+n_t, n_t_skip);
	saveDataToFile("3/a0_dt0.001.csv", a0_equi, timeA0, n_t_equi, n_t_skip);
	
	free(positions);
	free(momenta);
	free(temperature);
	free(pressure);
	
	printf("a0 = %.4e\n", a0_equi[n_t_equi-1]);
	free(a0_equi);
	free(time);
	free(timeTP);
	free(timeA0);
}

void runTask4() {
	double a0 = 4.03; double mass = 27.0 / 9649.0;
	int N = 4; int n_t_equi = 10000; int n_t = 30000; double dt = 1e-3; int natoms = N_ATOMS; int n_p_skip = 85; int n_t_skip = 1;
	double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi; double *a0_equi2;
	double *temperature; double *pressure; double *time;  double *timeTP; double *timeA0;
	double T_eq_init = 1000 + 273.15; double T_eq = 700 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
	int i; int j; int t;

	positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
	momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
	a0_equi = malloc((2*n_t_equi + 1) * sizeof(double));
	a0_equi2 = malloc((n_t_equi) * sizeof(double));
	temperature = malloc((2*n_t_equi + n_t) * sizeof(double));
	pressure = malloc((2*n_t_equi + n_t) * sizeof(double));
	
	init_fcc(positions_equi[0], N, a0);
	
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
			momenta_equi[0][i][j] = 0.0;
		}
	}
		
	velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq_init, P_eq, tau_T, tau_P, 
			positions_equi, momenta_equi, a0_equi);
	
	a0_equi[n_t] = a0_equi[n_t-1];
	
	for ( t = 0; t < n_t_equi; t++ ) {
		temperature[t] = getInstantaneousTemperature(momenta_equi[t+1], natoms, mass);
		pressure[t] = getInstantaneousPressure(positions_equi[t+1], momenta_equi[t+1], natoms, N*a0_equi[t], mass);
	}
	
	init_fcc(positions_equi[0], N, a0);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
			momenta_equi[0][i][j] = 0.0;
		}
	}
	
	velocity_verlet_equi(n_t_equi, natoms, a0_equi[n_t_equi - 1], dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
			positions_equi, momenta_equi, a0_equi2);
	
	for (t = 0; t < n_t_equi; t++){
		a0_equi[n_t_equi + 1 + t] = a0_equi2[t];
	}
	
	for ( t = 0; t < n_t_equi; t++ ) {
		temperature[t + n_t_equi] = getInstantaneousTemperature(momenta_equi[t+1], natoms, mass);
		pressure[t + n_t_equi] = getInstantaneousPressure(positions_equi[t+1], momenta_equi[t+1], natoms, N*a0_equi2[t], mass);
	}
	
	
	double (*positions)[natoms][3]; double (*momenta)[natoms][3];
	
	positions = malloc((n_t+1) * sizeof *positions);
	momenta = malloc((n_t+1) * sizeof *momenta);
			
	
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions[0][i][j] = positions_equi[n_t_equi][i][j];
			momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
		}
	}
	
	free(positions_equi);
	free(momenta_equi);
	
	velocity_verlet(n_t, natoms, a0_equi2[n_t_equi - 1], positions, momenta, dt, mass, N);
	
	for ( t = 0; t < n_t + 1; t++ ) {
		temperature[t + 2*n_t_equi] = getInstantaneousTemperature(momenta[t+1], natoms, mass);
		pressure[t + 2*n_t_equi] = getInstantaneousPressure(positions[t+1], momenta[t+1], natoms, N*a0_equi2[n_t_equi-1], mass);
	}
	
	time = malloc((n_t+1) * sizeof(double));
	arange(time, 0.0, n_t+1, dt);
	
	timeTP = malloc((2*n_t_equi+n_t) * sizeof(double));
	arange(timeTP, 0.0, (2*n_t_equi+n_t), dt);
	
	timeA0 = malloc((2*n_t_equi+1) * sizeof(double));
	arange(timeA0, 0.0, (2*n_t_equi+1), dt);
	
	saveQPtoFile("4/QP_dt0.001.csv", time, positions, momenta, n_t, natoms, n_p_skip, n_t_skip);
	saveDataToFile("4/T_dt0.001.csv", temperature, timeTP, 2*n_t_equi+n_t, n_t_skip);
	saveDataToFile("4/P_dt0.001.csv", pressure, timeTP, 2*n_t_equi+n_t, n_t_skip);
	saveDataToFile("4/a0_dt0.001.csv", a0_equi, timeA0, 2*n_t_equi+1, n_t_skip);
	
	free(positions);
	free(momenta);
	free(temperature);
	free(pressure);
	
	//printf("a0 = %.4e\n", a0_equi[n_t_equi-1]);
	free(a0_equi);
	free(time);
	free(timeTP);
	free(timeA0);
}

void runTask5(char phase) {
	if (phase == 's') { // As task 3
		double a0 = 4.03; double mass = 27.0 / 9649.0;
		int N = 4; int n_t_equi = 10000; int n_t = 30000; double dt = 1e-3; int natoms = N_ATOMS;
		double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
		double T_eq = 500 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
		int i; int j; 
		
		positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
		momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
		a0_equi = malloc((n_t_equi) * sizeof(double));
		
		init_fcc(positions_equi[0], N, a0);
		
		gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
				momenta_equi[0][i][j] = 0.0;
			}
		}
			
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);

		
		double (*positions)[natoms][3]; double (*momenta)[natoms][3];
		
		positions = malloc((n_t+1) * sizeof *positions);
		momenta = malloc((n_t+1) * sizeof *momenta);
				
		
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		free(positions_equi);
		free(momenta_equi);
		
		velocity_verlet(n_t, natoms, a0_equi[n_t_equi-1], positions, momenta, dt, mass, N);
		
		free(momenta);
		
		double *MSD; double *time; int msd_intervall = 1000;
				
		MSD = malloc((n_t+1-msd_intervall) * sizeof(double));
		
		evalMeanSquaredDisplacement(n_t, natoms, positions, MSD, msd_intervall);
		time = malloc((n_t+1-msd_intervall) * sizeof(double));
		arange(time, 0.0, n_t+1-msd_intervall, dt);
		
		saveDataToFile("5/MSD_dt0.001_solid.csv", MSD, time,  n_t+1-msd_intervall, 1);
		
		free(positions);
		free(MSD);
		free(time);
	}
	
	else if (phase == 'l') { // As task 4
		double a0 = 4.03; double mass = 27.0 / 9649.0;
		int N = 4; int n_t_equi = 10000; int n_t = 10000; double dt = 1e-3; int natoms = N_ATOMS; 
		double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
		double T_eq_init = 1000 + 273.15; double T_eq = 700 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
		int i; int j;
	
		positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
		momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
		a0_equi = malloc((n_t_equi) * sizeof(double));
		
		init_fcc(positions_equi[0], N, a0);
		
		gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
				momenta_equi[0][i][j] = 0.0;
			}
		}
			
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq_init, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta_equi[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		double (*positions)[natoms][3]; double (*momenta)[natoms][3];
		
		positions = malloc((n_t+1) * sizeof *positions);
		momenta = malloc((n_t+1) * sizeof *momenta);
				
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		free(positions_equi);
		free(momenta_equi);
		
		velocity_verlet(n_t, natoms, a0_equi[n_t_equi - 1], positions, momenta, dt, mass, N);
		
		free(momenta);
				
		double *MSD; double *time; int msd_intervall = 1000;
				
		MSD = malloc((n_t+1-msd_intervall) * sizeof(double));
		
		evalMeanSquaredDisplacement(n_t, natoms, positions, MSD, msd_intervall);
		time = malloc((n_t+1-msd_intervall) * sizeof(double));
		arange(time, 0.0, n_t+1-msd_intervall, dt);
		
		saveDataToFile("5/MSD_dt0.001_liquid.csv", MSD, time,  n_t+1-msd_intervall, 1);
		
		printf("D=%.4e\n", MSD[n_t-msd_intervall]/(6*time[n_t-msd_intervall]));
		
		free(positions);
		free(MSD);
		free(time);
	}
	else {
		printf("Phase should be s for solid or l for liquid");
	}
}

void runTask6(char alg) {
	if (alg == 's') { // As task 4
		double a0 = 4.03; double mass = 27.0 / 9649.0;
		int N = 4; int n_t_equi = 10000; int n_t = 50000; double dt = 1e-3; int natoms = N_ATOMS; 
		double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
		double T_eq_init = 1000 + 273.15; double T_eq = 700 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
		int i; int j;
	
		positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
		momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
		a0_equi = malloc((n_t_equi) * sizeof(double));
		
		init_fcc(positions_equi[0], N, a0);
		
		gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
				momenta_equi[0][i][j] = 0.0;
			}
		}
			
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq_init, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta_equi[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		double (*positions)[natoms][3]; double (*momenta)[natoms][3];
		
		positions = malloc((n_t+1) * sizeof *positions);
		momenta = malloc((n_t+1) * sizeof *momenta);
				
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		free(positions_equi);
		free(momenta_equi);
		
		velocity_verlet(n_t, natoms, a0_equi[n_t_equi - 1], positions, momenta, dt, mass, N);
		
		free(positions);
				
		double *VCF; double *time; int vcf_intervall = 5000;
						
		VCF = malloc((n_t+1-vcf_intervall) * sizeof(double));
		
		evalVelocityCorrelationStandard(n_t, natoms, momenta, mass, VCF, vcf_intervall);
		time = malloc((n_t+1-vcf_intervall) * sizeof(double));
		arange(time, 0.0, n_t+1-vcf_intervall, dt);
		
		saveDataToFile("6/VCF_dt0.001_standard.csv", VCF, time,  n_t+1-vcf_intervall, 10);
		
		free(momenta);
		free(VCF);
		free(time);
	}
	else if (alg == 'f') { // As task 4
		double a0 = 4.03; double mass = 27.0 / 9649.0;
		int N = 4; int n_t_equi = 10000; int n_t = 50000; double dt = 1e-3; int natoms = N_ATOMS;
		double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
		double T_eq_init = 1000 + 273.15; double T_eq = 700 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
		int i; int j;
	
		positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
		momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
		a0_equi = malloc((n_t_equi) * sizeof(double));
		
		init_fcc(positions_equi[0], N, a0);
		
		gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
				momenta_equi[0][i][j] = 0.0;
			}
		}
			
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq_init, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions_equi[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta_equi[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
				positions_equi, momenta_equi, a0_equi);
		
		
		double (*positions)[natoms][3]; double (*momenta)[natoms][3];
		
		positions = malloc((n_t+1) * sizeof *positions);
		momenta = malloc((n_t+1) * sizeof *momenta);
				
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				positions[0][i][j] = positions_equi[n_t_equi][i][j];
				momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
			}
		}
		
		free(positions_equi);
		free(momenta_equi);
		
		velocity_verlet(n_t, natoms, a0_equi[n_t_equi - 1], positions, momenta, dt, mass, N);
		
		free(positions);
		
		double *VCF; double *time; int vcf_intervall = 2000;
								
		VCF = malloc((vcf_intervall) * sizeof(double));
		
		evalVelocityCorrelationFast(n_t, natoms, momenta, mass, VCF, vcf_intervall);
		time = malloc((vcf_intervall) * sizeof(double));
		arange(time, 0.0, vcf_intervall, dt);
		
		saveDataToFile("6/VCF_dt0.01_fast.csv", VCF, time,  vcf_intervall, 1);
		
		free(momenta);
		free(VCF);
		free(time);
	}
	else {
		printf("Phase should be s for standard or f for Fast Correlation Algorithm");
	}
}

void runTask7() {
	double a0 = 4.03; double mass = 27.0 / 9649.0;
	int N = 4; int n_t_equi = 10000; int n_t = 50000; double dt = 1e-3; int natoms = N_ATOMS;
	double (*positions_equi)[natoms][3]; double (*momenta_equi)[natoms][3]; double *a0_equi;
	double T_eq_init = 1000 + 273.15; double T_eq = 700 + 273.15; double P_eq = 6.24e-7; double tau_T = 400 * dt; double tau_P = 400 * dt;
	int i; int j;

	positions_equi = malloc((n_t_equi+1) * sizeof *positions_equi);
	momenta_equi = malloc((n_t_equi+1) * sizeof *momenta_equi);
	a0_equi = malloc((n_t_equi) * sizeof(double));
	
	init_fcc(positions_equi[0], N, a0);
	
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions_equi[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
			momenta_equi[0][i][j] = 0.0;
		}
	}
		
	velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq_init, P_eq, tau_T, tau_P, 
			positions_equi, momenta_equi, a0_equi);
	
	
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions_equi[0][i][j] = positions_equi[n_t_equi][i][j];
			momenta_equi[0][i][j] = momenta_equi[n_t_equi][i][j];
		}
	}
	
	velocity_verlet_equi(n_t_equi, natoms, a0, dt, mass, N, T_eq, P_eq, tau_T, tau_P, 
			positions_equi, momenta_equi, a0_equi);
	
	
	double (*positions)[natoms][3]; double (*momenta)[natoms][3];
	
	positions = malloc((n_t+1) * sizeof *positions);
	momenta = malloc((n_t+1) * sizeof *momenta);
			
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions[0][i][j] = positions_equi[n_t_equi][i][j];
			momenta[0][i][j] = momenta_equi[n_t_equi][i][j];
		}
	}
	
	free(positions_equi);
	free(momenta_equi);
	
	velocity_verlet(n_t, natoms, a0_equi[n_t_equi - 1], positions, momenta, dt, mass, N);
	
	free(positions);
	
	double *VCF; int vcf_intervall = 1000;
								
	VCF = malloc((vcf_intervall) * sizeof(double));
	
	
	evalVelocityCorrelationFast(n_t, natoms, momenta, mass, VCF, vcf_intervall);
	
	double *omegas; double *powspec; int n_omegas = 1000; double domega = 0.001;
	omegas = malloc((n_omegas) * sizeof(double));
	powspec = malloc((n_omegas) * sizeof(double));
	
	powerspectrum(VCF, powspec, vcf_intervall, dt, omegas, n_omegas, domega);
	
	saveDataToFile("7/powerspectrum.csv", powspec, omegas, n_omegas, 1);
	
	free(momenta);
	free(VCF);
	free(powspec);
	free(omegas);
}

/* Main program */
int main()
{
	//runTask1();
	//runTask2(0.001);
	//runTask2(0.01);
	//runTask2(0.02);
	runTask3();
	runTask4();
	//runTask5('s');
	//runTask5('l');
	//runTask6('s');
	//runTask6('f');
	//runTask7();
	
	return 0;   
}
