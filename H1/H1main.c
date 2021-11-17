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

#define N_ATOMS 256

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
void saveEpotsToFile(char *fname, double *yvals, double *tvals, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "lattice const, Energy\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f, %f\n", tvals[i], yvals[i]);
    }
    fclose(fp);
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

/* Main program */
int main()
{
    /*
     * Task 1
	
    double *a0s; int N = 4; int natoms = N_ATOMS; double X[natoms][3];
    double *E_pots; int na = 1001; double da = 2/((double)na-1);

    a0s = malloc(na * sizeof(double));
    arange(a0s, 3.0, na, da);

    E_pots = malloc(na * sizeof(double));

    for (int i=0; i < na; i++) {
        init_fcc(X, N, a0s[i]);
        E_pots[i] = get_energy_AL(X, N * a0s[i] , natoms);
    }

    saveEpotsToFile("1/Epots_a0.csv", E_pots, a0s, na);

    free(E_pots);
    free(a0s);
    */



    /*
     * Task 2
     */
    
	double a0 = 4.03; int N = 4; int n_t = 1000; double dt = 2e-2; int natoms = N_ATOMS; double mass = 27.0 / 9649.0;
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
	saveEpotsToFile("2/Epot_dt0.02.csv", E_pot, time, n_t);
	saveEpotsToFile("2/Ekin_dt0.02.csv", E_kin, time, n_t);
	
	
    free(positions);
    free(momenta);
    free(E_pot);
    free(E_kin);

    /*
     Descriptions of the different functions in the files H1lattice.c and
     H1potential.c are listed below.
    */
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
    */
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */
    
    
    return 0;   
}
