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
void velocity_verlet(int n_timesteps, int n_particles, int n_cells, double ***pos, double ***mom, double dt, double alpha, double mass)
{
    double q[n_particles][3];
    double p[n_particles][3];
    double f[n_particles][3];

    int i;
    int j;
    
    for (i=0; i< n_particles; i++) {
        for (j=0; j<3; j++) {
            q[i][j] = pos[0][i][j];
            p[i][j] = mom[0][i][j];
        }
    }
    get_forces_AL(f,q,n_cells,n_particles);
    
    for (int t = 1; t < n_timesteps + 1; t++) {
        /* v(t+dt/2) */
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
        
        /* q(t+dt) */
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                q[i][j] += dt * p[i][j]/mass;
            }
        }
        
        /* a(t+dt) */
        get_forces_AL(f,q,n_cells,n_particles);
        
        /* v(t+dt) */
        for (i = 0; i < n_particles; i++) {
            for (j=0; j<3; j++) {
                p[i][j] += dt * 0.5 * f[i][j];
            }
        }
		
        /* Save the displacement of the three atoms */
        for (i = 0; i<n_particles; i++) {
            for (j = 0; j<3; j++) {
                pos[t][i][j] = q[i][j];
                mom[t][i][j] = p[i][j];
            }
        }
    }
}

/* Main program */
int main()
{
    /*
     * Task 1

    double *a0s; int N = 4; int natoms = 256; double X[natoms][3];
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
    
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
    srand(time(NULL));
    double random_value;
    random_value = (double) rand() / (double) RAND_MAX;
    */
	double a0 = 4.03; int N = 4; int n_t = 1000; int natoms = 256; double ***positions; double ***momenta;
	int t; int i; int j;
	
	positions = malloc((n_t + 1) * sizeof **positions);
	for (t = 0; t < n_t + 1; t++){
		positions[t] = malloc(natoms * sizeof **positions[t]);
		for (i = 0; i < natoms; i++) {
			positions[t][i] = malloc(3 * sizeof *positions[t][i]);
		}
	}
	momenta = malloc((n_t + 1) * sizeof **momenta);
	for (t = 0; t < n_t + 1; t++){
		momenta[t] = malloc(natoms * sizeof **momenta[t]);
		for (i = 0; i < natoms; i++) {
			momenta[t][i] = malloc(3 * sizeof *momenta[t][i]);
		}
	}
	
	init_fcc(positions[0], N, a0);
	
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	double randomUni; 
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			positions[0][i][j] += a0 * (-0.065 + 0.13 * gsl_rng_uniform(r));
		}
	}
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			momenta[0][i][j] = 0;
		}
	}
	 
	
	for (t = 0; t < n_t + 1; t ++){
		for (i=0; i < n_t+1; i++){
			free(positions[t][i]);
		}
		free(positions[t]);
	}
	free(positions);
	for (t = 0; t < n_t + 1; t ++){
		for (i=0; i < n_t+1; i++){
			free(momenta[t][i]);
		}
		free(momenta[t]);
	}
	free(momenta);
    
    /*
     Descriptions of the different functions in the files H1lattice.c and
     H1potential.c are listed below.
    */
    
    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    /*
     init_fcc(pos, Nc, a0);
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
