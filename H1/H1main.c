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

/* Main program */
int main()
{
    double **X; double *a0s; int N = 4; int natoms = 256; int ndim = 3;
    double *E_pots; double na = 101; double da = 0.02;

    a0s = malloc(na * sizeof(double));
    arange(a0s, 3.0, na, da);

    E_pots = malloc(na * sizeof(double));

    X = malloc(natoms * sizeof *X);
	for (int i=0; i < natoms; i++){
		X[i] = malloc(ndim * sizeof *X[i]);
	}

    for (int i=0; i < na; i++) {
        init_fcc(X, N, a0s[i]);
        E_pots[i] = get_energy_AL(X, N * a0s[i] , natoms);
    }

    saveEpotsToFile("Epots_a0.csv", E_pots, a0s, na);

    free(E_pots);
    free(a0s);

    for (int i=0; i < natoms; i++){
		free(X[i]);
	}
	free(X);
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
     srand(time(NULL));
     double random_value;
     random_value = (double) rand() / (double) RAND_MAX;
    */
    
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
