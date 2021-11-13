/* 
 * Help routines for E2
 *
 * E2code.c
 * 
 * Compile me as:
 * clang E2code.c -o ./Executable_files/<executable name> -lm
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N_PARTICLES 32
#define PI 3.141592653589

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
void calc_acc(double *a, double *u, int size_of_u, double alpha)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = (- 2*u[0] + u[1]) 
    		+ alpha * ( u[1]*u[1] - 2*u[1]*u[0] );
    a[size_of_u - 1] = (u[size_of_u - 2] - 2*u[size_of_u - 1]) 
    		+ alpha * ( 2*u[size_of_u - 1]*u[size_of_u - 2] - u[size_of_u - 2]*u[size_of_u - 2] );
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = (u[i - 1] - 2*u[i] + u[i + 1]) 
        		+ alpha * ( u[i+1]*u[i+1] - u[i-1]*u[i-1] + 2*u[i-1]*u[i] - 2*u[i]*u[i+1] );
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
void velocity_verlet(int n_timesteps, int n_particles, double **p, double **q, double dt, double alpha)
{
    double u[n_particles];
    double v[n_particles];
    double a[n_particles];
    
    for (int i=0; i< n_particles; i++) {
    	u[i] = q[0][i];
    	v[i] = p[0][i];
    }
    calc_acc(a, u, n_particles, alpha);
    
    for (int t = 1; t < n_timesteps + 1; t++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            u[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, u, n_particles, alpha);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
		
        /* Save the displacement of the three atoms */
        for (int i = 0; i<n_particles; i++) {
        	q[t][i] = u[i];
        	p[t][i] = v[i];
        }
    }
}
/*
 * trans_matrix[N_PARTICLES][N_PARTICLES]: empty allocated array which
 * will be filled with sine transformation matrix
 * N_PARTICLES: number of particles in system
 */
void construct_transformation_matrix(
    double trans_matrix[N_PARTICLES][N_PARTICLES], int n_particles)
{
    double factor = 1 / ((double)n_particles + 1);
    for(int i = 0; i < n_particles; i++){
	for(int j = 0; j < n_particles; j++){
	    trans_matrix[i][j] = sqrt(2 * factor)
				 * sin((j + 1) * (i + 1) * PI * factor);
	}
    }
}

/*
 * Transformation matrix constucted in above function
 * q cartesian coordinate of paricles
 * Q output normal modes coordinate
 * N_PARTICLES is number of particles in system
 */
void transform_to_normal_modes(double trans_matrix[N_PARTICLES][N_PARTICLES],
			       int n_particles,
			       double *q, double *Q)
{
    for(int i = 0; i < n_particles; i++){
	double sum = 0;
	for(int j = 0; j < n_particles; j++){
	    sum += q[j] * trans_matrix[i][j];
	}
	Q[i] = sum;
    }
}

void saveqQpPtoFile(char *fname, double *time_array,
		   double **q, double **Q, double **p, double **P, int n_timesteps, int n_particles) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "time");
    for (int j = 1; j < n_particles+1; ++j) {
        fprintf(fp, ", q%i, Q%i, p%i, P%i", j, j, j, j);
    }
    fprintf(fp, "\n");

    for(int i = 0; i < n_timesteps; ++i){
	    fprintf(fp, "%e", time_array[i]);

        for(int j = 0; j < n_particles; ++j){
            fprintf(fp, ", %e, %e, %e, %e", q[i][j], Q[i][j], p[i][j], P[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
    }
}


void calculateEnergy(int n_timesteps, int n_particles, double **Q, double **P,
        double **kinE, double **potE, double **totE) 
{
    double wSq;
	for (int t = 0; t < n_timesteps + 1; t++) { 
		for (int p = 0; p < n_particles; p++) {
            wSq = 4 * sin( (p+1)*PI /  ((double) 2*(n_particles + 1)) )*sin( (p+1)*PI / ((double) 2*(n_particles + 1)) );
			kinE[t][p] = 0.5 * P[t][p]*P[t][p];
            potE[t][p] = 0.5 * wSq * Q[t][p]*Q[t][p];
            totE[t][p] = kinE[t][p] + potE[t][p];
		}
	}
}

/*
 * Saves t, kinE, potE and totE to a file (Usually .csv).
 * @fname - File name 
 * @time_array - array of time values
 * @Q - n_particles x n_timesteps array with displacement values
 * @V - n_particles x n_timesteps array with velocity values
 * @n_timesteps - number of timesteps
 * @n_particles - number of particles
*/
void saveKPEtoFile(char *fname, double *time_array,
		   double **kinE, double **potE, double **totE, int n_timesteps, int n_particles) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "time");
    for (int j = 1; j < n_particles + 1; ++j) {
        fprintf(fp, ", kinE%i, potE%i, totE%i", j, j, j);
    }
    fprintf(fp, "\n");

    for(int i = 0; i < n_timesteps + 1; ++i){
	    fprintf(fp, "%e", time_array[i]);

        for(int j = 0; j < n_particles; ++j){
            fprintf(fp, ", %e, %e, %e", kinE[i][j], potE[i][j], totE[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int main()
{
	double t_max = 25000.0; int n_t = 250000; double dt = t_max / (double) n_t; int n_p = N_PARTICLES; double alpha = 0.01;
	double trans_matrix[n_p][n_p];
    double **q;
    double **Q;
    double **p;
	double **P;
	
	construct_transformation_matrix(trans_matrix, n_p);
	
	q = malloc((n_t + 1) * sizeof *q);
	for (int i=0; i < n_t + 1; i++){
		q[i] = malloc(n_p * sizeof *q[i]);
	}
	Q = malloc((n_t + 1) * sizeof *Q);
	for (int i=0; i < n_t + 1; i++){
		Q[i] = malloc(n_p * sizeof *Q[i]);
	}
	p = malloc((n_t + 1) * sizeof *q);
	for (int i=0; i < n_t + 1; i++){
		p[i] = malloc(n_p * sizeof *q[i]);
	}
	P = malloc((n_t + 1) * sizeof *Q);
	for (int i=0; i < n_t + 1; i++){
		P[i] = malloc(n_p * sizeof *Q[i]);
	}
	
	Q[0][0] = 0.0;
	P[0][0] = 8.0;
	for (int i = 1; i < n_p; i++){
		Q[0][i] = 0.0;
		P[0][i] = 0.0;
	}
	
	transform_to_normal_modes(trans_matrix, n_p, Q[0], q[0]);
	transform_to_normal_modes(trans_matrix, n_p, P[0], p[0]);
	
	/*
	for (int pi = 0; pi < n_p; pi++) {
		printf("q_pi: %.2e\n", q[0][pi]);
	}
	for (int pi = 0; pi < n_p; pi++) {
			printf("p_pi: %.2e\n", p[0][pi]);
	}
	*/
	
	velocity_verlet(n_t, n_p, p, q, dt, alpha);
	
	for (int t = 1; t < n_t + 1; t++){
		transform_to_normal_modes(trans_matrix, n_p, q[t], Q[t]);
		transform_to_normal_modes(trans_matrix, n_p, p[t], P[t]);
	}
	
	double *time_array;
	time_array = malloc((n_t + 1) * sizeof (double));
	arange(time_array, 0, n_t+1, dt);
	saveqQpPtoFile("3/qQpP_0.01.csv", time_array, q, Q, p, P, n_t, n_p);
	
	for (int i=0; i < n_t+1; i++){
		free(q[i]);
	}
	free(q);
	for (int i=0; i < n_t+1; i++){
		free(p[i]);
	}
	free(p);
    
    // Energies per mode at each timestep
    double **kinE;
    double **potE;
    double **totE;

    kinE = malloc((n_t + 1) * sizeof *kinE);
	for (int i=0; i < n_t + 1; i++){
		kinE[i] = malloc(n_p * sizeof *kinE[i]);
	}
	potE = malloc((n_t + 1) * sizeof *potE);
	for (int i=0; i < n_t + 1; i++){
		potE[i] = malloc(n_p * sizeof *potE[i]);
	}
	totE = malloc((n_t + 1) * sizeof *totE);
	for (int i=0; i < n_t + 1; i++){
		totE[i] = malloc(n_p * sizeof *totE[i]);
	}

    calculateEnergy(n_t, N_PARTICLES, Q, P, kinE, potE, totE);
    saveKPEtoFile("3/KPE_0.01.csv", time_array, kinE, potE, totE, n_t, N_PARTICLES);

    free(time_array);

    for (int i=0; i < n_t+1; i++){
		free(kinE[i]);
	}
	free(kinE);
	for (int i=0; i < n_t+1; i++){
		free(potE[i]);
	}
	free(potE);
    for (int i=0; i < n_t+1; i++){
		free(totE[i]);
	}
	free(totE);

    for (int i=0; i < n_t+1; i++){
		free(Q[i]);
	}
	free(Q);
	for (int i=0; i < n_t+1; i++){
		free(P[i]);
	}
	free(P);
    return 0;
}
