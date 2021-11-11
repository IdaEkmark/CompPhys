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
void calc_acc(double *a, double *u, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = (- 2*u[0] + u[1]);
    a[size_of_u - 1] = (u[size_of_u - 2] - 2*u[size_of_u - 1]);
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = (u[i - 1] - 2*u[i] + u[i + 1]);
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
void velocity_verlet(int n_timesteps, int n_particles, double **p, double **q, double dt)
{
    double u[n_particles];
    double v[n_particles];
    double a[n_particles];
    
    for (int i=0; i< n_particles; i++) {
    	u[i] = q[0][i];
    	v[i] = p[0][i];
    }
    calc_acc(a, u, n_particles);
    
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
        calc_acc(a, u, n_particles);
        
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
	    fprintf(fp, "%f", time_array[i]);

        for(int j = 0; j < n_particles; ++j){
            fprintf(fp, ", %f, %f, %f, %f", q[i][j], Q[i][j], p[i][j], P[i][j]);
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

int main()
{
	double dt = 0.1; double t_max = 25000; int n_t = 250000; int n_p = N_PARTICLES;
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
	
	velocity_verlet(n_t, n_p, p, q, dt);
	
	for (int t = 1; t < n_t; t++){
		transform_to_normal_modes(trans_matrix, n_p, q[t], Q[t]);
		transform_to_normal_modes(trans_matrix, n_p, p[t], P[t]);
	}
	
	double *time_array;
	time_array = malloc((n_t + 1) * sizeof (double));
	arange(time_array, 0, n_t, dt);
	saveqQpPtoFile("2/qQpP.csv", time_array, q, Q, p, P, n_t, n_p);
	
	for (int i=0; i < n_t+1; i++){
		free(q[i]);
	}
	free(q);
	for (int i=0; i < n_t+1; i++){
		free(Q[i]);
	}
	free(Q);
	for (int i=0; i < n_t+1; i++){
		free(p[i]);
	}
	free(p);
	for (int i=0; i < n_t+1; i++){
		free(P[i]);
	}
	free(P);
    return 0;
}
