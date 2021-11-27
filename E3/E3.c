#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#define PI 3.141592653589

void saveDataToFile1D(char *fname, double *yvals, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Data\n");
    for(int i = 0; i < n_points; i++){
	    fprintf(fp, "%f,\n", yvals[i]);
    }
    fclose(fp);
}

void saveDataToFile(char *fname, double *yvals, double *tvals, int n_points, int n_skip)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "xData, yData\n");
    for(int i = 0; i < n_points; i += n_skip){
	    fprintf(fp, "%f, %f\n", tvals[i], yvals[i]);
    }
    fclose(fp);
}

void readData(char *fname, double *array, int N)
{
	int i;
	FILE *fp = fopen(fname, "r");
	for (i=0; i<N; i++) {
		if ( fscanf (fp," %lf", &(array[i])) == 0) {
			printf("no\n");
		}
	}
	fclose(fp);
}

double evalWeightFunction3(double x, double y, double z) {
	return pow(PI, -3.0/2.0) * exp( - x*x - y*y - z*z);
}

double evalFunction3(double x, double y, double z) {
	return x*x *(1 + y*y *(1 + z*z));// * evalWeightFunction3(x, y, z);
}

double monteCarloIntegrationTask1(double a, double b, int N, gsl_rng * r, double x_i_vec[], double f_i_vec[])  {
	double I = 0; double f_i; double x_i;
	for (int i = 0; i < N; i++) {
		x_i = a + (b-a)*gsl_rng_uniform(r);
		x_i_vec[i] = x_i;
		f_i = x_i * (1 - x_i);
		f_i_vec[i] = f_i;
		I += f_i;
	}
	I *= (b - a) / (double) N;
	return I;
}

double monteCarloIntegrationTask2(int N, gsl_rng * r, double x_i_vec[], double g_i_vec[])  {
	double I = 0; double g_i; double x_i; double xi_i;
	for (int i = 0; i < N; i++) {
		xi_i = 2*gsl_rng_uniform(r)/PI;
		x_i = acos(1 - PI * xi_i)/PI;
		x_i_vec[i] = x_i;
		g_i = x_i * (1 - x_i) / sin(PI * x_i) * 2 / PI;
		g_i_vec[i] = g_i;
		I += g_i;
	}
	I *= 1 / (double) N;
	return I;
}

void task12() {
	double a = 0; double b = 1; int N[4] = {10, 100, 1000, 10000}; double I[4];
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r, 3);
	double x_i_mat[4][N[3]]; double f_i_mat[4][N[3]]; double g_i_mat[4][N[3]];
	
	printf("Task 1\n");
	for (int i = 0; i < 4; i++){
		I[i] = monteCarloIntegrationTask1(a, b, N[i], r, x_i_mat[i], f_i_mat[i]);
		printf("I=%.4f for N=%i, should be 0.1667\n", I[i], N[i]);
	}
	saveDataToFile("1/f_x_i_10.csv", f_i_mat[0], x_i_mat[0], N[0], 1);
	saveDataToFile("1/f_x_i_100.csv", f_i_mat[1], x_i_mat[1], N[1], 1);
	saveDataToFile("1/f_x_i_1000.csv", f_i_mat[2], x_i_mat[2], N[2], 1);
	saveDataToFile("1/f_x_i_10000.csv", f_i_mat[3], x_i_mat[3], N[3], 1);
	
	printf("\n\nTask 2\n");
	for (int i = 0; i < 4; i++){
		I[i] = monteCarloIntegrationTask2(N[i], r, x_i_mat[i], g_i_mat[i]);
		printf("I=%.4f for N=%i, should be 0.1667\n", I[i], N[i]);
	}
	
	saveDataToFile("2/g_x_i_10.csv", g_i_mat[0], x_i_mat[0], N[0], 1);
	saveDataToFile("2/g_x_i_100.csv", g_i_mat[1], x_i_mat[1], N[1], 1);
	saveDataToFile("2/g_x_i_1000.csv", g_i_mat[2], x_i_mat[2], N[2], 1);
	saveDataToFile("2/g_x_i_10000.csv", g_i_mat[3], x_i_mat[3], N[3], 1);
}

double metropolisTask3(double delta, gsl_rng * r, int N, double fraction_use){
	double x_m = 0; double x_t = x_m + delta * (gsl_rng_uniform(r) - 0.5);
	double y_m = 0; double y_t = y_m + delta * (gsl_rng_uniform(r) - 0.5); 
	double z_m = 0; double z_t = z_m + delta * (gsl_rng_uniform(r) - 0.5);
	double p_m = evalWeightFunction3(x_m, y_m, z_m); double p_t;
	double I = 0; double normfactor = 1/(N * fraction_use);
	int counter = 0;
	for (int i = 0; i < N; i++) {
		x_t = x_m + delta * (gsl_rng_uniform(r) - 0.5);
		y_t = y_m + delta * (gsl_rng_uniform(r) - 0.5); 
		z_t = z_m + delta * (gsl_rng_uniform(r) - 0.5);
		
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
	printf("accaptancerate = %.4f\n", (double)counter / N);
	return I;
}

void task3() {
	printf("\n\nTask 3\n");
	double delta = 2; int N = 100000000; double fraction_use = 0.8;
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r, 17);
	
	double I = metropolisTask3(delta, r, N, fraction_use);

	printf("I = %.4e\n", I);
}

double evalCorrelationFunction1(double f[], int k, int N) {
	double f_mean2 = 0; double f2_mean = 0; double f_corr_mean = 0;
	int i;
	for (i = 0; i < N ; i++) {
		f_mean2 += f[i];
		f2_mean += f[i]*f[i];
	}
	f_mean2 = f_mean2 * f_mean2 /((double) N * N);
	f2_mean /= (double) N;
	for (i = 0; i < N-k; i++) {
		f_corr_mean += f[i]*f[i+k];
	}
	f_corr_mean /= (double) (N - k);
	
	return (f_corr_mean - f_mean2)/(f2_mean - f_mean2);
}

double evalCorrelationFunction2(double f[], int k, int N) {
	double f_mean = 0; double denom = 0; double num = 0;
	int i;
	for (i = 0; i < N ; i++) {
		f_mean += f[i];
	}
	f_mean /= (double) N;
	for (i = 0; i < N ; i++) {
		denom += (f[i] - f_mean)*(f[i] - f_mean);
	}
	for (i = 0; i < N-k; i++) {
		num += (f[i] - f_mean)*(f[i+k] - f_mean);
	}
	
	return num / denom;
}

int evalStatisticalInefficiency1(double f[], int N, double phi_vec[]) {
	int k = 0; double phi = 100;
	while (phi > 0.135) {
		k++; 
		phi = evalCorrelationFunction1(f, k, N);
		phi_vec[k-1] = phi;
	}
	return k;
}

int evalStatisticalInefficiency2(double f[], int N) {
	int k = 0; double phi = 100;
	while (phi > 0.135) {
		phi = evalCorrelationFunction2(f, k, N);
		k++;
	}
	
	return k;
}

void evalBlockVariables(double f[], double F[], int N, int N_B, int B) {
	for (int j = 0; j < N_B; j++) {
		F[j] = 0;
		for (int i = 0; i < B; i++) {
			F[j] += f[i + j*B];
		}
		F[j] /= (double) B;
	}
}

double evalVariance(double a[], int n) {
	double a_mean2=0; double a2_mean=0;
	for (int i=0; i<n; i++) {
		a_mean2 += a[i];
		a2_mean += a[i]*a[i];
	}
	a_mean2 /= (double) n;
	a_mean2 = a_mean2 * a_mean2;
	a2_mean /= (double) n;
	return a2_mean - a_mean2;
}

double evalStatisticalInefficiencyBlockAtB(double f[], int N, int B) {
	int N_B = N / B; double *F; 
	F = malloc(N_B * sizeof(double));
	evalBlockVariables(f, F, N, N_B, B);
	double varf = evalVariance(f, N);
	double varF = evalVariance(F, N_B);
	double s = B * varF / varf;
	return s;
}

void evalStatisticalInefficiencyBlock(double f[], int N, double s_vec[], double B_vec[]) {
	double s; int B = 1; int i = 0;
	while (B <= N/2) {
		if ( N % B == 0 ) {
			s = evalStatisticalInefficiencyBlockAtB(f, N, B);
			s_vec[i] = s;
			B_vec[i] = B;
			i++;
			//printf("Kan dela en miljon med %i heltal, B=%i, s=%.1f\n", i, B, s);
		}
		/*if (i % 100 == 0) {
			printf("i=%i\n", i);
		}*/
		B += 1;
		
	}
}

void evalMeanedA(double a[], double a_meaned[], int n){
	double a_mean=0;
	for (int i=0; i<n; i++) {
		a_mean += a[i];
	}
	a_mean /= (double) n;
	for (int i=0; i<n; i++){
		a_meaned[i] = a[i] - a_mean;
	}
}

void task4() {
	int N = 1000000; double *f; double *f_meaned; double s_vec[48]; double B_vec[48]; double phi_vec[181];
	f = malloc((N) * sizeof(double));
	f_meaned = malloc((N) * sizeof(double));
	readData("MC.txt", f, N);
	
	evalMeanedA(f, f_meaned, N);

	printf("\n\nTask 4\n");
	int s = evalStatisticalInefficiency1(f_meaned, N, phi_vec);
	printf("s = %i\n", s);
	saveDataToFile1D("4/phi.csv", phi_vec, 181);
	//s = evalStatisticalInefficiency2(f, N);
	//printf("s_2 = %i\n", s);
	
	//evalStatisticalInefficiencyBlock(f_meaned, N, s_vec, B_vec);
	//saveDataToFile("4/block.csv", s_vec, B_vec, 48, 1);
	
	free(f);
	free(f_meaned);
}


int main() {
	//task12();
	//task3();
	task4();
	
	return 0;
}