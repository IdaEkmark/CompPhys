#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>

#define PI 3.141592653589
#define HBAR 6.582119569e-1 // eV * fs
#define HMASS 104.453702 // Mass of hydrogen atom, eV / (Å/fs)^2

void saveDataToFile(char *fname, double *xvals, double *yvals, int n_points, int n_skip)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "xData, yData\n");
    for(int i = 0; i < n_points; i = i + n_skip){
	    fprintf(fp, "%f, %f\n", xvals[i], yvals[i]);
    }
    fclose(fp);
}

void saveDataToFile1D(char *fname, double *yvals, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Data\n");
    for(int i = 0; i < n_points; i++){
	    fprintf(fp, "%f,\n", yvals[i]);
    }
    fclose(fp);
}

void readDataFromFile(char *fname, double *xvals, double *yvals)
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
    int i = 0;
    while(fgets(line, sizeof(line), fp) != NULL){
	token = strtok(line, ",");
	xvals[i] = strtod(token, NULL);
	token = strtok(NULL, ",");
	yvals[i] = strtod(token, NULL);
	i++;
	memset(line, 0, sizeof(line));
	token = NULL;
    }
    fclose(fp);
}


double psiAbsSq(double x, double d, double x0) {
	return 1/sqrt(PI * d*d) * exp(- (x-x0)*(x-x0)/(d*d));
}


void runTask1() {
	double kinE = 0.1; // p0^2/2m in eV
	double d = 0.5; // Å
	double xMax = 5.0; int nX = 201; double dx = 2*xMax / (nX - 1.0);
	double *xVec; double *pDens;
	xVec = malloc(nX * sizeof(double));
	pDens = malloc(nX * sizeof(double));

	for (int i = 0; i < nX; i++) {
		xVec[i] = -xMax + i * dx;
		pDens[i] = psiAbsSq(xVec[i], d, 0);
	}

	saveDataToFile("1/pDensPos.csv", xVec, pDens, nX, 1);

	free(xVec);
	free(pDens);
}

int main() {
	runTask1();
	//runTask2();
	
	return 0;
}