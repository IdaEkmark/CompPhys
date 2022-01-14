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

#define REAL(z,i) ((z)[2*(i)]) 
#define IMAG(z,i) ((z)[2*(i)+1])

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


// Modulus squared of position space wave function
double psiAbsSq(double x, double d, double x0) {
	return 1/sqrt(PI * d*d) * exp(- (x-x0)*(x-x0)/(d*d));
}

// Modulus squared of momentum space wave function
double phiAbsSq(double p, double d, double p0) {
	return sqrt(d*d / (PI * HBAR*HBAR)) * exp(- d*d/(HBAR*HBAR) * (p-p0)*(p-p0));
}


gsl_complex psi(double x, double d, double x0, double p0) {
	return gsl_complex_mul_real(gsl_complex_exp( gsl_complex_rect(0.0, p0/HBAR * (x-x0)) ), pow(PI * d*d,-0.25) * exp(- (x-x0)*(x-x0)/(2 * d*d)));
}

void psiTest() {
	double kinE = 0.1; // p0^2/2m in eV
	double d = 0.5; // Å
	double p0 = sqrt(2*HMASS * kinE);

	// Calculated position space wavefunction (modulus squared)
	double xRange = 20.0; int nX = 401; double dx = xRange / (nX - 1.0);
	double *xVec; double *pDensPos;
	xVec = malloc(nX * sizeof(double));
	pDensPos = malloc(nX * sizeof(double));

	for (int i = 0; i < nX; i++) {
		xVec[i] = -xRange/2.0 + i * dx;
		pDensPos[i] = gsl_complex_abs2(psi(xVec[i], d, 0, p0));
	}

	saveDataToFile("1/pDensPosCalc.csv", xVec, pDensPos, nX, 1);

	free(xVec);
	free(pDensPos);
}


void runTask1() {
	double kinE = 0.1; // p0^2/2m in eV
	double d = 0.5; // Å

	// Theoretical position space wavefunction (modulus squared)
	double xRange = 20.0; int nX = 401; double dx = xRange / (nX - 1.0);
	double *xVec; double *pDensPos;
	xVec = malloc(nX * sizeof(double));
	pDensPos = malloc(nX * sizeof(double));

	for (int i = 0; i < nX; i++) {
		xVec[i] = -xRange/2.0 + i * dx;
		pDensPos[i] = psiAbsSq(xVec[i], d, 0);
	}

	saveDataToFile("1/pDensPos.csv", xVec, pDensPos, nX, 1);

	//free(xVec);
	free(pDensPos);


	// Theoretical momentum space wavefunction (modulus squared)
	double pRange = 20.0; int nP = 401; double dp = pRange / (nP - 1.0); double p0 = sqrt(2*HMASS * kinE);
	double *pVec; double *pDensMom;
	pVec = malloc(nP * sizeof(double));
	pDensMom = malloc(nP * sizeof(double));

	for (int i = 0; i < nP; i++) {
		pVec[i] = sqrt(2*HMASS * kinE) - pRange/2.0 + i * dp; // Assumes p0 is positive
		pDensMom[i] = phiAbsSq(pVec[i], d, p0);
	}

	saveDataToFile("1/pDensMomTheory.csv", pVec, pDensMom, nP, 1);

	free(pVec);
	free(pDensMom);


	// Calculated momentum space wavefunction (modulus squared)
	xRange = 40.0; nX = 401; dx = xRange / (nX - 1.0);
	double *psiData; gsl_complex psiHere; nP = nX; dp = 2*PI*HBAR/(dx*nX); //gsl_complex psiDataS[nX];
	psiData = malloc(2*nX * sizeof(double));
	pDensMom = malloc(nX * sizeof(double));
	pVec = malloc(nX * sizeof(double));

	for (int i=0; i<nX; i++) {
		xVec[i] = -xRange/2.0 + i * dx;
		psiHere = psi(xVec[i], d, 0, p0);
		REAL(psiData,i) = GSL_REAL(psiHere);
		IMAG(psiData,i) = GSL_IMAG(psiHere);

		pVec[i] = (i - (nX-1)/2) * dp;
	}

	/*Declare wavetable and workspace for fft*/
	gsl_fft_complex_wavetable *comp;
	gsl_fft_complex_workspace *work;
	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_complex_workspace_alloc(nX);
	comp = gsl_fft_complex_wavetable_alloc(nX);

	gsl_fft_complex_forward(psiData, 1, nX, comp, work);

	gsl_fft_complex_wavetable_free (comp);
  	gsl_fft_complex_workspace_free (work);

	for (int i=0; i<nX; i++) {
		if (i<(nX+1)/2) {
			pDensMom[i+(nX-1)/2] = REAL(psiData,i) * REAL(psiData,i) + IMAG(psiData,i) * IMAG(psiData,i);
			pDensMom[i+(nX-1)/2] *= dx*dx/(2*PI*HBAR);
		} else {
			pDensMom[i-(nX+1)/2] = REAL(psiData,i) * REAL(psiData,i) + IMAG(psiData,i) * IMAG(psiData,i);
			pDensMom[i-(nX+1)/2] *= dx*dx/(2*PI*HBAR);
		}
	}

	saveDataToFile("1/pDensMomCalc.csv", pVec, pDensMom, nP, 1);
	
	free(psiData);
	free(pDensMom);
	free(xVec);
	free(pVec);
}

void runTask2() {
	double kinE = 0.1; // p0^2/2m in eV
	double d = 0.5; // Å
	double p0 = sqrt(2*HMASS*kinE);

	int nX = 401; double xRange = 40.0; double dx = xRange / (nX - 1.0);
	int nT = 2; double tMax = 40.0; double dt = tMax / (nT - 1.0); double xp; double k; //double dp = 2*PI*HBAR/(dx*nX);
	double **psiList; double *xList; gsl_complex psiHere; //double *pList; double *tList;

	xList = malloc(nX * sizeof (double));
	psiList = malloc(nT * sizeof *psiList);
	for (int i=0; i < nT; i++) {
		psiList[i] = malloc(2*nX * sizeof *psiList);
	}

	// Evaluate psi for all relevant x at t=0
	for (int j=0; j < nX; j++) {
		xList[j] = -xRange/2.0 + j * dx;
		psiHere = psi(xList[j], d, 0, p0);
		REAL(psiList[0],j) = GSL_REAL(psiHere);
		IMAG(psiList[0],j) = GSL_IMAG(psiHere);

		//pList[j] = (j - (nX-1)/2) * dp;
	}

	/*Declare wavetable and workspace for fft*/
	gsl_fft_complex_wavetable *comp;
	gsl_fft_complex_workspace *work;
	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_complex_workspace_alloc(nX);
	comp = gsl_fft_complex_wavetable_alloc(nX);

	double *psiData; double newR; double newI;
	psiData = malloc(2*nX * sizeof(double));

	// Initialise psiData as Psi(x,t=0)
	for (int j=0; j < nX; j++) {
		REAL(psiData,j) = REAL(psiList[0],j);
		IMAG(psiData,j) = IMAG(psiList[0],j);
	}

	// Perform time evolution using FFTs
	for (int i=1; i < nT; i++) {
		// Potential is zero everywhere, and can be ignored
		// (Otherwise we would multiply by exp(-i V_j stuff) here)

		// Perform FFT into momentum space
		gsl_fft_complex_forward(psiData, 1, nX, comp, work);

		// Multiply each entry by exp(-iE(p)dt) to get evolution +dt
		for (int l=0; l < nX; l++) {
			k = (l-(nX-1)/2.0) * 2*PI/(nX*dx); // (l - (nX-1)/2) * dk
			xp = HBAR * k*k/(2.0*HMASS) * dt;
			newR = REAL(psiData,l)*cos(xp) - IMAG(psiData,l)*sin(xp);
			newI = IMAG(psiData,l)*cos(xp) + REAL(psiData,l)*sin(xp);
			REAL(psiData,l) = newR;
			IMAG(psiData,l) = newI;
		}

		// Perform IFFT back to position space
		gsl_fft_complex_inverse(psiData, 1, nX, comp, work);

		// Set current row of psiList to psiData
		for (int j=0; j < nX; j++) {
			REAL(psiList[i],j) = REAL(psiData,j);
			IMAG(psiList[i],j) = IMAG(psiData,j);
		}
	}

	double *pDensPos;
	pDensPos = malloc(nX * sizeof(double));
	for (int j=0; j < nX; j++) {
		pDensPos[j] = REAL(psiList[nT-1],j)*REAL(psiList[nT-1],j) + IMAG(psiList[nT-1],j)*IMAG(psiList[nT-1],j);
	}


	saveDataToFile("2/pDensPos40fs.csv", xList, pDensPos, nX, 1);

	free(pDensPos);

	gsl_fft_complex_wavetable_free (comp);
	gsl_fft_complex_workspace_free (work);
	free(psiData);

	for (int i=0; i < nT; i++) {
		free(psiList[i]);
	}
	free(psiList);
	free(xList);
}

int main() {
	//psiTest();
	//runTask1();
	runTask2();
	//printf("This complex number is %.2f + %.2f i\n", GSL_REAL(psi(1.0, 0.5, 0.0, sqrt(2*HMASS * 0.1))), GSL_IMAG(psi(1.0, 0.5, 0.0, sqrt(2*HMASS * 0.1))));
	
	return 0;
}