#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>

#define PI 3.141592653589
#define K_B 8.617342e-5

void saveDataToFile(char *fname, double *yvals, double *xvals, int n_points, int n_skip)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "xData, yData\n");
    for(int i = 0; i < n_points; i = i + n_skip){
	    fprintf(fp, "%f, %f\n", xvals[i], yvals[i]);
    }
    fclose(fp);
}

void saveDataBlockToFile(char *fname, double (*yvals)[6], double *xvals, int n_points, int n_skip) {
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "xData, y10000, y16000, y20000, y25000, y32000, y40000\n");

    for(int t = 0; t < n_points; t = t + n_skip){
	    fprintf(fp, "%f", xvals[t]);
	    
	    for (int i = 0; i < 6; i++) {
			fprintf(fp, ", %f", yvals[t][i]);
	    }
        fprintf(fp, "\n");
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

double evalE0(double E_AA, double E_BB, double E_AB, int N) {
	return 2*N * (E_AA + E_BB + 2*E_AB);
}

double evalDeltaE(double E_AA, double E_BB, double E_AB, int N) {
	return E_AA + E_BB - 2*E_AB;
}

double evalU(double P, double E_AA, double E_BB, double E_AB, int N){
	double E0 = evalE0(E_AA, E_BB, E_AB, N); double deltaE = evalDeltaE(E_AA, E_BB, E_AB, N);
	return E0 - 2 * N * P*P * deltaE;
}
double evalUderiv(double P, double E_AA, double E_BB, double E_AB, int N){
	double deltaE = evalDeltaE(E_AA, E_BB, E_AB, N);
	return -4 * N * P * deltaE;
}

double evalS(double P, int N) {
	double term1 = 2 * N * K_B * log(2);
	double term2 = N * K_B * ((1+P)*log(1+P) + (1-P)*log(1-P));
	return term1 - term2;
}

double evalSderiv(double P, int N) {
	double term1 = 0;
	double term2 = N * K_B * (log(1+P) - log(1-P));
	return term1 - term2;
}

double evalFderiv(double T, double P, double E_AA, double E_BB, double E_AB, int N){
	double dU = evalUderiv(P, E_AA, E_BB, E_AB, N); double dS = evalSderiv(P, N);
	return dU - T * dS;
}


double minFreeEnergBisection(double T, double E_AA, double E_BB, double E_AB, int N, double tol){
	double P_min = -1; double P_max = 1; double P_crit;
	double dF_crit ;
	while (fabs(P_min - P_max) > tol) {
		P_crit = (P_min + P_max)/2;
		dF_crit = evalFderiv(T, P_crit, E_AA, E_BB, E_AB, N);
		if (dF_crit > 0) {
			P_max = P_crit;
		} else {
			P_min = P_crit;
		}
	}
	return (P_min + P_max)/2;
}

void runTask1() {
	double E_CuCu = -436e-3; double E_ZnZn = -113e-3; double E_CuZn = -294e-3;
	
	int N = 1000; double tol = 1e-6;
	double T = 200; double deltaT = 0.5; int nT = 1601; 
	double *T_vec; double *P_vec; double *U_vec;
	
	T_vec = malloc( nT * sizeof(double));
	P_vec = malloc( nT * sizeof(double));
	U_vec = malloc( nT * sizeof(double));
	for (int t=0; t<nT; t++){
		T_vec[t] = T;
		P_vec[t] = minFreeEnergBisection(T, E_CuCu, E_ZnZn, E_CuZn, N, tol);
		U_vec[t] = evalU(P_vec[t], E_CuCu, E_ZnZn, E_CuZn, N);
		// printf("Critical temperature: %.2f\n", 2.0/K_B * evalDeltaE(E_CuCu, E_ZnZn, E_CuZn, N) - 273.15);
		T += deltaT;
	}
	
	saveDataToFile("1/TP.csv", P_vec, T_vec, nT, 1);
	saveDataToFile("1/TU.csv", U_vec, T_vec, nT, 1);
		
	free(T_vec);
	free(P_vec);
	free(U_vec);
}

void initializeNeighbourMatrix(int nUnitCellLengths, int (*neighbourMatrix)[8], int nAtoms){
	int sc1_mat3d[nUnitCellLengths+1][nUnitCellLengths+1][nUnitCellLengths+1]; 
	int sc2_mat3d[nUnitCellLengths+1][nUnitCellLengths+1][nUnitCellLengths+1];
	
	/*
	 * Building first of two sc-lattices with nUnitCellLengths^3 unit cells
	 * Giving vertex in the lattice a index, periodicity included 
	 * Putting the indexes in either sc-lattice into the corresponding sc-list
	 */
	int i; int j; int k; int index=0; 
	for (i=0; i<nUnitCellLengths; i++) {
		for (j=0; j<nUnitCellLengths; j++) {
			for (k=0; k<nUnitCellLengths; k++) {
				sc1_mat3d[i][j][k] = index;
				index++;
			}
		}
	}
	for (j=0; j<nUnitCellLengths; j++) {
		for (k=0; k<nUnitCellLengths; k++) {
			sc1_mat3d[nUnitCellLengths][j][k] = sc1_mat3d[0][j][k];
		}
	}
	for (i=0; i<nUnitCellLengths; i++) {
		for (k=0; k<nUnitCellLengths; k++) {
			sc1_mat3d[i][nUnitCellLengths][k] = sc1_mat3d[i][0][k];
		}
	}
	for (i=0; i<nUnitCellLengths; i++) {
		for (j=0; j<nUnitCellLengths; j++) {
			sc1_mat3d[i][j][nUnitCellLengths] = sc1_mat3d[i][j][0];
		}
	}
	for (k=0; k<nUnitCellLengths; k++) {
		sc1_mat3d[nUnitCellLengths][nUnitCellLengths][k] = sc1_mat3d[0][0][k];
	}
	for (j=0; j<nUnitCellLengths; j++) {
		sc1_mat3d[nUnitCellLengths][j][nUnitCellLengths] = sc1_mat3d[0][j][0];
	}
	for (i=0; i<nUnitCellLengths; i++) {
		sc1_mat3d[i][nUnitCellLengths][nUnitCellLengths] = sc1_mat3d[i][0][0];
	}
	sc1_mat3d[nUnitCellLengths][nUnitCellLengths][nUnitCellLengths] = sc1_mat3d[0][0][0];
	
	/*
	 * Building second of two sc-lattices with nUnitCellLengths^3 unit cells
	 * Giving vertex in the lattice a index, periodicity included 
	 * Putting the indexes in either sc-lattice into the corresponding sc-list
	 */
	for (i=1; i<nUnitCellLengths+1; i++) {
		for (j=1; j<nUnitCellLengths+1; j++) {
			for (k=1; k<nUnitCellLengths+1; k++) {
				sc2_mat3d[i][j][k] = index;
				index++;
			}
		}
	}
	for (j=1; j<nUnitCellLengths+1; j++) {
		for (k=1; k<nUnitCellLengths+1; k++) {
			sc2_mat3d[0][j][k] = sc2_mat3d[nUnitCellLengths][j][k];
		}
	}
	for (i=1; i<nUnitCellLengths+1; i++) {
		for (k=1; k<nUnitCellLengths+1; k++) {
			sc2_mat3d[i][0][k] = sc2_mat3d[i][nUnitCellLengths][k];
		}
	}
	for (i=1; i<nUnitCellLengths+1; i++) {
		for (j=1; j<nUnitCellLengths+1; j++) {
			sc2_mat3d[i][j][0] = sc2_mat3d[i][j][nUnitCellLengths];
		}
	}
	for (k=1; k<nUnitCellLengths+1; k++) {
		sc2_mat3d[0][0][k] = sc2_mat3d[nUnitCellLengths][nUnitCellLengths][k];
	}
	for (j=1; j<nUnitCellLengths+1; j++) {
		sc2_mat3d[0][j][0] = sc2_mat3d[nUnitCellLengths][j][nUnitCellLengths];
	}
	for (i=1; i<nUnitCellLengths+1; i++) {
		sc2_mat3d[i][0][0] = sc2_mat3d[i][nUnitCellLengths][nUnitCellLengths];
	}
	sc2_mat3d[0][0][0] = sc2_mat3d[nUnitCellLengths][nUnitCellLengths][nUnitCellLengths];
	
	/*
	 * For each index in the sc-lattices, checking putting the neighbouring indices into the neighbour matrix
	 */
	int indexMatris = 0; 
	for (i=0; i<nUnitCellLengths; i++) {
		for (j=0; j<nUnitCellLengths; j++) {
			for (k=0; k<nUnitCellLengths; k++) {
				neighbourMatrix[indexMatris][0] = sc2_mat3d[i][j][k];
				neighbourMatrix[indexMatris][1] = sc2_mat3d[i+1][j][k];
				neighbourMatrix[indexMatris][2] = sc2_mat3d[i][j+1][k];
				neighbourMatrix[indexMatris][3] = sc2_mat3d[i][j][k+1];
				neighbourMatrix[indexMatris][4] = sc2_mat3d[i+1][j+1][k];
				neighbourMatrix[indexMatris][5] = sc2_mat3d[i+1][j][k+1];
				neighbourMatrix[indexMatris][6] = sc2_mat3d[i][j+1][k+1];
				neighbourMatrix[indexMatris][7] = sc2_mat3d[i+1][j+1][k+1];
				indexMatris++;
			}
		}
	}
	for (i=1; i<nUnitCellLengths+1; i++) {
		for (j=1; j<nUnitCellLengths+1; j++) {
			for (k=1; k<nUnitCellLengths+1; k++) {
				neighbourMatrix[indexMatris][0] = sc1_mat3d[i][j][k];
				neighbourMatrix[indexMatris][1] = sc1_mat3d[i-1][j][k];
				neighbourMatrix[indexMatris][2] = sc1_mat3d[i][j-1][k];
				neighbourMatrix[indexMatris][3] = sc1_mat3d[i][j][k-1];
				neighbourMatrix[indexMatris][4] = sc1_mat3d[i-1][j-1][k];
				neighbourMatrix[indexMatris][5] = sc1_mat3d[i-1][j][k-1];
				neighbourMatrix[indexMatris][6] = sc1_mat3d[i][j-1][k-1];
				neighbourMatrix[indexMatris][7] = sc1_mat3d[i-1][j-1][k-1];
				indexMatris++;
			}
		}
	}
}

void newConfiguration(int atomList[], int *newAtomList, int nAtoms, gsl_rng * r){
	int atom1 = floor(nAtoms * gsl_rng_uniform(r)); int atom2 = floor(nAtoms * gsl_rng_uniform(r));
	
	while(atomList[atom1] == atomList[atom2]){
		atom2 = floor(nAtoms * gsl_rng_uniform(r));
	}
	
	for(int i=0; i<nAtoms; i++) {
		if (i == atom1) {
			newAtomList[i] = atomList[atom2];
		} else if (i == atom2) {
			newAtomList[i] = atomList[atom1];
		} else {
			newAtomList[i] = atomList[i];
		}
	}
}

double evalEnergyForState(int (*neighbourMatrix)[8], int atomList[], int nAtoms, double E_AA, double E_BB, double E_AB){
	double E = 0; int i; int j;
	for(i=0; i<nAtoms; i++) {
		if (atomList[i] == 0) {
			for(j=0; j<8; j++) {
				if(atomList[neighbourMatrix[i][j]] == 0) {
					E += E_AA/2.0;
				} else {
					E += E_AB/2.0;
				}
			}
		} else {
			for(j=0; j<8; j++) {
				if(atomList[neighbourMatrix[i][j]] == 1) {
					E += E_BB/2.0;
				} else {
					E += E_AB/2.0;
				}
			}
		}
	}
	return E;
}

double evalLongRangeOrder(int atomList[], int nAtoms) {
	int N_A = 0;
	for (int i=0; i<nAtoms/2; i++){
		if (atomList[i] == 0) {
			N_A++;
		}
	}
	return ((double) N_A /(double) nAtoms) * 4.0 - 1.0; 
}

double evalShortRangeOrder(int (*neighbourMatrix)[8], int atomList[], int nAtoms){
	double q = 0; int i; int j; double N = (double) nAtoms/2.0;
	for(i=0; i<nAtoms; i++) {
		for(j=0; j<8; j++) {
			if(atomList[i] != atomList[neighbourMatrix[i][j]]) {
				q += 1.0/2.0;
			}
		}
	}
	return 1/(4.0*N)*(q - 4.0*N);
}

void evalMeanedA(double *a, double *a_meaned, int n){
	double a_mean=0;
	for (int i=0; i<n; i++) {
		a_mean += a[i];
	}
	a_mean /= (double) n;
	for (int i=0; i<n; i++){
		a_meaned[i] = a[i] - a_mean;
	}
}

double evalCorrelationFunction(double f[], int k, int N, double f_mean2, double f2_mean) {
	double f_corr_mean = 0; 
	int i;
	for (i = 0; i < N; i++) {
		f_corr_mean += f[i]*f[i+k];
	}
	f_corr_mean /= (double) (N - k);
	double num = (f_corr_mean - f_mean2);
	double denom = (f2_mean - f_mean2);
	
	/*if(fabs(num) < tol && fabs(denom) < tol) {
		return -1.0;
	} else if(fabs(denom) < tol) {
		return 1/tol;
	} else {*/
	return num/denom; 
	//}
}

int evalStatisticalInefficiency(double f[], int N) {
	double f_mean2 = 0; double f2_mean = 0; double tol = 1.0e-25;
	int i;
	for (i = 0; i < N ; i++) {
		f_mean2 += f[i];
		f2_mean += f[i]*f[i];
	}
	f_mean2 = f_mean2 * f_mean2 /((double) N * N);
	f2_mean /= (double) N;
	if (fabs(f2_mean-f_mean2) < tol) {
		return -1;
	} 
	
	double phi = 1; int k = 0; 
	while (phi > 0.135 && k < N/2) { 
		if(k<10000) {
			k++;
		} else if(k<100000) {
			k += 10;
		} else {
			k += 100;
		}
		phi = evalCorrelationFunction(f, k, N/2, f_mean2, f2_mean);
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
	int N_B = N / B; double *F; double tol = 1e-25;
	F = malloc(N_B * sizeof(double));
	evalBlockVariables(f, F, N, N_B, B);
	double varf = evalVariance(f, N);
	double varF = evalVariance(F, N_B);
	if( varf < tol) {
		varF = 1.0;
		varf = -B;
	}
	double s = B * varF / varf;
	return s;
}

void metropolis(double T, int nAtoms, int nUnitCellLengths, int (*neighbourMatrix)[8], int N_tot, int N_eq, 
		double E_AA, double E_BB, double E_AB, gsl_rng * rng, double *UCPr, double *se_corr, double *se_block, double *var){
	
	int atomList_t[nAtoms]; int atomList_m[nAtoms];
	for (int i=1; i<nAtoms; i++) {
		if (i<nAtoms/2) {
			atomList_m[i] = 0;
		} else {
			atomList_m[i] = 1;
		}
	}	
	double E_m = evalEnergyForState(neighbourMatrix, atomList_m, nAtoms, E_AA, E_BB, E_AB); 
	double E_t;
	
	double beta = 1/(K_B * T);
	double meanE=0; double meanE2=0; double meanP=0; double meanr=0;
	
	double *E_m_vec; double *P_vec; double *r_vec; double *E_m_vec_meaned; double *P_vec_meaned; double *r_vec_meaned;
	double P=1; double rr=1;
	
	//double s_E_m_corr; double s_P_corr; double s_r_corr; double s_E_m_block; double s_P_block; double s_r_block; 
	
	E_m_vec = malloc((N_tot - N_eq) * sizeof(double));
	P_vec = malloc((N_tot - N_eq) * sizeof(double));
	r_vec = malloc((N_tot - N_eq) * sizeof(double));
	E_m_vec_meaned = malloc((N_tot - N_eq) * sizeof(double));
	P_vec_meaned = malloc((N_tot - N_eq) * sizeof(double));
	r_vec_meaned = malloc((N_tot - N_eq) * sizeof(double));
	
	int t; int i;
	for (t=0; t<N_tot; t++) {
		/*if (t%100000==0) {
			printf("%i\n", t);
		}*/
		newConfiguration(atomList_m, atomList_t, nAtoms, rng);
		
		E_t = evalEnergyForState(neighbourMatrix, atomList_t, nAtoms, E_AA, E_BB, E_AB); 
		
		if ( gsl_rng_uniform(rng) < exp(-(E_t - E_m) * beta) ) {
			for (i=0; i<nAtoms; i++) {
				atomList_m[i] = atomList_t[i]; 
			}
			E_m = E_t;
		}
		if ( t >= N_eq ) {
			meanE += E_m;
			meanE2 += E_m * E_m;
			E_m_vec[t-N_eq] = E_m;
			P = fabs(evalLongRangeOrder(atomList_m, nAtoms));
			P_vec[t-N_eq] = P;
			rr = fabs(evalShortRangeOrder(neighbourMatrix, atomList_m, nAtoms));
			r_vec[t-N_eq] = rr;
			meanP += P;
			meanr += rr;
		}
	}
	meanE /= ((double) N_tot - N_eq);
	meanE2 /= ((double) N_tot - N_eq); 
	meanP  /= ((double) N_tot - N_eq);
	meanr  /= ((double) N_tot - N_eq);
	
	UCPr[0] = meanE;
	UCPr[1] = (meanE2 - meanE*meanE) * (beta*beta) * K_B;
	UCPr[2] = P;
	UCPr[3] = rr;
	UCPr[4] = meanP;
	UCPr[5] = meanr;
	
	var[0] = evalVariance(E_m_vec, N_tot - N_eq);
	var[1] = evalVariance(P_vec, N_tot - N_eq);
	var[2] = evalVariance(r_vec, N_tot - N_eq);
	
	evalMeanedA(E_m_vec, E_m_vec_meaned, N_tot - N_eq);
	evalMeanedA(P_vec, P_vec_meaned, N_tot - N_eq);
	evalMeanedA(r_vec, r_vec_meaned, N_tot - N_eq);
	
	se_corr[0] = evalStatisticalInefficiency(E_m_vec_meaned, N_tot - N_eq);
	se_corr[1] = evalStatisticalInefficiency(P_vec_meaned, N_tot - N_eq);
	se_corr[2] = evalStatisticalInefficiency(r_vec_meaned, N_tot - N_eq);
	
	se_block[0] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 10000);
	se_block[1] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 16000);
	se_block[2] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 20000);
	se_block[3] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 25000);
	se_block[4] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 32000);
	se_block[5] = evalStatisticalInefficiencyBlockAtB(E_m_vec_meaned, N_tot - N_eq, 40000);
	
	se_block[6] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 10000);
	se_block[7] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 16000);
	se_block[8] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 20000);
	se_block[9] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 25000);
	se_block[10] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 32000);
	se_block[11] = evalStatisticalInefficiencyBlockAtB(P_vec_meaned, N_tot - N_eq, 40000);

	se_block[12] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 10000);
	se_block[13] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 16000);
	se_block[14] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 20000);
	se_block[15] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 25000);
	se_block[16] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 32000);
	se_block[17] = evalStatisticalInefficiencyBlockAtB(r_vec_meaned, N_tot - N_eq, 40000);
	
	
	printf("U=%.0f, P=%.2f, r=%.2f\n", meanE, meanP, meanr);
	
	free(E_m_vec);
	free(P_vec);
	free(r_vec);
	free(E_m_vec_meaned);
	free(P_vec_meaned);
	free(r_vec_meaned);
}

void runTask2() {
	double E_CuCu = -436e-3; double E_ZnZn = -113e-3; double E_CuZn = -294e-3;
	int nUnitCellLengths = 10; int nAtoms = 2 * nUnitCellLengths * nUnitCellLengths* nUnitCellLengths; 
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r, 48);
	
	
	int N = 4000000; int N_eq = 1000000; int N_tot = N_eq + N;
	double T = 200; double deltaT = 4; int nT = 201; double *UCPr;  //double deltaT = 2; int nT = 501;
	double *T_vec; double *U_vec; double *C_vec; double *P_vec; double *r_vec; double *P_mean_vec; double *r_mean_vec;
	
	int t;
	UCPr = malloc(6 * sizeof(double));
	T_vec = malloc(nT * sizeof(double));
	U_vec = malloc(nT * sizeof(double));
	C_vec = malloc(nT * sizeof(double));
	P_vec = malloc(nT * sizeof(double));
	r_vec = malloc(nT * sizeof(double));
	P_mean_vec = malloc(nT * sizeof(double));
	r_mean_vec = malloc(nT * sizeof(double));
	
	double *se_corr; double *s_E_m_corr_vec; double *s_P_corr_vec; double *s_r_corr_vec;
	se_corr = malloc(3 * sizeof(double));
	s_E_m_corr_vec = malloc(nT * sizeof(double));
	s_P_corr_vec = malloc(nT * sizeof(double));
	s_r_corr_vec = malloc(nT * sizeof(double));
	
	double *se_block; double (*s_E_m_block_vec)[6]; double (*s_P_block_vec)[6]; double (*s_r_block_vec)[6]; 
	se_block = malloc(18 * sizeof(double));
	s_E_m_block_vec = malloc(nT * sizeof *s_E_m_block_vec);
	s_P_block_vec = malloc(nT * sizeof *s_P_block_vec);
	s_r_block_vec = malloc(nT * sizeof *s_r_block_vec);
	
	double *var_E_m_vec; double *var_P_vec; double *var_r_vec; double *var;
	var_E_m_vec = malloc(nT * sizeof (double));
	var_P_vec = malloc(nT * sizeof (double));
	var_r_vec = malloc(nT * sizeof (double));
	var = malloc(3 * sizeof(double));
	
	int (*neighbourMatrix)[8]; 
	neighbourMatrix = malloc(nAtoms * sizeof *neighbourMatrix);
	
	initializeNeighbourMatrix(nUnitCellLengths, neighbourMatrix, nAtoms);
	
	int j;
	for (t=0; t<nT; t++){
		printf("T=%.0f\n", T);
		T_vec[t] = T;
		metropolis(T, nAtoms, nUnitCellLengths, neighbourMatrix, N_tot, N_eq, E_CuCu, E_ZnZn, E_CuZn, r, UCPr, se_corr, se_block, var);

		printf("\n");
		U_vec[t] = UCPr[0];
		C_vec[t] = UCPr[1];
		P_vec[t] = UCPr[2];
		r_vec[t] = UCPr[3];
		P_mean_vec[t] = UCPr[4];
		r_mean_vec[t] = UCPr[5];
		
		s_E_m_corr_vec[t] = se_corr[0];
		s_P_corr_vec[t] = se_corr[1];
		s_r_corr_vec[t] = se_corr[2];
		
		for(j=0; j<6; j++) {
			s_E_m_block_vec[t][j] = se_block[j];
			s_P_block_vec[t][j] = se_block[j+6];
			s_r_block_vec[t][j] = se_block[j+12];
		}
		
		var_E_m_vec[t] = var[0];
		var_P_vec[t] = var[1];
		var_r_vec[t] = var[2];
		
		T += deltaT;
	}
	
	saveDataToFile("2/TU.csv", U_vec, T_vec, nT, 1);
	saveDataToFile("2/TC.csv", C_vec, T_vec, nT, 1);
	saveDataToFile("2/TP.csv", P_vec, T_vec, nT, 1);
	saveDataToFile("2/Tr.csv", r_vec, T_vec, nT, 1);
	saveDataToFile("2/TPmean.csv", P_mean_vec, T_vec, nT, 1);
	saveDataToFile("2/Trmean.csv", r_mean_vec, T_vec, nT, 1);
	
	saveDataToFile("2/s_E_m_corr.csv", s_E_m_corr_vec, T_vec, nT, 1);
	saveDataToFile("2/s_P_corr.csv", s_P_corr_vec, T_vec, nT, 1);
	saveDataToFile("2/s_r_corr.csv", s_r_corr_vec, T_vec, nT, 1);
	
	saveDataBlockToFile("2/s_E_m_block.csv", s_E_m_block_vec, T_vec, nT, 1);
	saveDataBlockToFile("2/s_P_block.csv", s_P_block_vec, T_vec, nT, 1);
	saveDataBlockToFile("2/s_r_block.csv", s_r_block_vec, T_vec, nT, 1);
	
	saveDataToFile("2/var_E_m.csv", var_E_m_vec, T_vec, nT, 1);
	saveDataToFile("2/var_P.csv", var_P_vec, T_vec, nT, 1);
	saveDataToFile("2/var_r.csv", var_r_vec, T_vec, nT, 1);
	
	free(neighbourMatrix);
	free(UCPr);
	free(se_corr); free(s_E_m_corr_vec); free(s_P_corr_vec); free(s_r_corr_vec);
	free(se_block); free(s_E_m_block_vec); free(s_P_block_vec); free(s_r_block_vec);
	free(var_E_m_vec); free(var_P_vec); free(var_r_vec); free(var);
	free(T_vec);
	free(U_vec);
	free(C_vec);
	free(P_vec);
	free(r_vec);
	free(P_mean_vec);
	free(r_mean_vec);
}

int main() {
	//runTask1();
	runTask2();
	
	return 0;
}