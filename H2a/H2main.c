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
	double T = -100+273; double deltaT = 2; int nT = 451; 
	double *T_vec; double *P_vec; double *U_vec;
	
	T_vec = malloc( nT * sizeof(double));
	P_vec = malloc( nT * sizeof(double));
	U_vec = malloc( nT * sizeof(double));
	for (int t=0; t<nT; t++){
		T_vec[t] = T - 273.15;
		P_vec[t] = minFreeEnergBisection(T, E_CuCu, E_ZnZn, E_CuZn, N, tol);
		U_vec[t] = evalU(P_vec[t], E_CuCu, E_ZnZn, E_CuZn, N);
		T += deltaT;
	}
	
	saveDataToFile("1/TP.csv", P_vec, T_vec, nT, 1);
	saveDataToFile("1/TU.csv", U_vec, T_vec, nT, 1);
		
	free(T_vec);
	free(P_vec);
	free(U_vec);
}

void initializeNeighbourMatrix(int nUnitCellLengths, int (*neighbourMatrix)[8], int *sc1, int *sc2){
	int sc1_mat3d[nUnitCellLengths+1][nUnitCellLengths+1][nUnitCellLengths+1]; 
	int sc2_mat3d[nUnitCellLengths+1][nUnitCellLengths+1][nUnitCellLengths+1];
	
	/*
	 * Building first of two sc-lattices with nUnitCellLengths^3 unit cells
	 * Giving vertex in the lattice a index, periodicity included 
	 * Putting the indexes in either sc-lattice into the corresponding sc-list
	 */
	int i; int j; int k; int index=0; int iSc = 0;
	for (i=0; i<nUnitCellLengths; i++) {
		for (j=0; j<nUnitCellLengths; j++) {
			for (k=0; k<nUnitCellLengths; k++) {
				sc1_mat3d[i][j][k] = index;
				sc1[iSc] = index;
				index++;
				iSc++;
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
	iSc = 0;
	for (i=1; i<nUnitCellLengths+1; i++) {
		for (j=1; j<nUnitCellLengths+1; j++) {
			for (k=1; k<nUnitCellLengths+1; k++) {
				sc2_mat3d[i][j][k] = index;
				sc2[iSc] = index;
				index++;
				iSc++;
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
	
	/*for (i=0; i<nUnitCellLengths+1; i++) {
		for (j=0; j<nUnitCellLengths+1; j++) {
			for (k=0; k<nUnitCellLengths+1; k++) {
				printf("%i   ", sc1_mat3d[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	printf("\n\n\n");
	for (i=0; i<nUnitCellLengths+1; i++) {
		for (j=0; j<nUnitCellLengths+1; j++) {
			for (k=0; k<nUnitCellLengths+1; k++) {
				printf("%i   ", sc2_mat3d[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n\n");
	}
	//printf("index=%i\n", index);
	for (i=0; i<indexMatris; i++) {
		for (j=0; j<8; j++){
			printf("%i ", neighbourMatrix[i][j]);
		}
		printf("\n\n");
	}*/
	
}

void newConfiguration(int (*neighbourMatrix)[8], int (*newNeighbourMatrix)[8], int *sc1, int *sc2, int *newSc1, int *newSc2, int nAtoms, gsl_rng * r){
	/*int iSc1 = floor(nAtoms/2 * gsl_rng_uniform(r)); int iSc2 = floor(nAtoms/2 * gsl_rng_uniform(r));
	int atom1 = sc1[iSc1]; int atom2 = sc2[iSc2];
	
	if (atom1 < nAtoms/2) {
		while (atom2 < nAtoms/2) {
			iSc2 = floor(nAtoms/2 * gsl_rng_uniform(r));
			atom2 = sc2[iSc2];
		}
	} else {
		while (atom2 > nAtoms/2) {
			iSc2 = floor(nAtoms/2 * gsl_rng_uniform(r));
			atom2 = sc2[iSc2];
		}
	}
	int i; int j; 
	
	for (i=0; i<nAtoms/2; i++){
		newSc1[i] = sc1[i];
		newSc2[i] = sc2[i];
	}
	newSc1[iSc1] = atom2; newSc2[iSc2] = atom1;*/
	int i; int j; 
	int atom1 = floor(nAtoms/2 * gsl_rng_uniform(r)); int atom2 = floor(nAtoms/2 * gsl_rng_uniform(r)) + nAtoms/2;
	for (i=0; i<nAtoms/2; i++){
		if (sc1[i] == atom1) {
			newSc1[i] = atom2;
		} else if (sc1[i] == atom2) {
			newSc1[i] = atom1;
		} else {
			newSc1[i] = sc1[i];
		}

		if (sc2[i] == atom1) {
			newSc2[i] = atom2;
		} else if (sc2[i] == atom2) {
			newSc2[i] = atom1;
		} else {
			newSc2[i] = sc2[i];
		}
	}
	
	for (i=0; i<nAtoms; i++){
		if (i == atom1) {
			for (j=0; j<8; j++) {
				if (neighbourMatrix[atom2][j] == atom1) {
					newNeighbourMatrix[atom1][j] = atom2;
				} else {
					newNeighbourMatrix[atom1][j] = neighbourMatrix[atom2][j]; 
				}
			}
		} else if (i == atom2) {
			for (j=0; j<8; j++) {
				if (neighbourMatrix[atom1][j] == atom2) {
					newNeighbourMatrix[atom2][j] = atom1;
				} else {
					newNeighbourMatrix[atom2][j] = neighbourMatrix[atom1][j]; 
				}
			}
		} else {
			for (j=0; j<8; j++) {
				if (neighbourMatrix[i][j] == atom1) {
					newNeighbourMatrix[i][j] = atom2;
				} else if (neighbourMatrix[i][j] == atom2) {
					newNeighbourMatrix[i][j] = atom1;
				} else {
					newNeighbourMatrix[i][j] = neighbourMatrix[i][j];
				}
			}
		}
	}
}

double evalEnergyForState(int (*neighbourMatrix)[8], int nAtoms, double E_AA, double E_BB, double E_AB){
	double E = 0; int i; int j;
	for(i=0; i<nAtoms/2; i++) {
		for(j=0; j<8; j++) {
			if(neighbourMatrix[i][j] < nAtoms/2) {
				E += E_AA;
			} else {
				E += E_AB;
			}
		}
	}
	for(i=nAtoms/2; i<nAtoms; i++) {
		for(j=0; j<8; j++) {
			if(neighbourMatrix[i][j] < nAtoms/2) {
				E += E_AB;
			} else {
				E += E_BB;
			}
		}
	}
	return E/2.0;
}

double evalLongRangeOrder(int *sc_a, int nAtoms) {
	int N_A = 0;
	for (int i=0; i<nAtoms/2; i++){
		if (sc_a[i] < nAtoms/2) {
			N_A++;
		}
	}
	return (double) N_A / nAtoms * 4.0 - 1.0; 
}

double evalShortRangeOrder(int (*neighbourMatrix)[8], int nAtoms){
	double q = 0; int i; int j; double N = (double) nAtoms/2.0;
	for(i=0; i<nAtoms/2; i++) {
		for(j=0; j<8; j++) {
			if(neighbourMatrix[i][j] >= nAtoms/2) {
				q++;
			}
		}
	}
	return 1/(4.0*N)*(q - 4.0*N);
}

void metropolis(double T, int nAtoms, int nUnitCellLengths, int N_tot, int N_eq, double E_AA, double E_BB, double E_AB, gsl_rng * r, double *UCPr){
	int (*neighbourMatrix_m)[8]; int (*neighbourMatrix_t)[8];
	int *sc1_m; int *sc2_m; int *sc1_t; int *sc2_t;
	
	neighbourMatrix_m = malloc(nAtoms * sizeof *neighbourMatrix_m);
	neighbourMatrix_t = malloc(nAtoms * sizeof *neighbourMatrix_t);
	sc1_m = malloc(nAtoms * sizeof(int));
	sc2_m = malloc(nAtoms * sizeof(int));
	sc1_t = malloc(nAtoms * sizeof(int));
	sc2_t = malloc(nAtoms * sizeof(int));
	
	initializeNeighbourMatrix(nUnitCellLengths, neighbourMatrix_m, sc1_m, sc2_m);
	double E_m = evalEnergyForState(neighbourMatrix_m, nAtoms, E_AA, E_BB, E_AB); double E_t;
	
	double beta = 1/(K_B * T);
	
	double meanE=0; double meanE2=0; double meanP=0; double meanr=0;
	
	int t; int i; int j;
	for (t=0; t<N_tot; t++) {
		newConfiguration(neighbourMatrix_m, neighbourMatrix_t, sc1_m, sc2_m, sc1_t, sc2_t, nAtoms, r);
		
		E_t = evalEnergyForState(neighbourMatrix_t, nAtoms, E_AA, E_BB, E_AB); 
		
		if ( gsl_rng_uniform(r) < exp(-(E_t - E_m) * beta) ) {
			for (i=0; i<nAtoms; i++) {
				for (j=0; j<8; j++) {
					neighbourMatrix_m[i][j] = neighbourMatrix_t[i][j];
				}
			}
			for (i=0; i<nAtoms/2; i++) {
				sc1_m[i] = sc1_t[i];
				sc2_m[i] = sc2_t[i];
			}
			E_m = E_t;
		}
		if ( t > N_eq ) {
			meanE += E_m;
			meanE2 += E_m * E_m;
			meanP += fabs(evalLongRangeOrder(sc1_m, nAtoms));
			meanr += fabs(evalShortRangeOrder(neighbourMatrix_m, nAtoms));
		}
	}
	meanE /= ((double) N_tot - N_eq);
	meanE2 /= ((double) N_tot - N_eq); 
	meanP  /= ((double) N_tot - N_eq);
	meanr  /= ((double) N_tot - N_eq);
	UCPr[0] = meanE;
	UCPr[1] = (meanE2 - meanE*meanE) / (beta*beta) * K_B;
	UCPr[2] = fabs(evalLongRangeOrder(sc1_m, nAtoms));
	UCPr[3] = fabs(evalShortRangeOrder(neighbourMatrix_m, nAtoms));
	UCPr[4] = meanP;
	UCPr[5] = meanr;
	
	free(neighbourMatrix_m);
	free(neighbourMatrix_t);
	free(sc1_m);
	free(sc2_m);
	free(sc1_t);
	free(sc2_t);
}

void runTask2Metropolis() {
	double E_CuCu = -436e-3; double E_ZnZn = -113e-3; double E_CuZn = -294e-3;
	int nUnitCellLengths = 10; int nAtoms = 2 * nUnitCellLengths * nUnitCellLengths* nUnitCellLengths; 
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r, 34);
	
	
	int N = 1000000; int N_eq = 750000; int N_tot = N_eq + N;
	double T = -200+273.15; double deltaT = 2; int nT = 501; double *UCPr; 
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
	
	for (t=0; t<nT; t++){
		T_vec[t] = T - 273.15;
		metropolis(T, nAtoms, nUnitCellLengths, N_tot, N_eq, E_CuCu, E_ZnZn, E_CuZn, r, UCPr);
		U_vec[t] = UCPr[0];
		C_vec[t] = UCPr[1];
		P_vec[t] = UCPr[2];
		r_vec[t] = UCPr[3];
		P_mean_vec[t] = UCPr[4];
		r_mean_vec[t] = UCPr[5];
		T += deltaT;
	}
	
	saveDataToFile("2/TU.csv", U_vec, T_vec, nT, 1);
	saveDataToFile("2/TC.csv", C_vec, T_vec, nT, 1);
	saveDataToFile("2/TP.csv", P_vec, T_vec, nT, 1);
	saveDataToFile("2/Tr.csv", r_vec, T_vec, nT, 1);
	saveDataToFile("2/TPmean.csv", P_mean_vec, T_vec, nT, 1);
	saveDataToFile("2/Trmean.csv", r_mean_vec, T_vec, nT, 1);
	free(UCPr);
	free(T_vec);
	free(U_vec);
	free(C_vec);
	free(P_vec);
	free(r_vec);
	free(P_mean_vec);
	free(r_mean_vec);
	/* double T = 800+273.15; int N = 1000000; int N_eq = 500000; int N_tot = N_eq + N;  double *UCPr;
	UCPr = malloc(4 * sizeof(double));
	metropolis(T, nAtoms, nUnitCellLengths, N_tot, N_eq, E_CuCu, E_ZnZn, E_CuZn, r, UCPr);
	// 50000; 
	printf("U=%.4e\n", UCPr[0]);
	//printf("C=%.2e\n", UCPr[1]);
	//printf("P=%.2e\n", UCPr[2]);
	//printf("r=%.2e\n", UCPr[3]);
	free(UCPr);*/
	
	
	/*int (*neighbourMatrix_m)[8]; int (*neighbourMatrix_t)[8];
	int *sc1_m; int *sc2_m; int *sc1_t; int *sc2_t;
	
	neighbourMatrix_m = malloc(nAtoms * sizeof *neighbourMatrix_m);
	neighbourMatrix_t = malloc(nAtoms * sizeof *neighbourMatrix_t);
	sc1_m = malloc(nAtoms * sizeof(int));
	sc2_m = malloc(nAtoms * sizeof(int));
	sc1_t = malloc(nAtoms * sizeof(int));
	sc2_t = malloc(nAtoms * sizeof(int));
	
	initializeNeighbourMatrix(nUnitCellLengths, neighbourMatrix_m, sc1_m, sc2_m);
	newConfiguration(neighbourMatrix_m, neighbourMatrix_t, sc1_m, sc2_m, sc1_t, sc2_t, nAtoms, r);
	
	free(neighbourMatrix_m);
	free(neighbourMatrix_t);
	free(sc1_m);
	free(sc2_m);
	free(sc1_t);
	free(sc2_t);*/
}
/*
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

double evalCorrelationFunction(double f[], int k, int N) {
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

int evalStatisticalInefficiency(double f[], int N, double phi_vec[]) {
	int k = 0; double phi = 100;
	while (phi > 0.135 || k > 1000) {
		k++; 
		phi = evalCorrelationFunction(f, k, N);
		phi_vec[k-1] = phi;
	}
	return k;
}

void task2StatisticalInefficiency() {
	int N = 1401; int N_max = 1000; double *T; double s_vec[N_max]; double B_vec[N_max]; double phi_vec[N_max];
	double *U; double *U_meaned; 
	T = malloc((N) * sizeof(double));
	U = malloc((N) * sizeof(double));
	U_meaned = malloc((N) * sizeof(double));
	readDataFromFile("2/TU.csv", T, U);
	
	evalMeanedA(U, U_meaned, N);
	
	int s = evalStatisticalInefficiency(U_meaned, N, phi_vec);
	printf("s = %i\n", s);
	saveDataToFile1D("4/phi.csv", phi_vec, N_max);
	
	evalStatisticalInefficiencyBlock(U_meaned, N, s_vec, B_vec);
	saveDataToFile("4/block.csv", s_vec, B_vec, N_max, 1);
	
	free(T);
	free(U);
	free(U_meaned);
}
*/
int main() {
	//runTask1();
	runTask2Metropolis();
	
	return 0;
}