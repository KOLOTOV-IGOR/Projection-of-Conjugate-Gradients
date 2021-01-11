#include "auxelirable.h"

//Auxelirable functions ============================================================================
long double* descending_direction(long double *a, long double *b, long double p, int N) { //PTICR1
	long double *descDirecrion = (long double*)malloc(N*sizeof(long double));
	for (int i = 0; i < N; i++) {
		descDirecrion[i] = a[i] + p*b[i];
	}

	return descDirecrion;
}

void descending_direction_v2(long double a[], long double b[], long double p, int N) { //PTICR1 variant with adding to first vector 
	for (int i = 0; i < N; i++) {
		a[i] += p*b[i];
	}
}

long double* matrix_vector_mult(long double **A, long double *Z, int M/*rows*/, int N/*cols*/) {//PTICR3
	long double *U = (long double*)malloc(M*sizeof(long double));
	for (int i = 0; i < M; i++) {
		U[i] = 0;
		for (int j = 0; j < N; j++) {
			U[i] += A[i][j]*Z[j];
		}
	}
	
	return U;
}

long double scalar_product(long double *A, long double *B, int size) {//PTICR6
	long double s = 0;
	for (int i = 0; i < size; i++) {
		s += A[i]*B[i];
	}

	return s;
}

long double weighted_scalar_prod(long double *A, long double *B, int *weights, int size) {//PTICR7
	long double s = 0;
	for (int i = 0; i < size; i++) {
		s += A[i]*B[i]*weights[i];
	}

	return s;
}

long double discrepancy(long double A[], long double B[], int size) { //PTICR5
	long double s = 0;
	for (int i = 0; i < size; i++) {
		s += pow(A[i] - B[i], 2);
	}
	
	return s;
}

long double stabilizer_calc(long double Z[], int size, long double discr, long double alpha, long double ro) { //PTICR9
	if (alpha == 0) {
		return discr;
	}

	long double s = 0;
	for (int i = 1; i < size; i++) {
		s += ro*pow(Z[i] - Z[i-1], 2) + pow(Z[i], 2);
	}
	s += pow(Z[0],2);

	return discr + alpha*s;
}

void reinitialize_of_vector_real(long double Z[], int size, long double value) {//PTICR2
	for (int i = 0; i < size; i++) {
		Z[i] = value;
	}
}

void reinitialize_of_vector_int(int Z[], int size, int value) {//PTICI2
	for (int i = 0; i < size; i++) {
		Z[i] = value;
	}
}

