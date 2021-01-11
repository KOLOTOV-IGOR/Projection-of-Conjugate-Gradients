#ifndef AUXELIRABLE
#define AUXELIRABLE
//#include "IntegralEquation/IntegralEquation.h"
#include <stdlib.h>
#include <math.h>

//Auxelirable functions ============================================================================
long double* descending_direction(long double *a, long double *b, long double p, int N);

void descending_direction_v2(long double a[], long double b[], long double p, int N);

long double* matrix_vector_mult(long double **A, long double *Z, int M/*rows*/, int N/*cols*/);

long double scalar_product(long double *A, long double *B, int size);

long double weighted_scalar_prod(long double *A, long double *B, int *weights, int size);

long double discrepancy(long double A[], long double B[], int size);

long double stabilizer_calc(long double Z[], int size, long double discr, long double alpha, long double ro);

void reinitialize_of_vector_real(long double Z[], int size, long double value);

void reinitialize_of_vector_int(int Z[], int size, int value);

#endif
