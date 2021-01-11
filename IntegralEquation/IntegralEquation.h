#ifndef INTEGRAL_EQUATION
#define INTEGRAL_EQUATION
#include "./../auxelirable.h"
#include "kernel_funcs.h"
#include "./../memory.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265359

typedef struct Grid {
	//Integral domain
	long double *S;
	//Observe domain
	long double *X;
} Grid;

typedef struct IntegralEquation { //AZ = U	
	long double *U;
	long double *Z;
	//Matrix A is an approximation of integral operator.
	//Dimensions
	int N; //N - Cols  
	int M; //M - Rows
	long double **A;
} IntegralEquation;

//long double ** allocate_memory_for_matrix_floats(int rows, int cols);

//void free_matrix_memory(long double **Matrix);

//Prototype =============================================================
int one_D_minimization(struct IntegralEquation *InEq, 
	long double G[], long double H[], 
	long double *al, long double almax, long double alpha, long double ro);
//======================================================================

void InitGridAndEquation(Grid *Grid, IntegralEquation *InEq, long double a, long double b, long double c, long double d, int N, int M);


void InitGridAndEquation1(Grid *Grid, IntegralEquation *InEq, long double a, long double b, long double c, long double d, int N, int M);

void destructor(Grid *Grid, IntegralEquation *InEq);

void MapAandZToPIPlus(IntegralEquation *InEq);

long double* gradient_norm(IntegralEquation *InEq, long double *U);

void stab_grad_weighted(long double *G, long double *Z, int N, long double alpha, long double ro);

void form_active_constr(IntegralEquation *InEq, int IPLUS[], long double G[], int *idim, int *idimo, int *ich/*, int *jch*/);

void form_new_desc_dir(long double H[], long double G[], int IPLUS[], long double *angri, long double *angri0, int *ich, int N);

long double search_of_max_step(IntegralEquation *InEq, long double H[], int *iprim, int N);

void PTISR1(struct IntegralEquation *InEq, 
	long double *Z0,//size = N
	long double alpha,
	long double ro,
	long double dl2);

int one_D_minimization(struct IntegralEquation *InEq, 
	long double G[], long double H[], 
	long double *al, long double almax, long double alpha, long double ro);

long double* from_pi_plus_to_initial_space(long double Z[], int N, long double Y[]);

#endif

