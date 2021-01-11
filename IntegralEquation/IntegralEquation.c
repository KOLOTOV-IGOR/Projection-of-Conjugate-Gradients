#include "IntegralEquation.h"
/*
long double ** allocate_memory_for_matrix_floats(int rows, int cols) {
	long double **Matrix = (long double**)malloc(rows*sizeof(long double));
	Matrix[0] = (long double*)malloc(rows*cols*sizeof(long double)); 
	
	for (int i = 1; i < rows; ++i) {
		Matrix[i] = Matrix[i-1] + cols;
	}
	
	return Matrix;
}

void free_matrix_memory(long double **Matrix) {
	free(Matrix[0]);
	free(Matrix);
}
*/


void InitGridAndEquation(Grid *Grid, IntegralEquation *InEq, long double a, long double b, long double c, long double d, int N, int M) { //Not Parallelizable
	//Grid initialization
	Grid->X = (long double*)malloc(M*sizeof(long double)); 
	Grid->S = (long double*)malloc(N*sizeof(long double));
	//Equation initialization
	InEq->U = (long double*)malloc(M*sizeof(long double)); 
	InEq->Z = (long double*)malloc(N*sizeof(long double));
	InEq->N = N; InEq->M = M;

	long double hx = (d-c)/(M-1);
	for (int i = 0; i < M; i++) {
		InEq->U[i] = 0.0;
		Grid->X[i] = c + hx*i;
	}
	long double hs = (b-a)/(N-1);
	for (int j = 0; j < N; j++) {
		Grid->S[j] = a + hs*j;
		InEq->Z[j] = 4.0*Grid->S[j]*(1-Grid->S[j]);
	}

	//Aproximate Operator by Matrix A
	InEq->A = allocate_memory_for_matrix_floats(M, N);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			InEq->A[i][j] = hs*kernel(Grid->X[i], Grid->S[j]);
			if (j == 0 || j == N-1) {
				InEq->A[i][j] = 0.5*hs*kernel(Grid->X[i], Grid->S[j]);
			}
			InEq->U[i] += InEq->A[i][j]*InEq->Z[j];
		}
	}

	for (int j = 0; j < N; j++) {
		InEq->Z[j] = 0.0;
	}

	return;
}

void InitGridAndEquation1(Grid *Grid, IntegralEquation *InEq, long double a, long double b, long double c, long double d, int N, int M) { //Not Parallelizable
	//Grid initialization
	Grid->X = (long double*)malloc(M*sizeof(long double)); 
	Grid->S = (long double*)malloc(N*sizeof(long double));
	//Equation initialization
	InEq->U = (long double*)malloc(M*sizeof(long double)); 
	InEq->Z = (long double*)malloc(N*sizeof(long double));
	InEq->N = N; InEq->M = M;

	long double hx = (d-c)/(M-1);
	for (int i = 0; i < M; i++) {
		InEq->U[i] = 0.0;
		Grid->X[i] = c + hx*i;
	}
	long double hs = (b-a)/(N-1);
	for (int j = 0; j < N; j++) {
		Grid->S[j] = a + hs*j;
		InEq->Z[j] = pow(sin(2*PI*Grid->S[j]), 2);
		//printf("+++ %LF --- %LF === %LF\n", Grid->S[j], InEq->Z[j], (long double)8*pow(PI,2)*cos(4*PI*Grid->S[j]));
	}

	//Aproximate Operator by Matrix A
	InEq->A = allocate_memory_for_matrix_floats(M, N);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			InEq->A[i][j] = hs*kernel(Grid->X[i], Grid->S[j]);
			if (j == 0 || j == N-1) {
				InEq->A[i][j] = 0.5*hs*kernel(Grid->X[i], Grid->S[j]);
			}
			InEq->U[i] += InEq->A[i][j]*InEq->Z[j];
		}
	}

	for (int j = 0; j < N; j++) {
		InEq->Z[j] = 0.0;
	}


	return;
}

void destructor(Grid *Grid, IntegralEquation *InEq) {
	//1. Free memory for Grid 
	free(Grid->X); free(Grid->S);

	//2. Free memory for solution, right part and other 1D arr
	free(InEq->U); free(InEq->Z);

	//3.Free memory for Operator Matrix
	free_matrix_memory(InEq->A);
		
	return;
}

void MapAandZToPIPlus(IntegralEquation *InEq) {//PTISR2
	int N = InEq->N;
	int M = InEq->M;
	//long double *S = (long double*)malloc(N*sizeof(long double));
	long double S[N];

	for (int i = 1; i <= N-2; i++) {
		S[i] = (InEq->Z[i-1]-2*InEq->Z[i]+InEq->Z[i+1])*(long double)(N-i-1)*i/(1-N);	
	}
	S[0] = InEq->Z[0]; S[N-1] = InEq->Z[N-1];
	
	long double *Z = descending_direction(S, S, 0.0, N);
	free(InEq->Z); InEq->Z = Z;

	long double T; //Element of transformation matrix
	for (int k = 0; k < M; k++) {
		for (int j = 0; j < N; j++) {
			S[j] = 0.0;
			for (int i = 0; i < N; i++) {
				//PTISR4=====================================
				if (i <= j && j != 0) {
					T = (long double)(i)/(j);
				} else if (i >= j && j != N-1) {
					T = (long double)(N-i-1)/(N-j-1);
				}	
				//printf("%LF\n", T);
				//===========================================
				S[j] += InEq->A[k][i]*T;
			}
			//printf("\n");
			//InEq->A[k][j] = S[j]; //WRONG!!!
		}
		for (int i = 0; i < N; i++) {
			InEq->A[k][i] = S[i];
		}
	}

	return;
}

long double* gradient_norm(IntegralEquation *InEq, long double *U) { //PTICR4
	int N = InEq->N;
	int M = InEq->M;
	long double *G = (long double*)malloc(N*sizeof(long double)); 

	for (int j = 0; j < N; j++) {
		long double s = 0;
		for (int i = 0; i < M; i++) {
			s += InEq->A[i][j]*(U[i] - InEq->U[i]);	
		}
		G[j] = s*2;
	}

	return G;
}

void stab_grad_weighted(long double *G, long double *Z, int N, long double alpha, long double ro) {//PTICR8
	if (alpha == 0) {
		return;
	}
	for (int i = 0; i < N; i++) {
		G[i] += 2.0*alpha*Z[i];
		if (i == 0) {
			G[i] += 2*alpha*(Z[0] - Z[1])*ro;
		} else if (i == N-1) {
			G[i] += 2*alpha*(Z[N-1] - Z[N-2])*ro;
		} else {
			G[i] += 2*alpha*(2*Z[i] - Z[i-1] - Z[i+1])*ro;
		}
	}

	return;
}

void form_active_constr(IntegralEquation *InEq, int IPLUS[], long double G[], int *idim, int *idimo, int *ich/*, int *jch*/) {
	int N = InEq->N;
	//int M = InEq->M;
	printf("Z[N-1] = %LE --- G[N-1] = %LE\n", InEq->Z[N-1], G[N-1]);
	printf("Z[0] = %LE --- G[0] = %LE\n", InEq->Z[0], G[0]);
	for (int i = 0; i < N; i++) {
		if (InEq->Z[i] > 0 || G[i] < 0) {
			continue;
		}
		IPLUS[i] = 0;
		InEq->Z[i] = 0.0;
	}
	*idim = 0; //size of current verge(грань)	
	//*idimo = -1; //Error!
	for (int i = 0; i < N; i++) {
		*idim += IPLUS[i];
	}
	*ich = 1;

	return;
}

void form_new_desc_dir(long double H[], long double G[], int IPLUS[], long double *angri, long double *angri0, int *ich, int N) {
	long double bt = 0.0;
	if (*ich == 0) {
		//printf("angri0 = %LE\n", *angri0);
		bt = (*angri)/(*angri0);
	}
	*angri0 = *angri;
	for (int i = 0; i < N; i++) {
		H[i] = (bt*H[i] + G[i])*IPLUS[i];
	}
	
	return;
}

long double search_of_max_step(IntegralEquation *InEq, long double H[], int *iprim, int N) {
	long double almax = 1.0E+18;
	*iprim = 0;
	for (int i = 0; i < N; i++) {
		long double al;
		if (H[i] <= 0.0) {
			continue;
		}
		al = InEq->Z[i]/H[i]; 
		if (almax < al) {
			continue;
		}
		almax = al;
		*iprim = i; 
	}

	return almax;
}

void PTISR1(struct IntegralEquation *InEq, 
	long double *Z0,//size = N
	long double alpha,
	long double ro,
	long double dl2) {

	long double stab, stab1;
	int N = InEq->N; int M = InEq->M;
	int iterbn, iend, iprim;
	long double angrd, angri, angri0, al;
	long double *H = (long double*)malloc(N*sizeof(long double)); 
	reinitialize_of_vector_real(H, InEq->N, 0.0);

	free(InEq->Z);
	//Preparation =============================================================
	InEq->Z = descending_direction(Z0, Z0, 0.0, InEq->N);
	long double *U = matrix_vector_mult(InEq->A, InEq->Z, InEq->M, InEq->N);
	long double discr = discrepancy(U, InEq->U, InEq->M);
	stab = stabilizer_calc(InEq->Z, InEq->N, discr, alpha, ro);
	long double *G = gradient_norm(InEq, U);
	stab_grad_weighted(G, InEq->Z, InEq->N, alpha, ro);
	//=========================================================================
	
	int *IPLUS = (int*)malloc(N*sizeof(int)); 
	//reinitialize_of_vector_int(IPLUS, InEq->N, 1);
	iterbn = -1; iend = 0;//iterbn = iend = 0;
	angrd/* = angri0 */= 0.0;
	int jch = -1;
	int idim; //size of current verge(грань)	
	int idimo = -1; 
	int iter_max = 8000;
	long double almax;

	for (int iter = 0; iter < iter_max; iter++) {
		//for (int i = 0; i < N; i++) {  
			//printf("------ %d \n", idimo);
		//}
		//printf("jch = %d\n", jch);
		//printf("iter = %d\n", iter);
		int ich = 0;//ich - set changed flag. If ich = 1 active constraint changed
		if (jch == -1) {
			//1. Forming of the set of active constraints 
			reinitialize_of_vector_int(IPLUS, InEq->N, 1);
			form_active_constr(InEq, IPLUS, G, &idim, &idimo, &ich);
		} else if (jch == 1) {
			IPLUS[iprim] = 0;
			InEq->Z[iprim] =0.0;
			ich = 1;
			idim -= 1;	
		}
		angri = weighted_scalar_prod(G, G, IPLUS, N);
		iterbn++;
		printf("iterbn = %d --- idim = %d --- idim0 = %d --- ich = %d --- jch = %d\n", iterbn, idim, idimo, ich, jch);
		if (ich == 1) {
			iterbn = 0;//1;
		}
		
		if (angri > angrd && iterbn != idim+1 && iend != 1) {
			//2.Forming new descending direction
			form_new_desc_dir(H, G, IPLUS, &angri, &angri0, &ich, N);
			//3.Searching of max step
			almax = search_of_max_step(InEq, H, &iprim, N);

			int ied = one_D_minimization(InEq, G, H, &al, almax, alpha, ro);
			descending_direction_v2(InEq->Z, H, -al, N);
			free(U); free(G);
			U = matrix_vector_mult(InEq->A, InEq->Z, M, N);
			G = gradient_norm(InEq, U);
			stab_grad_weighted(G, InEq->Z, N, alpha, ro);
			angri = weighted_scalar_prod(G, G, IPLUS, N);
			discr = discrepancy(U, InEq->U, M);
			stab1 = stabilizer_calc(InEq->Z, N, discr, alpha, ro);
			jch = 0;
			if (ied == 1) {
				jch = 1;
			} else if (ied == 2) { //Шаг отрицателен. Ошибка.
				printf("Шаг отрицателен. Ошибка\n");
				break;
			} //else 
			if (discr <= dl2 || iter >= iter_max) {//Выход по невязке или числу итераций
				printf("Выход по невязке или числу итераций\n");
				break;
			} else if (stab1 >= stab) {//Невязка не уменьшилась
				printf("Невязка не уменьшилась\n");
				break;
			} 
			stab = stab1;
			continue;
		} //else if (iend == 1 && idim == idimo) {//Идёт проверка на достижение точного минимума.
		//	break; 
		//} else {
			///////////////////////////////////////////////////////////////////////////////
			if (iend == 1 && idim == idimo) {//Идёт проверка на достижение точного минимума.
				printf("точный минимум\n");
				break;
			}
			///////////////////////////////////////////////////////////////////////////////
			iend = 1 - iend;
			idimo = idim;
			if (iend == 0) {
				//2.Forming new descending direction
				form_new_desc_dir(H, G, IPLUS, &angri, &angri0, &ich, N);
				//3.Searching of max step
				almax = search_of_max_step(InEq, H, &iprim, N);

				int ied = one_D_minimization(InEq, G, H, &al, almax, alpha, ro);
				descending_direction_v2(InEq->Z, H, -al, N);
				free(U); free(G);
				U = matrix_vector_mult(InEq->A, InEq->Z, M, N);
				G = gradient_norm(InEq, U);
				stab_grad_weighted(G, InEq->Z, N, alpha, ro);
				angri = weighted_scalar_prod(G, G, IPLUS, N);
				discr = discrepancy(U, InEq->U, M);
				stab1 = stabilizer_calc(InEq->Z, N, discr, alpha, ro);
				jch = 0;
				if (ied == 1) {
					jch = 1;
				} else if (ied == 2) {
					printf("Шаг отрицателен. Ошибка(1)\n");
					break;
				} //else 
				if (discr <= dl2 || iter >= iter_max) {
					printf("Выход по невязке или числу итераций(1)\n");
					break;
				} else if (stab1 >= stab) {
					printf("Невязка не уменьшилась(1)\n");
					break;
				} 
				stab = stab1;
				continue;
			}
			jch = -1;
			iter -= 1;
		//}
	}

	for (int i = 0; i < N; i++) {  
		printf(">>>>> %LE\n", InEq->Z[i]);
	}

	printf("......> %LE\n", discr);
	free(U); free(G); free(IPLUS); free(H);
}

int one_D_minimization(struct IntegralEquation *InEq, 
	long double G[], long double H[], 
	long double *al, long double almax, long double alpha, long double ro) { //PTICRO

	long double h2 = 0;
	int N = InEq->N;
	int M = InEq->M;
	long double sp1 = scalar_product(G, H, N);

	if (sp1 < 0.0) {
		//The step is negative or equal zero!
		*al = 0; 
		return 2;
	}
	long double *U = matrix_vector_mult(InEq->A, H, M/*rows*/, N/*cols*/);
	long double sp2 = scalar_product(U, U, M);
	free(U);

	if (alpha == 0.0) {
		if (sp2 + h2*alpha == 0.0) {
			*al = 0;
			return 3;
		}
		*al = sp1/(sp2 + h2*alpha )*0.5;
		if (*al >= almax) {
			*al = almax;
			return 1;
		}
		return 0;
	}
	//alpha != 0.0
	for (int i = 1; i < N; i++) {
		h2 += pow(H[i] - H[i-1], 2)*ro + pow(H[i], 2);
	}
	h2 += pow(H[0], 2);
	if (sp2 + h2*alpha == 0.0) {
		*al = 0;
		return 3;
	}
	*al = sp1/(sp2 + h2*alpha )*0.5;
	if (*al >= almax) {
		*al = almax;
		return 1;
	}

	return 0;
}

long double* from_pi_plus_to_initial_space(long double Z[], int N, long double Y[]) { //PTISR3
	long double *res;
	for (int i = 0; i < N; i++) {
		long double S = 0;
		long double T;
		for (int j = 0; j < N; j++) {
			//PTISR4=====================================
			if (i <= j && j != 0) {
				T = (long double)(i)/(j);
			} else if (i >= j && j != N-1) {
				T = (long double)(N-i-1)/(N-j-1);
			}	
			//===========================================
			S += Z[j]*T;
		}
		Y[i] = S;
	}
	res = descending_direction(Y, Y, 0.0, N);

	return res;
}


