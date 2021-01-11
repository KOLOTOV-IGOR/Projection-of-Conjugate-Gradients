#include "memory.h"
    
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


