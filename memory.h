#ifndef MEMORY
#define MEMORY
#include <stdlib.h>
    
long double ** allocate_memory_for_matrix_floats(int rows, int cols);

void free_matrix_memory(long double **Matrix);

#endif
