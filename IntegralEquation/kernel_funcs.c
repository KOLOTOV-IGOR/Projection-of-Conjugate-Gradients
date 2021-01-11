#include "kernel_funcs.h"

long double kernel(long double x, long double s) {
	return 1/(1 + 100*pow(x-s,2));
}

