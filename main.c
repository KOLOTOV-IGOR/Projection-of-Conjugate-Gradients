#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "IntegralEquation/IntegralEquation.h"
//#include "memory.h"
#include <time.h> 

void copy(long double *A, long double *B, int size) {
	for (int i = 0; i < size; i++) {
		B[i] = A[i];
	}
	free(A);
}

int main() {
	//clock_t t; 
        //t = clock();
	Grid Grid;
	IntegralEquation InEq;

	InitGridAndEquation(&Grid, &InEq, 0, 1, -2, 2, 41, 60);
	//InitGridAndEquation(&Grid, &InEq, 0, 1, -2, 2, 5000, 5100);
	//for (int i = 0; i < InEq.M; i++) {
	//	for (int j = 0; j < InEq.N; j++) {
	//		printf("%LF ", InEq.A[i][j]);
	//	}
	//	printf("\n");
	//}
	MapAandZToPIPlus(&InEq);

	//t = clock() - t; 
    	//double time_taken = ((double)t)/CLOCKS_PER_SEC;
	//printf("%f s\n", time_taken);

	long double *H = descending_direction(InEq.Z, InEq.Z, 0.0, InEq.N);
	//for (int i = 0; i < InEq.M; i++) {
	//	for (int j = 0; j < InEq.N; j++) {
	//		printf("%LF ", InEq.A[i][j]);
	//	}
	//	printf("\n");
	//}
	PTISR1(&InEq, H, 0.0, 0.0, 0.0);
	free(H);
	H = descending_direction(InEq.Z, InEq.Z, 0.0, InEq.N);
	copy(from_pi_plus_to_initial_space(InEq.Z, InEq.N, H)/*H*/, InEq.Z, InEq.N);

	for (int i = 0; i < InEq.N; i++) {  
		printf("-----+> %LE\n", InEq.Z[i]);
	}
	FILE *fp;
	fp = fopen("vals_test.txt", "w+");
	for (int i = 0; i < InEq.N; i++) {  
		//printf("%LF -----+> %LE\n", Grid.S[i], InEq.Z[i]);
		fprintf(fp, "%LE\n", InEq.Z[i]);
	}
	fclose(fp);

	destructor(&Grid, &InEq);

	return 0;
}
