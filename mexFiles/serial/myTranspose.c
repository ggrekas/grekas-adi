#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
   int N, i, j;
   double temp;
   double *A;
   
   plhs[0] = prhs[0];
   
   A = mxGetPr(prhs[0]);
   N= mxGetM( prhs[0] );

	if( 1 == N)
       return;
   for( j= 0; j < N -1; ++j)
      for( i= j+1; i < N; ++i){
          temp = A[i+j*N];
          A[i+j*N] = A[i*N +j];
          A[i*N+ j] = temp;
      }
   return;
}