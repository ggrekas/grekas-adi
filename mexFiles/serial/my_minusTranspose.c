#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
   int N, i, j;
   double temp;
   double *A;
   
   plhs[0] = prhs[0];
   
   //plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
   A = mxGetPr(prhs[0]);
   N= mxGetM( prhs[0] );
     
   for( j= 0; j < N -1; ++j){
      A[j+j*N] = -A[j+j*N];
      for( i= j+1; i < N; ++i){
         temp = A[i+j*N];
         A[i+j*N] = -A[i*N +j];
         A[i*N+ j] = -temp;
      }
   }
   A[N*N-1] = -A[N*N-1];
   return;
   }