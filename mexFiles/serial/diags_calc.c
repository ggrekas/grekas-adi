#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"


typedef enum {IS_SCALAR, IS_MATRIX} varType;


varType find_var_type(mxArray *mAPtr);

void diags_calculation(mxArray *prhs[], double *diag_y_xSweep, double 
	*diag_x_xSweep, int N);
void initVariables(mxArray *prhs[], double **a, double **C, double *h2, double *k_t);

void d_calc_aMatrix_cScalar(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aScalar_cMatrix(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aMatrix_cMatrix(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aScalar_cScalar(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *diag_y_xSweep, *diag_x_xSweep;
	int N;

	if(mxGetM(prhs[0]) < mxGetM(prhs[1]) )
		N= (int)mxGetM(prhs[1]);
	else
		N= (int)mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
	diag_y_xSweep = mxGetPr(plhs[0]);
	diag_x_xSweep = mxGetPr(plhs[1]);
	
	diags_calculation(prhs, diag_y_xSweep, diag_x_xSweep, N);

	return;
}

void diags_calculation(mxArray *prhs[], double *diag_y_xSweep, double
	*diag_x_xSweep, int N){
	double *a, *C, h2, k_t;
	void (*fPtr[2][2])(double*, double*, double, double, double *, double *,
		int ) = {{d_calc_aScalar_cScalar, d_calc_aScalar_cMatrix},
		{d_calc_aMatrix_cScalar, d_calc_aMatrix_cMatrix}};
	
	initVariables( prhs, &a, &C, &h2, &k_t);
	fPtr[find_var_type(prhs[0])][find_var_type(prhs[1])](a, C, h2, k_t, diag_y_xSweep,
		diag_x_xSweep, N);
	
	
				
	return;
}


void initVariables(mxArray *prhs[], double **a, double **C, double *h2, double *k_t){
	*a = mxGetPr(prhs[0]);
	*C = mxGetPr(prhs[1]);
	*k_t = mxGetScalar(prhs[2]);
	*h2 = mxGetScalar(prhs[3]);
	
	return;
}



varType find_var_type(mxArray *mAPtr){
	if(1 == mxGetM(mAPtr))
		return IS_SCALAR;
		
	return IS_MATRIX;
}

void d_calc_aScalar_cMatrix(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j;
	double const_var1, const_var2, temp;
	
	const_var1 = k_t - *a/h2;
	const_var2 = k_t + *a/h2;
	for(i = 0; i <N; ++i)
		for(j = 0; j <N; ++j){
			temp = 0.25*C[i +j*N];
			diag_y_xSweep[i +j*N] =  const_var1 + temp;		
			diag_x_xSweep[i +j*N] =  const_var2 - temp;
		}
	return;
}

void d_calc_aMatrix_cScalar(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j;
	double const_var1, const_var2, temp;
	
	const_var1 = k_t + 0.25*(*C);
	const_var2 = k_t - 0.25*(*C);
	for(i = 0; i <N; ++i)
		for(j = 0; j <N; ++j){
			temp = a[i +j*N]/h2;
			diag_y_xSweep[i +j*N] =  const_var1 - temp;		
			diag_x_xSweep[i +j*N] =  const_var2 + temp;
		}
	return;
}

void d_calc_aMatrix_cMatrix(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j;
	
	for(i = 0; i <N; ++i)
		for(j = 0; j <N; ++j){
			diag_y_xSweep[i +j*N] =  k_t + 0.25*C[i +j*N] - a[i+j*N]/h2;		
			diag_x_xSweep[i +j*N] =  k_t - 0.25*C[i +j*N] + a[i+j*N]/h2;
		}
	return;
}


void d_calc_aScalar_cScalar(double *a, double *C, double h2, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j;
	double const_var1, const_var2;

	const_var1 = k_t + 0.25*(*C) - *a/h2;
	const_var2 = k_t - 0.25*(*C) + *a/h2;
	for(i = 0; i <N; ++i)
		for(j = 0; j <N; ++j){
			diag_y_xSweep[i +j*N] =  const_var1;		
			diag_x_xSweep[i +j*N] =  const_var2;
		}
	return;
}