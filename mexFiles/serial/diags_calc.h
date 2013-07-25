#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

#ifndef DIAGS_CALC_H
#define DIAGS_CALC_H

typedef enum {IS_SCALAR, IS_MATRIX} varType;


varType find_var_type(const mxArray *mAPtr);

void diags_calculation(mxArray *plhs[], const mxArray *prhs[], double *diag_y_xSweep,
	double *diag_x_xSweep, const int N);
void initVariables(const mxArray *prhs[], double **a, double **C, double *k_t);


/* hyp and sub diagonals calculation*/	
void hyp_sub_diag_aScalar(mxArray *plhs[], const double* a, const int N);
void hyp_sub_diag_aMatrix(mxArray *plhs[], const double* a, const int N);
void hyp_sub_scalar_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep);
void hyp_sub_Matrix_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep, 
		const int N);

double *derivative_x(const double *a, const N);
double *derivative_y(const double *a, const N);

/* main diagonals calculation*/
void d_calc_aMatrix_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aScalar_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aMatrix_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
void d_calc_aScalar_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N);
		
		
		
#endif