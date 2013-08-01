#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

#ifndef DIAGS_CALC_H
#define DIAGS_CALC_H

typedef enum {IS_SCALAR, IS_MATRIX} varType;


varType find_var_type(const mxArray *mAPtr);

void diags_calculation(mxArray *plhs[], const mxArray *prhs[], double *diag_y_xSweep,
	double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep,
	const int N);
void initVariables(const mxArray *prhs[], double **a, double **C, double *k_t);


/* hyp and sub diagonals calculation*/	
void hyp_sub_diag_aScalar(mxArray *plhs[], const double* a, const int N);
void hyp_sub_diag_aMatrix(mxArray *plhs[], const double* a, const int N);
void hyp_sub_scalar_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep);
void hyp_sub_Matrix_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep, 
		const int N);



/* main diagonals calculation*/
void d_calc_aMatrix_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N);
void d_calc_aScalar_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N);
void d_calc_aMatrix_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N);
void d_calc_aScalar_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N);
		
		
/* if arguments == 8*/
void diags_update(mxArray *plhs[], const mxArray *prhs[], int N);
void diags_additions(double *d_x_add, double *d_subhyp_x_add, double *d_y_add, 
	double *d_subhyp_y_add, double *Cg, double* g, double* Cphi,
	double* phi, const int N, const mxArray *prhs[]);
void update_all_diags(mxArray *plhs[], const double* d_x_add, const 
	double* d_y_add, const int N);
void update_all_hypsub_diags(mxArray *plhs[], const mxArray *prhs[], 
		const double *d_subhyp_x_add, const double *d_subhyp_y_add, const int N);
	
	
void update_all_hypsub_diags_aScalar(mxArray *plhs[], double y_sub_x_scalar,
	double y_hyp_x_scalar,	double x_sub_x_scalar, double x_hyp_x_scalar,
	const double *d_subhyp_x_add, const double *d_subhyp_y_add, const int N);

void mx_mem_allocate(mxArray **mx_y_hyp_x, mxArray ** mx_y_sub_x, mxArray ** mx_x_hyp_x,
		mxArray ** mx_x_sub_x, const int N);

void update_all_hypsub_diags_aMatrix(double *y_sub_x, double *y_hyp_x, double *x_sub_x,
	double *x_hyp_x, const double *d_subhyp_x_add, const double *d_subhyp_y_add,
	const int N);


		
void diags_additions_xy(double *d_x_add, double *d_subhyp_x_add,
	double *d_y_add, double *d_subhyp_y_add, double *Cg, double* g, double* Cphi,
	double* phi, const int N, const mxArray *prhs[]);
void init_var(const mxArray *prhs[], double **Cg, double **g, double **Cphi,
	double **phi);
	
void restore_diags( mxArray *plhs[], double **diag_x_x, double ** diag_y_x,
	double** diag_x_y, double** diag_y_y);
	void restore_subhyp_diags(mxArray *plhs[], double **y_hyp_x, double **y_sub_x,
	double **x_hyp_x, double **x_sub_x);
void mem_alloc(double **d_x_add, double **d_y_add, double** d_subhyp_x_add,
 double** d_subhyp_y_add, const int N);
void mem_free(double *d_x_add, double *d_y_add, double* d_subhyp_x_add,
	double* d_subhyp_y_add);
		
		
/* first and second order derivatives*/
double *derivative_x(const double *a, const N);
double *derivative_xx(const double *a, const N);
double *derivative_y(const double *a, const N);
double *derivative_yy(const double *a, const N);

		
#endif