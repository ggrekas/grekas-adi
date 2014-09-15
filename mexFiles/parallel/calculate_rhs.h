#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

#ifndef CALCULATE_RHS_H
#define CALCULATE_RHS_H

void boundaryStart_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs, const mxArray *prhs[],
        double bound_start_x, double bound_end_x);
void boundaryEnd_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs, const mxArray *prhs[],
        double bound_start_x, double bound_end_x);
void rhs_calc_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs,double bound_start_x, double bound_end_x);

void boundaryStart(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs,
        const mxArray *prhs[], double bound_start_x, double bound_end_x);
void boundaryEnd(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs,
        const mxArray *prhs[], double bound_start_x, double bound_end_x);
void rhs_calc(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs,
        double bound_start_x, double bound_end_x);


/* new added ----------------------------------------------  */
double *y_Neumann_conditions_start(const mxArray *prhs[], const int N);
double *y_Neumann_conditions_end(const mxArray *prhs[], const int N);

double *x_Neumann_conditions_start(const mxArray *prhs[], const int N);
double *x_Neumann_conditions_end(const mxArray *prhs[], const int N);
/* new added ----------------------------------------------  */


#endif
