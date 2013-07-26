#include"diags_calc.h"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *diag_y_xSweep, *diag_x_xSweep;
	int N;

	/*
	if(mxGetM(prhs[0]) < mxGetM(prhs[1]) )
		N= (int)mxGetM(prhs[1]);
	else
		N= (int)mxGetM(prhs[0]);
	*/	
	N = (int)mxGetScalar(prhs[3]);
	
	plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[2] = plhs[3] = plhs[4] = plhs[5] = NULL;
	diag_y_xSweep = mxGetPr(plhs[0]);
	diag_x_xSweep = mxGetPr(plhs[1]);
	
	
	diags_calculation(plhs, prhs, diag_y_xSweep, diag_x_xSweep, N);

	return;
}

void diags_calculation(mxArray *plhs[], const mxArray *prhs[], double *diag_y_xSweep,
	double *diag_x_xSweep, const int N)
	{
	double *a, *C, k_t;
	void (*fPtr[2][2])(double*, double*, double, double *, double *,
		int ) = {{d_calc_aScalar_cScalar, d_calc_aScalar_cMatrix},
		{d_calc_aMatrix_cScalar, d_calc_aMatrix_cMatrix}};
	void (*fSubHypDiagPtr[2])(mxArray **, const double*, const int ) =
				{hyp_sub_diag_aScalar,	hyp_sub_diag_aMatrix};
	
	
	initVariables( prhs, &a, &C, &k_t);
	fPtr[find_var_type(prhs[0])][find_var_type(prhs[1])](a, C, k_t, diag_y_xSweep,
		diag_x_xSweep, N);
	
	fSubHypDiagPtr[find_var_type(prhs[0])](plhs, a, N);
	
	
	return;
}


void hyp_sub_diag_aScalar(mxArray *plhs[], const double* a, const int N)
{
	double *y_hypDiag_xSweep, *y_subDiag_xSweep, *x_hypDiag_xSweep, *x_subDiag_xSweep;
	
	hyp_sub_scalar_mem_alloc(plhs, &y_hypDiag_xSweep, &y_subDiag_xSweep,
		&x_hypDiag_xSweep, &x_subDiag_xSweep);

	*y_hypDiag_xSweep = 0.5*(N-1)*(N-1) * *a;
	*y_subDiag_xSweep = 0.5*(N-1)*(N-1) * *a;
	*x_hypDiag_xSweep = -0.5*(N-1)*(N-1) * *a;
	*x_subDiag_xSweep = -0.5*(N-1)*(N-1) * *a;
	
	return;
}

void hyp_sub_scalar_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep)
{
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	 *y_hypDiag_xSweep = mxGetPr(plhs[2]);
	 *y_subDiag_xSweep = mxGetPr(plhs[3]);
	 *x_hypDiag_xSweep = mxGetPr(plhs[4]);
	 *x_subDiag_xSweep = mxGetPr(plhs[5]);	

	return;
}

void hyp_sub_diag_aMatrix(mxArray *plhs[], const double* a, const int N)
{
	int i, j, k=0;
	double  N2;
	double *y_hypDiag_xSweep, *y_subDiag_xSweep, *x_hypDiag_xSweep, *x_subDiag_xSweep;
	double *a_x, *a_y;
	
	hyp_sub_Matrix_mem_alloc(plhs, &y_hypDiag_xSweep, &y_subDiag_xSweep,
		&x_hypDiag_xSweep, &x_subDiag_xSweep, N);
		
	a_x = derivative_x(a, N);	
	a_y = derivative_y(a, N);
	N2 = 0.5*(N-1)*(N-1);
	
	#pragma omp parallel shared(a_y, a_x, a, N2, y_hypDiag_xSweep, y_subDiag_xSweep,\
	x_hypDiag_xSweep, x_subDiag_xSweep)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j)
		for(i = 0; i < N; ++i)
		{
			k = i +j*N;
			y_hypDiag_xSweep[k] = N2*( a[k] + 0.25*a_y[k] );
			y_subDiag_xSweep[k] = N2*( a[k] - 0.25*a_y[k] );
		
			x_hypDiag_xSweep[k] = -N2*( a[k] + 0.25*a_x[k] );
			x_subDiag_xSweep[k] = -N2*( a[k] - 0.25*a_x[k] );
		}

	free(a_x);
	free(a_y);
	return;
}

double *derivative_x(const double *a, const N){
	int i, j;
	double *a_x;
	
	a_x = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_x, a, N)\
	private(i, j)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j){
		a_x[j*N] =0;		
		for(i = 1; i < N-1; ++i)
			a_x[i+j*N] = a[(i+1) + j*N] - a[(i-1) + j*N];
		a_x[N-1 + j*N] =0;	
	}
				
	return a_x;
}

double *derivative_y(const double *a, const N){
	int i, j;
	double *a_y;
	
	a_y = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_y, a, N)\
	private(i, j)
	#pragma omp for schedule(guided) nowait
	for(j = 1; j < N-1; ++j){
		for(i = 0; i < N; ++i)
			a_y[i+j*N] = a[i+(j+1)*N] - a[i+(j-1)*N];
	}
	
	#pragma omp parallel shared(a_y, N)\
	private(i, j)
	#pragma omp for schedule(guided,10) nowait
	for(j =0; j < N; j+=N)
		for(i = 0; i < N; ++i)
			a_y[i+j*N] = 0;
		
	return a_y;
}

void hyp_sub_Matrix_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep,
		const int N)
{
	plhs[2] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(N, N, mxREAL);
	
	 *y_hypDiag_xSweep = mxGetPr(plhs[2]);
	 *y_subDiag_xSweep = mxGetPr(plhs[3]);
	 *x_hypDiag_xSweep = mxGetPr(plhs[4]);
	 *x_subDiag_xSweep = mxGetPr(plhs[5]);	

	return;
}


void initVariables(const mxArray *prhs[], double **a, double **C, double *k_t){
	
	
	*a = mxGetPr(prhs[0]);
	*C = mxGetPr(prhs[1]);
	*k_t = mxGetScalar(prhs[2]);
	
	return;
}



varType find_var_type(const mxArray *mAPtr){
	if(1 == mxGetM(mAPtr))
		return IS_SCALAR;
		
	return IS_MATRIX;
}

void d_calc_aScalar_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j, k;
	double const_var1, const_var2, temp;
	double N2;
	
	N2 = (N-1)*(N-1);
	const_var1 = k_t - *a*N2;
	const_var2 = k_t + *a*N2;
	
	#pragma omp parallel shared(diag_y_xSweep, diag_x_xSweep, C, const_var1, const_var2, N)\
	private(i, j, temp, k)
		#pragma omp for schedule(guided) nowait
	for(j = 0; j <N; ++j)
		for(i = 0; i <N; ++i){	
			k = i +j*N;
			temp = 0.25*C[k];
			diag_y_xSweep[k] =  const_var1 + temp;
			diag_x_xSweep[k] =  const_var2 - temp;
		}
	return;
}

void d_calc_aMatrix_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j, k;
	double const_var1, const_var2, temp;
	double N2;
	
	N2 = (N-1)*(N-1);
	const_var1 = k_t + 0.25*(*C);
	const_var2 = k_t - 0.25*(*C);
	
	#pragma omp parallel shared(diag_y_xSweep, diag_x_xSweep, a, const_var1, const_var2, N2, N)\
	private(i, j, temp, k)
		#pragma omp for schedule(guided) nowait
	for(j = 0; j <N; ++j)
		for(i = 0; i <N; ++i){
			k = i +j*N;
			temp = a[k]*N2;
			diag_y_xSweep[k] =  const_var1 - temp;		
			diag_x_xSweep[k] =  const_var2 + temp;
		}
	return;
}

void d_calc_aMatrix_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j, k;
	double N2;
	
	N2 = (N-1)*(N-1);
	#pragma omp parallel shared(diag_y_xSweep, diag_x_xSweep, k_t, C, a, N2, N)\
		private(i, j, k)
		#pragma omp for schedule(guided) nowait
		for(j = 0; j <N; ++j)
			for(i = 0; i <N; ++i){
				k = i +j*N;
				diag_y_xSweep[k] =  k_t + 0.25*C[k] - a[k]*N2;		
				diag_x_xSweep[k] =  k_t - 0.25*C[k] + a[k]*N2;
			}
	
	return;
}


void d_calc_aScalar_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, int N){
	int i, j, k;
	double const_var1, const_var2;
	double N2;

	N2 = (N-1)*(N-1);
	const_var1 = k_t + 0.25*(*C) - *a*N2;
	const_var2 = k_t - 0.25*(*C) + *a*N2;
	
	#pragma omp parallel shared(diag_y_xSweep, diag_x_xSweep, const_var1, const_var2)\
	private(i, j, k)
		#pragma omp for schedule(guided) nowait
	for(j = 0; j <N; ++j)
		for(i = 0; i <N; ++i){
		k = i + j*N;
			diag_y_xSweep[k] =  const_var1;		
			diag_x_xSweep[k] =  const_var2;
		}
	return;
}