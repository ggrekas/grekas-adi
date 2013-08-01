#include"diags_calc.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *diag_y_xSweep, *diag_x_xSweep, *diag_y_ySweep, *diag_x_ySweep;
	int N;

	/*
	if(mxGetM(prhs[0]) < mxGetM(prhs[1]) )
		N= (int)mxGetM(prhs[1]);
	else
		N= (int)mxGetM(prhs[0]);
	*/	
	N = (int)mxGetScalar(prhs[3]);
	
	/*plhs[1] = plhs[2]; 
	plhs[0] = plhs[3];*/
	omp_set_num_threads(18);
	
	plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(N, N, mxREAL);


	diag_y_xSweep = mxGetPr(plhs[0]);
	diag_x_xSweep = mxGetPr(plhs[1]);
	diag_y_ySweep = mxGetPr(plhs[2]);
	diag_x_ySweep = mxGetPr(plhs[3]);
	
	
	diags_calculation(plhs, prhs, diag_y_xSweep, diag_x_xSweep, diag_y_ySweep,
	 	diag_x_ySweep, N);
	if(8 == nrhs){ 
	/*	plhs[2] = mxCreateDoubleMatrix(N, N, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(N, N, mxREAL);
*/
		diags_update(plhs, prhs, N); 
	}
	
	return;
}


void diags_update(mxArray *plhs[], const mxArray *prhs[], int N)
{
	double *d_x_add, *d_y_add, *d_subhyp_x_add, *d_subhyp_y_add;
	double *Cg, *g, *Cphi, *phi;
	
	mem_alloc(&d_x_add, &d_y_add, &d_subhyp_x_add, &d_subhyp_y_add, N);
	init_var(prhs, &Cg, &g, &Cphi, &phi);
	diags_additions(d_x_add, d_subhyp_x_add, d_y_add, d_subhyp_y_add, Cg, g,
		Cphi, phi, N, prhs);		
	update_all_diags(plhs, d_x_add, d_y_add, N);
	update_all_hypsub_diags(plhs, prhs, d_subhyp_x_add, d_subhyp_y_add, N);

	mem_free(d_x_add, d_y_add, d_subhyp_x_add, d_subhyp_y_add);
	
	return;
}

void update_all_diags(mxArray *plhs[], const double* d_x_add, const 
	double* d_y_add, const int N)
{
	int i, j, k;
	double *diag_x_x, *diag_y_x, *diag_x_y, *diag_y_y; 
	
	restore_diags( plhs, &diag_x_x, &diag_y_x, &diag_x_y, &diag_y_y);
	#pragma omp parallel shared(diag_y_y, diag_x_y, diag_y_x, diag_x_x, d_y_add,\
		d_x_add)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j)
		for(i = 0; i < N; ++i){
			k = i + j*N;
			diag_y_y[k] = diag_x_x[k] - 0.5 * d_y_add[k];
			diag_x_y[k] = diag_y_x[k] + 0.5 * d_x_add[k];
			
			diag_y_x[k] += 0.5 * d_y_add[k];
			diag_x_x[k] -= 0.5 * d_x_add[k];
			
		}
	return;
}

void restore_diags( mxArray *plhs[], double **diag_x_x, double ** diag_y_x,
	double** diag_x_y, double** diag_y_y)
{
	*diag_y_x = mxGetPr(plhs[0]);
	*diag_x_x = mxGetPr(plhs[1]);
	*diag_y_y = mxGetPr(plhs[2]);
	*diag_x_y = mxGetPr(plhs[3]);
		
	return;
}

void update_all_hypsub_diags(mxArray *plhs[], const mxArray *prhs[], 
		const double *d_subhyp_x_add, const double *d_subhyp_y_add, const int N)
{
	double *y_sub_x, *y_hyp_x, *x_sub_x, *x_hyp_x;

	restore_subhyp_diags( plhs, &y_hyp_x, &y_sub_x, &x_hyp_x, &x_sub_x);
	if( N == mxGetM(prhs[0]) )
		update_all_hypsub_diags_aMatrix(y_sub_x, y_hyp_x, x_sub_x, x_hyp_x,
			d_subhyp_x_add, d_subhyp_y_add, N);
	else	
		update_all_hypsub_diags_aScalar(plhs, *y_sub_x, *y_hyp_x, *x_sub_x, *x_hyp_x,
			d_subhyp_x_add, d_subhyp_y_add, N);
			
	return;
}

void update_all_hypsub_diags_aScalar(mxArray *plhs[], double y_sub_x_scalar,
	double y_hyp_x_scalar,	double x_sub_x_scalar, double x_hyp_x_scalar,
	const double *d_subhyp_x_add, const double *d_subhyp_y_add, const int N)
{
	int i, j, N2, k;
	mxArray *mx_y_hyp_x, *mx_y_sub_x, *mx_x_hyp_x, *mx_x_sub_x;  
	double *y_sub_x, *y_hyp_x, *x_sub_x, *x_hyp_x;
	
	mx_mem_allocate(&mx_y_hyp_x, &mx_y_sub_x, &mx_x_hyp_x, &mx_x_sub_x, N);	
	
	y_hyp_x = mxGetPr(mx_y_hyp_x);
	y_sub_x = mxGetPr(mx_y_sub_x);
	x_hyp_x = mxGetPr(mx_x_hyp_x);
	x_sub_x = mxGetPr(mx_x_sub_x);
	
	
	N2 = (N-1)*(N-1);
	#pragma omp parallel shared(y_hyp_x, y_sub_x, y_hyp_x_scalar, y_sub_x_scalar\
		,d_subhyp_y_add, x_hyp_x, x_sub_x, x_hyp_x_scalar, x_sub_x_scalar,\
		d_subhyp_x_add, N2)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j)
		for(i = 0; i < N; ++i){
			k = i + j*N;
			y_hyp_x[k] = y_hyp_x_scalar + 0.125 * N2 * d_subhyp_y_add[k];
			y_sub_x[k] = y_sub_x_scalar - 0.125 * N2 * d_subhyp_y_add[k];
		
			x_hyp_x[k] = x_hyp_x_scalar - 0.125 * N2 * d_subhyp_x_add[k];
			x_sub_x[k] = x_sub_x_scalar + 0.125 * N2 * d_subhyp_x_add[k];
		}
		
	mxDestroyArray(plhs[4]);
	mxDestroyArray(plhs[5]);
	mxDestroyArray(plhs[6]);
	mxDestroyArray(plhs[7]);
	
	plhs[4] = mx_y_hyp_x;
	plhs[5] = mx_y_sub_x;
	plhs[6] = mx_x_hyp_x;
	plhs[7] = mx_x_sub_x;
	
	return;
	
}

void mx_mem_allocate(mxArray **mx_y_hyp_x, mxArray ** mx_y_sub_x, mxArray ** mx_x_hyp_x,
		mxArray ** mx_x_sub_x, const int N)
{
	*mx_y_hyp_x = mxCreateDoubleMatrix(N, N, mxREAL);
	*mx_y_sub_x = mxCreateDoubleMatrix(N, N, mxREAL);
	*mx_x_hyp_x = mxCreateDoubleMatrix(N, N, mxREAL);
	*mx_x_sub_x = mxCreateDoubleMatrix(N, N, mxREAL);

	return;
}


void update_all_hypsub_diags_aMatrix(double *y_sub_x, double *y_hyp_x, double *x_sub_x,
	double *x_hyp_x, const double *d_subhyp_x_add, const double *d_subhyp_y_add,
	const int N)
{
	int i, j, k;
	double N2;
	N2 = 0.125*(N-1)*(N-1);
	
	#pragma omp parallel shared(y_hyp_x, y_sub_x, d_subhyp_y_add, x_hyp_x,\
		x_sub_x, d_subhyp_x_add, N2)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j)
		for(i = 0; i < N; ++i){
			k = i + j*N;
			y_hyp_x[k] += N2 * d_subhyp_y_add[k];
			y_sub_x[k] -= N2 * d_subhyp_y_add[k];
		
			x_hyp_x[k] -= N2 * d_subhyp_x_add[k];
			x_sub_x[k] += N2 * d_subhyp_x_add[k];
		}
	return;
	
}


void restore_subhyp_diags(mxArray *plhs[], double **y_hyp_x, double **y_sub_x,
	double **x_hyp_x, double **x_sub_x)
{
	*y_hyp_x = mxGetPr(plhs[4]);
	*y_sub_x = mxGetPr(plhs[5]);
	*x_hyp_x = mxGetPr(plhs[6]);
	*x_sub_x = mxGetPr(plhs[7]);

	return;
} 


void diags_additions(double *d_x_add, double *d_subhyp_x_add, double *d_y_add, 
	double *d_subhyp_y_add, double *Cg, double* g, double* Cphi,
	double* phi, const int N, const mxArray *prhs[])
{
	int i, j, k;
	double *g_x, *g_xx, *phi_x, *phi_xx, *g_y, *g_yy, *phi_y, *phi_yy;
	
	g_x = derivative_x(g, N);
	g_xx = derivative_xx(g, N);
	phi_x = derivative_x(phi, N);
	phi_xx = derivative_xx(phi, N);
	g_y = derivative_y(g, N);
	g_yy = derivative_yy(g, N);
	phi_y = derivative_y(phi, N);
	phi_yy = derivative_yy(phi, N);
	
	
	if(1 == mxGetM(prhs[4]) )
	{
	#pragma omp parallel shared(d_subhyp_x_add, d_x_add, d_subhyp_y_add, d_y_add\
		,Cg, g_x, g_xx, g_y, g_yy)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
		for(j = 0; j < N; ++j)
			for(i = 0; i < N; ++i){
				k = i + j*N;
				d_subhyp_x_add[k] = *Cg *g_x[k];
				d_x_add[k] = *Cg *g_xx[k];
				d_subhyp_y_add[k] = *Cg *g_y[k];
				d_y_add[k] = *Cg *g_yy[k];
			}
	}
	else
	{
	#pragma omp parallel shared(d_subhyp_x_add, d_x_add, d_subhyp_y_add, d_y_add\
		,Cg, g_x, g_xx, g_y, g_yy)\
		private(i, j, k)
	#pragma omp for schedule(guided) nowait
		for(j = 0; j < N; ++j)
			for(i = 0; i < N; ++i){
				k = i + j*N;
				d_subhyp_x_add[k] = Cg[k] *g_x[k];
				d_x_add[k] = Cg[k] *g_xx[k];
				d_subhyp_y_add[k] = Cg[k] *g_y[k];
				d_y_add[k] = Cg[k] *g_yy[k];
			}
	}
	if(1 == mxGetM(prhs[6]) ){
		#pragma omp parallel shared(d_subhyp_x_add, d_x_add, d_subhyp_y_add, d_y_add\
			,Cphi, phi_x, phi_xx, phi_y, phi_yy)\
		private(i, j, k)
		#pragma omp for schedule(guided) nowait
		for(j = 0; j < N; ++j)
			for(i = 0; i < N; ++i){
				k = i + j*N;
				d_subhyp_x_add[k] += *Cphi *phi_x[k];
				d_x_add[k] += *Cphi *phi_xx[k];
				d_subhyp_y_add[k] += *Cphi *phi_y[k];
				d_y_add[k] += *Cphi *phi_yy[k];
			}
	}
	else
	{
		#pragma omp parallel shared(d_subhyp_x_add, d_x_add, d_subhyp_y_add, d_y_add\
			,Cphi, phi_x, phi_xx, phi_y, phi_yy)\
		private(i, j, k)
		#pragma omp for schedule(guided) nowait
		for(j = 0; j < N; ++j)
			for(i = 0; i < N; ++i){
				k = i + j*N;
				d_subhyp_x_add[k] += Cphi[k] *phi_x[k];
				d_x_add[k] += Cphi[k] *phi_xx[k];
				d_subhyp_y_add[k] += Cphi[k] *phi_y[k];
				d_y_add[k] += Cphi[k] *phi_yy[k];
			}	
	}
	
	
	
	free(g_x);
	free(g_xx);
	free(phi_x);
	free(phi_xx);
	free(g_y);
	free(g_yy);
	free(phi_y);
	free(phi_yy);

	return;
}	



void init_var(const mxArray *prhs[], double **Cg, double **g, double **Cphi,
	double **phi)
{
	*Cg = mxGetPr(prhs[4]);
	*g = mxGetPr(prhs[5]);
	*Cphi = mxGetPr(prhs[6]);
	*phi = mxGetPr(prhs[7]);
	
	return;
}


void mem_alloc(double **d_x_add, double **d_y_add, double** d_subhyp_x_add,
 double** d_subhyp_y_add, const int N){
	
	*d_x_add = (double *)malloc(N*N*sizeof(double));
	*d_y_add = (double *)malloc(N*N*sizeof(double));
	*d_subhyp_y_add = (double *)malloc(N*N*sizeof(double));
	*d_subhyp_x_add = (double *)malloc(N*N*sizeof(double));
 
	return;
}


void mem_free(double *d_x_add, double *d_y_add, double* d_subhyp_x_add,
	double* d_subhyp_y_add)
{
	free(d_x_add);
	free(d_y_add);
	free(d_subhyp_y_add);
	free(d_subhyp_x_add);
	
	return; 
}


void diags_calculation(mxArray *plhs[], const mxArray *prhs[], double *diag_y_xSweep,	double *diag_x_xSweep, double *diag_y_ySweep,	double *diag_x_ySweep, const int N)
	{
	double *a, *C, k_t;
	void (*fPtr[2][2])(double*, double*, double, double *, double *, double *,
		double *, int ) = {{d_calc_aScalar_cScalar, d_calc_aScalar_cMatrix},
		{d_calc_aMatrix_cScalar, d_calc_aMatrix_cMatrix}};
	void (*fSubHypDiagPtr[2])(mxArray **, const double*, const int ) =
				{hyp_sub_diag_aScalar,	hyp_sub_diag_aMatrix};
	
	
	initVariables( prhs, &a, &C, &k_t);
	fPtr[find_var_type(prhs[0])][find_var_type(prhs[1])](a, C, k_t, diag_y_xSweep,
		diag_x_xSweep, diag_y_ySweep, diag_x_ySweep, N);
	
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
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	 *y_hypDiag_xSweep = mxGetPr(plhs[4]);
	 *y_subDiag_xSweep = mxGetPr(plhs[5]);
	 *x_hypDiag_xSweep = mxGetPr(plhs[6]);
	 *x_subDiag_xSweep = mxGetPr(plhs[7]);	

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



void hyp_sub_Matrix_mem_alloc(mxArray *plhs[], double **y_hypDiag_xSweep, double**
		y_subDiag_xSweep, double** x_hypDiag_xSweep, double** x_subDiag_xSweep,
		const int N)
{
	plhs[4] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(N, N, mxREAL);
	plhs[7] = mxCreateDoubleMatrix(N, N, mxREAL);
	
	 *y_hypDiag_xSweep = mxGetPr(plhs[4]);
	 *y_subDiag_xSweep = mxGetPr(plhs[5]);
	 *x_hypDiag_xSweep = mxGetPr(plhs[6]);
	 *x_subDiag_xSweep = mxGetPr(plhs[7]);	

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
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N)
{
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
			diag_x_ySweep[k] = diag_y_xSweep[k] =  const_var1 + temp;
			diag_y_ySweep[k] = diag_x_xSweep[k] =  const_var2 - temp;
		}
	return;
}

void d_calc_aMatrix_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N)
{
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
			diag_x_ySweep[k] = diag_y_xSweep[k] =  const_var1 - temp;		
			diag_y_ySweep[k] =diag_x_xSweep[k] =  const_var2 + temp;
		}
	return;
}

void d_calc_aMatrix_cMatrix(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N)
{
	int i, j, k;
	double N2;
	
	N2 = (N-1)*(N-1);
	#pragma omp parallel shared(diag_y_xSweep, diag_x_xSweep, k_t, C, a, N2, N)\
		private(i, j, k)
		#pragma omp for schedule(guided) nowait
		for(j = 0; j <N; ++j)
			for(i = 0; i <N; ++i){
				k = i +j*N;
				diag_x_ySweep[k] = diag_y_xSweep[k] = k_t + 0.25*C[k] - a[k]*N2;
				diag_y_ySweep[k] = diag_x_xSweep[k] = k_t - 0.25*C[k] + a[k]*N2;
			}
	
	return;
}


void d_calc_aScalar_cScalar(double *a, double *C, double k_t, double *diag_y_xSweep,
		double *diag_x_xSweep, double *diag_y_ySweep, double *diag_x_ySweep, int N)
{
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
			diag_x_ySweep[k] = diag_y_xSweep[k] =  const_var1;		
			diag_y_ySweep[k] = diag_x_xSweep[k] =  const_var2;
		}
	return;
}




double *derivative_x(const double *a, const N){
	int i, j, k;
	double *a_x;
	
	a_x = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_x, a)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j){
		a_x[j*N] =0;		
		for(i = 1; i < N-1; ++i){
			k = i+j*N;
			a_x[k] = a[k+1] - a[k-1];
		}
		a_x[N-1 + j*N] =0;	
	}
				
	return a_x;
}


double *derivative_xx(const double *a, const N){
	int i, j, k;
	double *a_xx;
	
	a_xx = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_xx, a)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 0; j < N; ++j){
		a_xx[j*N] =2*(a[1 + j*N] - a[j*N]);		
		for(i = 1; i < N-1; ++i){
			k = i+j*N;
			a_xx[k] = a[k+1] + a[k-1] - 2*a[k];		
		}
		a_xx[N-1 + j*N] =2*(a[N-2 + j*N] - a[N-1 + j*N]);	
	}
				
	return a_xx;
}

double *derivative_y(const double *a, const N){
	int i, j, k;
	double *a_y;
	
	a_y = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_y, a)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 1; j < N-1; ++j){
		for(i = 0; i < N; ++i){
			k = i+j*N;
			a_y[k] = a[k + N] - a[k-N];
		}
	}
	
	#pragma omp parallel shared(a_y)\
	private(i, j)
	#pragma omp for schedule(guided,10) nowait
	for(j =0; j < N; j+=N-1)
		for(i = 0; i < N; ++i)
			a_y[i+j*N] = 0;
		
	return a_y;
}


double *derivative_yy(const double *a, const N){
	int i, j, k, N2;
	double *a_yy;
	
	a_yy = (double *)malloc(N*N*sizeof(double) );
	
	#pragma omp parallel shared(a_yy, a)\
	private(i, j, k)
	#pragma omp for schedule(guided) nowait
	for(j = 1; j < N-1; ++j){
		for(i = 0; i < N; ++i){
			k = i+j*N;
			a_yy[k] = a[k+N] + a[k-N] -2*a[k];
		}
	}
	
	N2 = N*N;
	#pragma omp parallel shared(a_yy, N2)\
	private(i)
	#pragma omp for schedule(guided,10) nowait
		for(i = 0; i < N; ++i){
				a_yy[i] = 2*(a[i+N] - a[i]);
				a_yy[i+N2-N] = 2*(a[i+N2-2*N] - a[i+N2-N]);
		}
	return a_yy;
}
