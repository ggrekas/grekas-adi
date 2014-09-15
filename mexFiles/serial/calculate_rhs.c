#include"calculate_rhs.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *rhs, *u, *sub_diag, *diag, *hyp_diag, *f, *b;
    double *bound_start_x, *bound_end_x;
    int N, i;
    
    if( 0 == nlhs || nrhs < 7)
        mexErrMsgTxt("(calculate_rhs.c): not enough input arguments");
    plhs[0]= prhs[0];
    
    rhs= mxGetPr(plhs[0]);
    u= mxGetPr(prhs[1]);
    sub_diag= mxGetPr(prhs[2]);
    diag= mxGetPr(prhs[3]);
    hyp_diag= mxGetPr(prhs[4]);
    f= mxGetPr(prhs[5]);
    
    N= mxGetM(prhs[1]);
    bound_start_x = x_Neumann_conditions_start(prhs, N); 
    bound_end_x = x_Neumann_conditions_end(prhs, N);  
/*     if( mxGetM(prhs[2])> 1){*/
         if( 1==mxGetM(prhs[2]) ){/*then sub_diag and hyp_diag are scalar values*/
            if ( 1!=mxGetM(prhs[4]) )
                mexErrMsgTxt("sub-diagonal and hyp-diagonal must have the same size");
            if( N == mxGetM(prhs[3]) ){
                boundaryStart_scalar_diffusion(u, diag, f, *sub_diag, N, 0, rhs,
                prhs, bound_start_x[0], bound_end_x[0]);
                boundaryEnd_scalar_diffusion(u, diag, f, *sub_diag, N, N-1, rhs,
                prhs, bound_start_x[N-1], bound_end_x[N-1]);
                for(i =1; i <N-1; ++i)
                    rhs_calc_scalar_diffusion(u, diag, f, *sub_diag, N, i, rhs,
                    bound_start_x[i], bound_end_x[i]);
            }
            #if 0
            else{/*diag is a scalar*/
                if( 1!=mxGetM(prhs[3]) )
                    mexErrMsgTxt("the coefficient of u (C) must be scalar");
                /*fix it................*/
            }
			#endif
         }
         else{ /*case of a:[NxN] and C:[NxN].....*/
            boundaryStart(u, sub_diag,diag, hyp_diag, f, N, 0, rhs, prhs,
            	bound_start_x[0], bound_end_x[0]);
            boundaryEnd(u, sub_diag, diag, hyp_diag, f, N, N-1, rhs, prhs,
            	bound_start_x[N-1], bound_end_x[N-1]);
            for(i =1; i <N-1; ++i)
                rhs_calc(u, sub_diag, diag, hyp_diag, f, N, i, rhs,
                 	bound_start_x[i], bound_end_x[i]);     
         }
    free(bound_start_x);
    free(bound_end_x);
    return;   
}

void boundaryStart_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs, const mxArray *prhs[],
        double bound_start_x, double bound_end_x)
    {
    int i;
    double *bound;
    
    bound = y_Neumann_conditions_start(prhs, N);
    for(i=0; i<N; ++i)
        rhs[i+ j*N]= 2*b*u[ i+(j+1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N] + bound[i];
    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;
    
    free(bound);
    return;
}


void boundaryEnd_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs, const mxArray *prhs[],
        double bound_start_x, double bound_end_x){
    int i;
    double *bound;
    
    bound = y_Neumann_conditions_end(prhs, N);    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= 2*b*u[ i+(j-1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N] + bound[i];
    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;
    
    free(bound);
    return;
}


void rhs_calc_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs, double bound_start_x, double bound_end_x){
    int i;
    for(i=0; i<N; ++i)
        rhs[i + j*N]= b*u[ i+(j-1)*N ]+ b*u[ i+(j+1)*N ] +
                d[i+j*N]*u[i+j*N]+f[i+j*N];
 
    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;           
                
   return;
}



void boundaryStart(const double *u, const double *sub_d, const double *d,
    const double *hyp_d, const double *f, int N, int j, double *rhs,const mxArray *prhs[],
    double bound_start_x, double bound_end_x)
    {
    int i;
    double *bound;
    
    bound = y_Neumann_conditions_start(prhs, N);
    for(i=0; i<N; ++i)
        rhs[i+ j*N]= (sub_d[i+j*N] + hyp_d[i+j*N])*u[ i+(j+1)*N ]
         + d[i+j*N]*u[i+j*N]+f[i+j*N] + bound[i];

    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;         
    
    free(bound);
	return;   
}
void boundaryEnd(const double *u, const double *sub_d, const double *d,
	const double *hyp_d, const double *f, int N, int j, double *rhs,
	const mxArray *prhs[], double bound_start_x, double bound_end_x)
    {
    int i;
	double *bound;
    
    bound = y_Neumann_conditions_end(prhs, N);
    for(i=0; i<N; ++i)
        rhs[i + j*N]= (hyp_d[i+j*N]+sub_d[i+j*N])*u[ i+(j-1)*N ] +
         d[i+j*N]*u[i+j*N]+f[i+j*N] + bound[i];
    
    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;

	free(bound);
    return;
}

void rhs_calc(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs,
        double bound_start_x, double bound_end_x){
    int i;
    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= sub_d[i+j*N]*u[ i+(j-1)*N ]+ hyp_d[i+j*N]*u[ i+(j+1)*N ] +
                d[i+j*N]*u[i+j*N]+f[i+j*N];
    
    rhs[0+ j*N] += bound_start_x;
    rhs[N-1+ j*N] += bound_end_x;

    return;
}

double *y_Neumann_conditions_start(const mxArray *prhs[], const int N){
    double *a, *bound, *u_y, *a_y;
    int *indxPr, i;
    
    
    indxPr = (int *)mxGetData(prhs[8]);    
    a= mxGetPr(prhs[6]);
    bound = (double *)malloc(N*sizeof(double));
    u_y = mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[0]-1) );  
        
    if( mxGetM(prhs[6]) == 1){
    	for(i=0; i<N; ++i)
    		bound[i] = -(*a) *u_y[i]*(N-1);
    }
    else{
    	a_y = mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[1]-1) );
    	for(i=0; i<N; ++i)
    		bound[i] =0.5*a_y[i]*u_y[i] -(N-1)*a[i] *u_y[i];    
    } 
	return bound;
}

double *y_Neumann_conditions_end(const mxArray *prhs[], const int N){
    double *a, *bound, *u_y, *a_y;
    int *indxPr, i;
    
    
    indxPr = (int *)mxGetData(prhs[8]);    
    a= mxGetPr(prhs[6]);
    bound = (double *)malloc(N*sizeof(double));
    u_y = N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[0]-1) );  
        
    if( mxGetM(prhs[6]) == 1){
    	for(i=0; i<N; ++i)
    		bound[i] = (*a) *u_y[i]*(N-1);
    }
    else{
    	a_y = N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[1]-1) );
    	for(i=0; i<N; ++i)
    		bound[i] =0.5*a_y[i]*u_y[i] +(N-1)*a[i+N*(N-1)] *u_y[i];    
    } 
	return bound;
}

double *x_Neumann_conditions_start(const mxArray *prhs[], const int N){
    double *a, *bound, *u_x, *a_x;
    int *indxPr, i;
    
    
    indxPr = (int *)mxGetData(prhs[8]);    
    a= mxGetPr(prhs[6]);
    bound = (double *)malloc(N*sizeof(double));
    u_x = 2*N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[0]-1) );
      
    if( mxGetM(prhs[6]) == 1){
    	for(i=0; i<N; ++i)
    		bound[i] = -(*a) *u_x[i]*(N-1);
    }
    else{
    	a_x = 2*N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[1]-1) );
    	for(i=0; i<N; ++i)
    		bound[i] = 0.5*a_x[i]*u_x[i] -(N-1)*a[i*N] *u_x[i];    
    } 
	return bound;
}

double *x_Neumann_conditions_end(const mxArray *prhs[], const int N){
    double *a, *bound, *u_x, *a_x;
    int *indxPr, i;
    
    
    indxPr = (int *)mxGetData(prhs[8]);    
    a= mxGetPr(prhs[6]);
    bound = (double *)malloc(N*sizeof(double));
    u_x = 3*N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[0]-1) );  
           
    if( mxGetM(prhs[6]) == 1){
    	for(i=0; i<N; ++i)
    		bound[i] = (*a) *u_x[i]*(N-1);
    }
    else{
    	a_x = 3*N + mxGetPr( mxGetFieldByNumber(prhs[7], 0, indxPr[1]-1) );
    	for(i=0; i<N; ++i)
    		bound[i] =0.5*a_x[i]*u_x[i] +(N-1)*a[i*N +N-1] *u_x[i];    
    } 
	return bound;
}
