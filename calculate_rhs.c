#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

void boundaryStart_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs);
void boundaryEnd_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs);
void rhs_calc_scalar_diffusion(const double *u, const double *d, const double *f,
        double  b, int N, int j, double *rhs);

void boundaryStart(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs);
void boundaryEnd(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs);
void rhs_calc(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *rhs, *u, *sub_diag, *diag, *hyp_diag, *f, *b;
    int N, i;
    
    if( 0 == nlhs || nrhs < 5)
        mexErrMsgTxt("not enough input arguments");
    plhs[0]= prhs[0];
    
    rhs= mxGetPr(plhs[0]);
    u= mxGetPr(prhs[1]);
    sub_diag= mxGetPr(prhs[2]);
    diag= mxGetPr(prhs[3]);
    hyp_diag= mxGetPr(prhs[4]);
    f= mxGetPr(prhs[5]);
    
    N= mxGetM(prhs[1]) ;
    
/*     if( mxGetM(prhs[2])> 1){*/
         if( 1==mxGetM(prhs[2]) ){/*then sub_diag and hyp_diag are scalar values*/
            if ( 1!=mxGetM(prhs[4]) )
                mexErrMsgTxt("sub-diagonal and hyp-diagonal must be of the same size");
            if( N ==mxGetM(prhs[3]) ){
                boundaryStart_scalar_diffusion(u, diag, f, *sub_diag, N, 0, rhs);
                boundaryEnd_scalar_diffusion(u, diag, f, *sub_diag, N, N-1, rhs);
                for(i =1; i <N-1; ++i)
                    rhs_calc_scalar_diffusion(u, diag, f, *sub_diag, N, i, rhs);
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
            boundaryStart(u, sub_diag,diag, hyp_diag, f, N, 0, rhs);
            boundaryEnd(u, sub_diag, diag, hyp_diag, f, N, N-1, rhs);
            for(i =1; i <N-1; ++i)
                rhs_calc(u, sub_diag, diag, hyp_diag, f, N, i, rhs);     
         }
    return;   
}

void boundaryStart_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs){
    int i;
    for(i=0; i<N; ++i)
        rhs[i+ j*N]= 2*b*u[ i+(j+1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N];
      
    return;
}

void boundaryEnd_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs){
    int i;
    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= 2*b*u[ i+(j-1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N];
    return;
}

void rhs_calc_scalar_diffusion(const double *u, const  double * d, const double * f,
        double  b, int N, int j, double *rhs){
    int i;
    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= b*u[ i+(j-1)*N ]+ b*u[ i+(j+1)*N ] +
                d[i+j*N]*u[i+j*N]+f[i+j*N];
	return;
}



void boundaryStart(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs){
    int i;
    for(i=0; i<N; ++i)
        rhs[i+ j*N]= sub_d[i+j*N]*u[ i+(j+1)*N ] + hyp_d[i+j*N]*u[ i+(j+1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N];
	return;   
}
void boundaryEnd(const double *u, const double *sub_d, const double *d,
         const double *hyp_d, const double *f, int N, int j, double *rhs){
    int i;
    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= hyp_d[i+j*N]*u[ i+(j-1)*N ] + sub_d[i+j*N]*u[ i+(j-1)*N ] + d[i+j*N]*u[i+j*N]+f[i+j*N];
    return;
}
void rhs_calc(const double *u, const double *sub_d, const double *d,
        const double *hyp_d, const double *f, int N, int j, double *rhs){
    int i;
    
    for(i=0; i<N; ++i)
        rhs[i + j*N]= sub_d[i+j*N]*u[ i+(j-1)*N ]+ hyp_d[i+j*N]*u[ i+(j+1)*N ] +
                d[i+j*N]*u[i+j*N]+f[i+j*N];
    return;
}



