mex serial/calculate_rhs.c
mex serial/TDMAsolver.c
mex serial/my_minusTranspose.c
mex serial/myTranspose.c
mex serial/diags_calc.c
if isunix
	mex paraller/calculate_rhs_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/TDMAsolver_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/my_minusTranspose_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/myTranspose_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/diags_calc_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
	mex paraller/calculate_rhs_par.c
	mex paraller/TDMAsolver_par.c
	mex paraller/my_minusTranspose_par.c
	mex paraller/myTranspose_par.c
	mex paraller/diags_calc_par.c
end
