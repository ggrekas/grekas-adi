mex serial/calculate_rhs.c
mex serial/TDMAsolver.c
mex serial/my_minusTranspose.c
mex serial/myTranspose.c
mex serial/diags_calc.c

if isunix
	mex paraller/calculate_rhs.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/TDMAsolver.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/my_minusTranspose.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/myTranspose.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex paraller/diags_calc.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
	mex paraller/calculate_rhs.c
	mex paraller/TDMAsolver.c
	mex paraller/my_minusTranspose.c
	mex paraller/myTranspose.c
	mex paraller/diags_calc.c
end
