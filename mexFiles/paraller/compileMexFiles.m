if isunix
	mex calculate_rhs.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex TDMAsolver.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex my_minusTranspose.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex myTranspose.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
	mex diags_calc.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
	mex calculate_rhs.c
	mex TDMAsolver.c
	mex my_minusTranspose.c
	mex myTranspose.c
	mex diags_calc.c
end
