mex calculate_rhs.c
mex TDMAsolver.c
if isunix
	mex TDMAsolver_par.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
	mex TDMAsolver_par.c
end
mex my_minusTranspose.c
mex myTranspose.c