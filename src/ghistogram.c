#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP ghistogram(SEXP a, SEXP b,SEXP skip)
{
	
	int i, j, l , na, nb, nab, Is;
    int *xa, *xb, *xab;
    SEXP ab;
     int k;
	
	Is = INTEGER_POINTER(skip)[0];
    
    na = Rf_length(a);
	nb = Rf_length(b);
	
	nab = na/Is;
	 
	xa = INTEGER(a);
	xb = INTEGER(b);
   
	xab = (int *) R_alloc(nab, sizeof(int));
	
	i=j=k=0;
	
	for(l = 0; l < nab; l++) xab[l] = 0.0;
    
	while ( i < na) {
		for(j = i; j < i+Is; j++)
			xab[k] += xb[j];
	        k=k+1;
	i=i+Is;
	}
    
	ab = Rf_allocVector(INTSXP, nab);
	memcpy(INTEGER(ab), xab, sizeof(int) * nab);	

	return(ab);
	
}


