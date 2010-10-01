#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP gcoverage(SEXP a, SEXP b,SEXP c)
{
    int  start, i, j, na, nb, nc, nwyn;
    int	*xa, *xb, *xc, *xwyn;
    SEXP wyn;
	
	na = Rf_length(a); 
	nb = Rf_length(b);
	nc = Rf_length(c);
	nwyn = Rf_length(a);
    
	
    xa = INTEGER(a); 
    xb = INTEGER(b);
    xc = INTEGER(c);
     
	xwyn = (int *) R_alloc(nwyn, sizeof(int));
	
	i=j=0;
	
	start = xa[i];
	
	for (i=0;i<na;i++) xwyn[i]=0.0;
	
	for (i=0;i<nb;i++) 
		for(j=xb[i]-start; j <= xc[i]-start; j++)  xwyn[j]++;
	
	wyn = Rf_allocVector(INTSXP, nwyn);
        memcpy(INTEGER(wyn), xwyn, sizeof(int) * nwyn);	

    return(wyn);
	
}
  

