#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP splicingind(SEXP a, SEXP b, SEXP c)
{
	
	int j, na, nb, nab;
	float con,si;
	
    double *xa, *xb, *xab;
    SEXP ab;
	
    PROTECT(a = AS_NUMERIC(a));
    PROTECT(b = AS_NUMERIC(b));
	
	PROTECT(c = AS_NUMERIC(c));
    con = NUMERIC_POINTER(c)[0];
	
	na = LENGTH(a);
	nb = LENGTH(b);
	nab = na;

	PROTECT(ab = NEW_NUMERIC(nab));
	xa = NUMERIC_POINTER(a);
	xb = NUMERIC_POINTER(b);
    xab = NUMERIC_POINTER(ab);
	
	for (j = 0; j < na; j++)
	{
		if (xa[j]==0 && xb[j]==0) xab[j]=0;
		else 
			if (xa[j]==0 && xb[j]!=0) xab[j]=1;
		    else
				if (xa[j]!=0 && xb[j]==0) xab[j]=-1;
				else 
				{ 
					si=xa[j]/(xb[j]*con);
					if (si>2) xab[j]=1;
				      else 
						  if (si<0.5) xab[j]=-1;
		      			  else { xab[j]=log2(si);}
					}
		}
			

	UNPROTECT(4);
    return(ab);
}

