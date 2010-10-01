#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

 double avg (int *tab,int s,int e)
{
	int i, n;
        n=1;
	double avg=0;
	for (i=s;i<=e;i++){
		avg	= (((avg *(n-1))/n+(double)tab[i]/n));
		n=n++;
	}
	return avg;
} 



SEXP regionmining(SEXP n, SEXP nr,SEXP p, SEXP ms)
{
    int a, b, i,j , last,nn, nnr, nwyn;
    int current, start, param, minsup,*xwyn; 
	int *xn, *xnr; 
	
	SEXP wyn;
	    
	param = INTEGER_POINTER(p)[0];
	minsup = INTEGER_POINTER(ms)[0];
	
	nn = Rf_length(n);
	nnr = Rf_length(nr);
	nwyn = Rf_length(n);
	
	
        xn = INTEGER(n); 
	xnr = INTEGER(nr);
	
	xwyn = (int *) R_alloc(2*nwyn,sizeof(int));
	
	
	a=b=i=j=0;
	last=nnr-1;
	current=0;
	start=xn[current];
   
	xwyn[i]=0;
	
	while (current <= last) 
	    {
			while ((xnr[current]<param) && current <last)  
				current++; 
			 if (current <= last)
			 {
				a=current;
				b=a;
				 while ((avg(xnr,a,current)>=param) && (current<last)) 
				 {
					 current++;
					 if (avg(xnr,b+1,current)>=param) //zmienione!!!
					    {
						 b=current;
					    }
				
				  }
                                    if (b-a>=minsup)
                                      {
					 xwyn[j]=xn[a];
					 j++;
				         if  (b==nn) xwyn[j]=xn[b-1];
					 else xwyn[j]=xn[b];
					 j++;
                                      }
				    current=b+1;
			 }
		 }
	if (xwyn[i]==0) j=1;
	wyn = Rf_allocVector(INTSXP,j);
	memcpy(INTEGER(wyn),xwyn,sizeof(int) * j);
    return(wyn);
}


