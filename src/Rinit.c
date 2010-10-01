#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>

SEXP gcoverage(SEXP a, SEXP b,SEXP c);
SEXP ghistogram(SEXP a, SEXP b,SEXP skip);
SEXP regionmining(SEXP n, SEXP nr,SEXP p);
SEXP splicingind(SEXP a, SEXP b, SEXP c);

 static const R_CallMethodDef R_CallDef[] = {
	{ "gcoverage", (DL_FUNC) &gcoverage, 3 },
	{ "ghistogram", (DL_FUNC) &ghistogram, 3 },
	{ "regionmining", (DL_FUNC) &regionmining, 4 },
	{ "splicingind", (DL_FUNC) &splicingind, 3 },
        {NULL, NULL, 0}
    };

void R_init_rnaSeqMap(DllInfo *info)
{
  R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
}

