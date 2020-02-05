#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mean_var_wrapper(SEXP, SEXP);
extern SEXP vars_wrapper(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mean_var_wrapper", (DL_FUNC) &mean_var_wrapper, 2},
    {"vars_wrapper",     (DL_FUNC) &vars_wrapper,     2},
    {NULL, NULL, 0}
};

void R_init_ergm_multi(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
