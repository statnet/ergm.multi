/*  File src/init.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
/* This is a dummy list, since ergm.multi is a library package.
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP mean_var_wrapper(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mean_var_wrapper",   (DL_FUNC) &mean_var_wrapper,   2},
    {NULL, NULL, 0}
};

void R_init_ergm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
