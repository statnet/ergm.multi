#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

static inline double sum(double *x, unsigned int n){
  double s = 0;
  double *end = x+n;
  for(; x!=end; x++) s+=*x;
  return s;
}

static inline double sumdev2(double *x, unsigned int n){
  double s2 = 0;
  double *end = x+n;
  double mean = sum(x,n)/n;
  for(; x!=end; x++){
    double d = (*x)-mean;
    s2+=d*d;
  }
  return s2;
}


/*
  For x a concatenation of m vectors of length n, calculate the mean sample variance of those vectors.
*/
static inline double mean_var(double *x, unsigned int n, unsigned int m){
  double s2 = 0;
  double *end = x + n*m;
  for(; x!=end; x+=n){
    s2 += sumdev2(x,n);
  }
  return s2/((n-1)*m);
}

SEXP mean_var_wrapper(SEXP xe, SEXP ne){
  xe = PROTECT(coerceVector(xe, REALSXP));
  ne = PROTECT(coerceVector(ne, INTSXP));
  
  SEXP oe = PROTECT(allocVector(REALSXP, 1));
  REAL(oe)[0] = mean_var(REAL(xe), INTEGER(ne)[0], length(xe)/INTEGER(ne)[0]);
  UNPROTECT(3);
  return oe;
}

/*
  For x a concatenation of m vectors of length n, calculate the sample
  variances of those vectors. (This can be used to compute column
  variances of matrices.)
*/

SEXP vars_wrapper(SEXP xe, SEXP ne){
  xe = PROTECT(coerceVector(xe, REALSXP));
  ne = PROTECT(coerceVector(ne, INTSXP));
  unsigned int n = asInteger(ne), m = length(xe)/n;
  
  SEXP oe = PROTECT(allocVector(REALSXP, m));
  double *o = REAL(oe), *x=REAL(xe);
  for(unsigned int i=0; i<m; i++, x+=n, o++)
    *o = sumdev2(x,n)/(n-1);
  UNPROTECT(3);
  return oe;
}
