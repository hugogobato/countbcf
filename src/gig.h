#ifndef gig_h
#define gig_h

#include <Rinternals.h> 
#include <R_ext/Rdynload.h>

//lambda = p, psi=a, chi=b

inline double logsumexp(const double &a, const double &b){
  return a < b ? b + log(1.0 + exp(a - b)) : a + log(1.0 + exp(b - a));
}

SEXP rgig0(int n, double lambda, double chi, double psi);
double rgig(double lambda, double chi, double psi);
double gig_norm(double lambda, double chi, double psi);

#endif
