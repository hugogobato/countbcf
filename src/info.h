#ifndef GUARD_info_h
#define GUARD_info_h

#include "arma_config.h"
#include <RcppArmadillo.h>
using namespace arma;

//============================================================
//data
//============================================================

class dinfo {
public:

   // Methods
   size_t p;  //number of vars
   size_t n;  //number of observations
   size_t basis_dim;
   double *x; // jth var of ith obs is *(x + p*i+j).  *x is first element of the data - is a pointer.
   double *omega; //matrix of basis fcns in same form as x
   //arma::mat *omega_mat; //arma mat version of same dim
   double *y; // 'residual' vector

   bool group; //Whether observations should be grouped by di.group when counting "datapoints"
   int *groups; // If obs are grouped for the purpose of counting unique "datapoints" in BART prior

   double *t; // ith t is *(t+i) or t[i].  (*t)[i] to reference value.  (Points to first element of t vector.)
   size_t tlen; // size_t T;  // Length of vector t_ref; total number of time points.

   //-------------------//
   // COUNT BART UPDATE //
   //-------------------//
   double *u_i; // u_i from expression (5), log-linear BART paper
   double *offset; // log(offset) vector 

   // Constructor
   dinfo() {
     p = 0;
     n = 0;
     basis_dim = 0;
     x = 0;
     y = 0;
     u_i = 0;
     offset = 0;
     groups = 0;
     omega = 0;
     tlen = 1;
     t = 0;
     group = false;
   }
};

//============================================================
//prior and mcmc
//============================================================

class pinfo
{
public:

   //----------------------------
   // Declare properties.
   //----------------------------

   //mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
   //prior info
   double alpha;
   double beta;
   // leaf hyperparameters
   double c;
   double d;
   // sigma
   double sigma;

   bool vanilla; // "standard" BART model, basis functions == 1

   // stand-in for now
   double tau;

   pinfo() {
     pbd=1.0; 
     pb=.5; 
     alpha=.95; 
     beta=.5;
     c = 5.0;
     d = 3.0;
     tau = 1.0;
     eta = 1.0;
     sigma = 1.0;
     dart = false;
     vanilla = false;
   }

   arma::vec mu0;       // Vector of prior means; length tlen.
   mat Sigma0;          // Prior cov matrix for mu_l ~iid N(mu0, Sigma0 = K^-1)
   mat Prec0;           // Prior precision matrix for mu_l, ie inv(Sigma0).
   double logdetSigma0; // log(det(Sigma0))
   double scale_df;     // df for half t hyperprior on function scale
   
   //DART
   bool dart;
   std::vector<double> var_probs; //store splitting probabilities for each variable
   std::vector<double> var_sizes; //store sizes, for scaling dummies
   double dart_alpha; //parameter for dirichlet prior with params dart_alpha/(p*var_sizes[j])

   // For augmented model to induce C+(0,sig2) hyperprior on var.
   // y = eta * f(x) + e, e ~ N(0,sig2)
   // eta | gamma ~ N(0,gamma2)
   double eta;
   double gamma;

   pinfo(size_t tlen) { 
    pbd=1.0; pb=.5; alpha=.95; beta=0.5; sigma=1.0;
    tau = 1.0;
    mu0 = zeros(tlen);
    eta = 1.0; gamma = 1.0;
    Sigma0 = eye(tlen,tlen);
    Prec0 = eye(tlen,tlen);
    vanilla = false;
   }

};

//============================================================
//sufficient statistics for 1 node
//============================================================

class sinfo
{
public:
   bool cache;
   double n0; // unweighted sample counts for het case.
   double n;
   double n_unique;
   double sy;
   double sy2;
   vec n0_vec; // unweighted sample counts for het case at each t.
   vec n_vec;
   vec sy_vec; // will be W'y for basis
   mat WtW;    // W'W

   // Original constructor; no length argument.
   sinfo() {n0=0.0;n=0;sy=0.0;sy2=0.0; n0_vec = zeros(1); n_vec = zeros(1); sy_vec = zeros(1);WtW = zeros(1,1); cache=true;}

   // New constructor, with length argument.
   sinfo(size_t tlen) { 
     n0=0.0;n=0;sy=0.0;sy2=0.0; 
     n0_vec = zeros(tlen); n_vec = zeros(tlen); sy_vec = zeros(tlen);
     WtW = zeros(tlen, tlen);
     cache=true;
     }
};

#endif
