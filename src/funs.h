#ifndef GUARD_funs_h
#define GUARD_funs_h

#include "arma_config.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rng.h"

using namespace arma;
/*
inline double logsumexp(const double &a, const double &b){
  return a < b ? b + log(1.0 + exp(a - b)) : a + log(1.0 + exp(b - a));
}
*/
using std::cout;
using std::endl;

//pi and log(2*pi)
//#define PI 3.1415926535897931
#define LTPI 1.83787706640934536


typedef std::vector<std::vector<int> > lookup_t;

void getbadvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& badvars);
void update_scale(double* r, double* fits, size_t n, double sigma, pinfo& pi, RNG& gen);
void update_scale(double* r, std::vector<double>& fits, size_t n, double sigma, pinfo& pi, RNG& gen);
void update_dart(std::vector<tree>& trees, pinfo& pi, dinfo& di, xinfo& xi, RNG gen);

lookup_t make_lookup(Rcpp::IntegerMatrix lookup_table, Rcpp::IntegerVector cx);

// functions for basis bart
void fit_basis(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla);
double fit_i_basis(size_t& i, std::vector<tree>& t, xinfo& xi, dinfo& di, bool vanilla);
void allsuff_basis(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, std::vector<tree::tree_cp>& node_pointers);
void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv,
                         std::vector<sinfo>& sv, sinfo& sl, sinfo& sr, std::vector<tree::tree_cp>& node_pointers);
void getsuff_basis(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
void getsuff_basis(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
void drmu_basis(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
//double lil_basis(double n, vec sy_vec, mat WtW, double sigma, vec mu0, mat Prec0, mat Sigma0);
double lil_basis(sinfo& s, pinfo& pi);
double lil_basis(sinfo& sl, sinfo& sr, pinfo& pi);

mat coef_basis(tree &t, xinfo &xi, dinfo &di);

// functions for count bart
double lil_loglinear(sinfo& s, pinfo& pi);
double lil_loglinear(sinfo& sl, sinfo& sr, pinfo& pi);
void drmu_loglinear(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
void drxi_loglinear(std::vector<double>& log_xi, double kappa, std::vector<double>& log_mu, std::vector<double>& y, std::vector<int>& z, RNG& gen);
void drz_loglinear(std::vector<int>& z, std::vector<double>& log_w, double kappa, std::vector<double>& log_mu, int n_zero, int model, RNG& gen);
void drphi_loglinear(std::vector<double>& log_phi, std::vector<double>& log_f_sum, RNG& gen);
double ll_loglinear(std::vector<double>& y_vals, std::vector<double>& log_mean, double kappa, int model, bool log_scale, std::vector<double>& log_w, int n_zero);
void fit_loglinear(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla);
void allsuff_loglinear(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, 
                   std::vector<tree::tree_cp>& node_pointers);
void allsuff_loglinear_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv,
                         std::vector<sinfo>& sv, sinfo& sl, sinfo& sr, std::vector<tree::tree_cp>& node_pointers);

//--------------------------------------------------
//normal density
double pn(
   double x,    //variate
   double m,    //mean
   double v     //variance
);
//--------------------------------------------------
//draw from a discrete distribution
int rdisc(
   double *p,   //vector of probabilities
   RNG& gen     //random number generator
);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);

//--------------------------------------------------
//get counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di);
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv);
//--------------------------------------------------
//update counts (inc or dec) to reflect observation i
// deprecated:
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, int sign);
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, int sign);

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp,size_t>& bnmap, int sign);
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp,size_t>& bnmap, int sign, tree::tree_cp &tbn);

//--------------------------------------------------
//check minimum leaf size
bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di);

//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv);
//--------------------------------------------------
// draw all the bottom node mu's

void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
void drmuhet(tree& t, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen);
void drphi(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//get min/max for p predictors needed to make cutpoints.
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx);
//make xinfo = cutpoints given minx/maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx);

//--------------------------------------------------
// Check if a vector is sorted.  For checking z and zpred for causal funbart.
bool is_sort(arma::vec x);

// //--------------------------------------------------
// //log of the integreted likelihood
// double lil(double n, double sy, double sy2, double sigma, double tau);
// //lilhet drops constants that cancel in the mh ratio, lil/lilprec don't (legacy)
// double lilhet(double n, double sy, double sy2, double sigma, double tau);
// //sy isn't needed, but convenient to maintain fcn sig
// double lilprec(double n, double sy, double sy2, double sigma, double tau);


// //--------------------------------------------------
// //get sufficients stats for all bottom nodes
// void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv);
// void allsuffhet(tree& x, xinfo& xi, dinfo& di, double* phi, tree::npv& bnv, std::vector<sinfo>& sv);

// //--------------------------------------------------
// //get sufficient stats for children (v,c) of node nx in tree x
// void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
// void getsuffhet(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);
// //--------------------------------------------------
// //get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
// void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
// void getsuffhet(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);

#endif
