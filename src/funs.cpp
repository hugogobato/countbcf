#include "arma_config.h"
#include <RcppArmadillo.h>

#include <cmath>
#include "funs.h"
#include "rng.h"
#include <map>
#ifdef MPIBART
#include "mpi.h"
#endif

#include "slice.h"

using Rcpp::Rcout;
using namespace arma;
using namespace Rcpp;

// find variables n can split on, put their indices in goodvars
void getbadvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& badvars)
{
  int L,U;
  for(size_t v=0;v!=xi.size();v++) { // try each variable
    L=0; U = xi[v].size()-1;
    n->rg(v,&L,&U);
    if(U<L) badvars.push_back(v);
  }
}


void update_dart(std::vector<tree>& trees, pinfo& pi, dinfo& di, xinfo& xi, RNG gen) {
  std::vector<double> var_use_counts(di.p, 0);
  
  for(size_t j=0;j<trees.size();j++) {
    std::vector<tree*> tj_nodes;
    trees[j].getnobots(tj_nodes);
    for(size_t b=0; b<tj_nodes.size(); ++b) {
      var_use_counts[tj_nodes[b]->getv()] += 1.0;
      std::vector<size_t> badvars; //variables this tree couldn't split on
      getbadvars(tj_nodes[b], xi, badvars);
      if(badvars.size()>0) {
        double total_prob = 0;
        vector<double> sample_probs(badvars.size());
        for(size_t r = 0; r<badvars.size(); ++r) {
          total_prob += pi.var_probs[badvars[r]];
          sample_probs[r] = pi.var_probs[badvars[r]];
        }
        for(size_t r = 0; r<badvars.size(); ++r) sample_probs[r]/=total_prob;
        double numaug = R::rnbinom(1, 1-total_prob);
        int tt = 0;
        while( (numaug>0.5) & (tt<50000) ) {
          var_use_counts[badvars[rdisc(&sample_probs[0], gen)]] += 1;
          numaug -= 1; tt += 1;
        }
      }
    }
  }

  double tot = 0;
  for(size_t v = 0; v < di.p; ++v) {
    pi.var_probs[v] = R::rgamma(var_use_counts[v] + pi.alpha/(di.p * pi.var_sizes[v]), 1);
    tot += pi.var_probs[v];
  }
  for(size_t v = 0; v < di.p; ++v) pi.var_probs[v] /= tot;

  ld_dart_alpha logdens(var_use_counts, pi.var_sizes, pi.var_probs);

  double a0 = pi.dart_alpha;
  pi.dart_alpha = slice(a0, &logdens, 1.0, INFINITY, 0.0, INFINITY);

}

void update_scale(double* r, double* fits, size_t n, double sigma, pinfo& pi, RNG& gen) {

  // assumes that the scale multiplier is contained in fits!

  double a = 0.0; double b = 0.0;
  // double sigma = pi.sigma*fabs(pi.eta);
  for(size_t k=0;k<n;k++) {
    // double r = y[k] - allfit_con[k];
    // if(randeff) r -= allfit_random[k];

    a += fits[k] * r[k]/pi.eta;
    b += fits[k] * fits[k]/(pi.eta*pi.eta);
  }

  double var = 1 / (1 / (pi.gamma * pi.gamma) + b/(sigma*sigma));
  double mean = var * (0.0/(pi.gamma*pi.gamma) + a/(sigma*sigma));

  pi.eta = mean + gen.normal(0.,1.) * sqrt(var);

  if(pi.scale_df>0) {
    pi.gamma = sqrt(0.5*(pi.scale_df+pi.eta*pi.eta)/gen.gamma((1.0 + pi.scale_df)/2.0, 1.0));
  }

}

//todo: make this a template or change types of allfits elsewhere
void update_scale(double* r, std::vector<double> &fits, size_t n, double sigma, pinfo& pi, RNG& gen) {
  
  // assumes that the scale multiplier is contained in fits!
  
  double a = 0.0; double b = 0.0;
  //double sigma = pi.sigma*fabs(pi.eta);
  for(size_t k=0;k<n;k++) {
    //double r = y[k] - allfit_con[k];
    //if(randeff) r -= allfit_random[k];
    
    a += fits[k] * r[k]/pi.eta;
    b += fits[k] * fits[k]/(pi.eta*pi.eta);
  }
  
  double var = 1 / (1 / (pi.gamma * pi.gamma) + b/(sigma*sigma));
  double mean = var * (0.0/(pi.gamma*pi.gamma) + a/(sigma*sigma));
  
  pi.eta = mean + gen.normal(0.,1.) * sqrt(var);
  
  if(pi.scale_df>0) {
    pi.gamma = sqrt(0.5*(pi.scale_df+pi.eta*pi.eta)/gen.gamma((1.0 + pi.scale_df)/2.0, 1.0));
  }
  
}

//-------------------------------------------------------------
// Generates realizations from multivariate normal.
//-------------------------------------------------------------
mat rmvnormArma(int n, vec mu, mat sigma) {
   //-------------------------------------------------------------
   // INPUTS:	   n = sample size
   //				   mu = vector of means
   //				   sigma = covariance matrix
   //-------------------------------------------------------------
   // OUTPUT:	n realizations of the specified MVN.
   //-------------------------------------------------------------
   int ncols = sigma.n_cols;
   mat Y = randn(n, ncols);
   mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
   return result;
}

struct cmpdouble {
  bool operator()(const double &a, const double &b) const {
    return fabs(a-b)<1e-8;
  }
};

//Basis

void allsuff_basis(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, 
                   std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;  
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    //sv[i].n_vec.resize(di.tlen);
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0.0;
    sv[i].n_unique = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }
  
  std::vector<std::map<int, int> > unique_groups(sv.size());

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
    y = di.y[i];           // Resolves to r.

    //tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    ++(sv[ni].n);
    if(di.group) unique_groups[ni][di.groups[i]] += 1;

    if(di.basis_dim == 1) {
      omega = di.omega[i];
      sv[ni].sy += omega*y;
      sv[ni].n0 += omega*omega;
    } else {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        sv[ni].WtW.at(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          //sv[ni].WtW(j,g) += a; //this is faster than below?
          //sv[ni].WtW[i,j] = sv[ni].WtW[i,j] + a;
          //sv[ni].WtW[j,i] = sv[ni].WtW[j,i] + a; //get this outside obs loop later
        }
      }
    }

  }
  
  
  for(size_t q=0; q<sv.size(); ++q) {
    if(di.group) {   
      sv[q].n_unique = unique_groups[q].size();
    } else {
      sv[q].n_unique = sv[q].n;
    }
  }

  if(di.basis_dim>1) {
    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW.at(j,g) = sv[q].WtW.at(g,j);
        }
      }
    }
  }

  if(di.basis_dim<2) {
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n0;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }


}

void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv,
                         std::vector<sinfo>& sv, sinfo& sl, sinfo& sr, std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0;
    sv[i].n_unique = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }
  
  std::vector<std::map<int, int> > unique_groups(sv.size());
  std::vector<std::map<int, int> > unique_groups_lr(2); //1st elt is l node, second is r node

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  double *omega_i_tmp;
  arma::vec omega_i;
  bool in_candidate_nog, left, right;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    y = di.y[i];

    //tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    left = false; right = false;

    in_candidate_nog = (tbn == nx);
    if(in_candidate_nog) {
      left  = (xx[v] < xi[v][c]);
      right = !(xx[v] < xi[v][c]);
    }

    ++(sv[ni].n);
    if(di.group) unique_groups[ni][di.groups[i]] += 1;
    
    if(left) {
      sl.n += 1;
      if(di.group) unique_groups_lr[0][di.groups[i]] += 1;
    }
    if(right) {
      sr.n += 1;
      if(di.group) unique_groups_lr[1][di.groups[i]] += 1;
    }

    if(di.basis_dim == 1) {
      omega = di.omega[i];
      sv[ni].sy += omega*y;
      sv[ni].n0 += omega*omega;
      if(left) {
        sl.sy += omega*y;
        sl.n0 += omega*omega;
      }
      if(right) {
        sr.sy += omega*y;
        sr.n0 += omega*omega;
      }

    } else {
      omega_i_tmp = di.omega + i*di.basis_dim;
      //get design vector
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;
      if(left) sl.sy_vec += y*omega_i;
      if(right) sr.sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        double a = omega_i_tmp[j]*omega_i_tmp[j];
        sv[ni].WtW.at(j,j) += a;
        if(left) sl.WtW.at(j,j) += a;
        if(right) sr.WtW.at(j,j) += a;
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          if(left) sl.WtW.at(g,j) += a;
          if(right) sr.WtW.at(g,j) += a;
        }
      }
    }
  }
  
  for(size_t q=0; q<sv.size(); ++q) {
    if(di.group) {   
      sv[q].n_unique = unique_groups[q].size();
    } else {
      sv[q].n_unique = sv[q].n;
    }
  }

  if(di.group) {   
    sl.n_unique = unique_groups_lr[0].size();
    sr.n_unique = unique_groups_lr[1].size();
  } else {
    sl.n_unique = sl.n;
    sr.n_unique = sr.n;
  }
  
  
  
  if(di.basis_dim>1) {

    //Rcout << "L" << endl << sl.WtW << endl << sl.sy_vec << endl;
    //Rcout << "R" << endl << sr.WtW << endl << sr.sy_vec << endl;

    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW(j,g) = sv[q].WtW(g,j);
          sl.WtW(j,g) = sl.WtW(g,j);
          sr.WtW(j,g) = sr.WtW(g,j);
        }
      }
    }

    //Rcoutt << "L" << endl << sl.WtW << endl;
    //Rcoutt << "R" << endl << sr.WtW << endl;
  }

  if(di.basis_dim<2) {
    sl.WtW.at(0,0) = sl.n0;
    sl.sy_vec.at(0) = sl.sy;
    sr.WtW.at(0,0) = sr.n0;
    sr.sy_vec.at(0) = sr.sy;
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n0;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }
}

void getsuff_basis(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  bool orig = false;

  double omega;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;


    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

      y = di.y[i];   // extract current yi.  resolves to r.

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      if(xx[v] < xi[v][c]) { // Update left.

        ++(sl.n);

        if(di.basis_dim==1) {
          omega = di.omega[i];
          sl.sy += omega*y;
          sl.n0 += omega*omega;

        } else {

          sl.sy_vec += y*omega_i;
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sl.WtW[i,j] = sl.WtW[i,j] + a;
              //sl.WtW(j,i) += a;
            }
          }

        }


      } else { //Update right.

        ++(sr.n);

        if(di.basis_dim==1) {
          omega = di.omega[i];
          sr.sy += omega*y;
          sr.n0 += omega*omega;

        } else {

          sr.sy_vec += y*omega_i;
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sr.WtW[i,j] = sr.WtW[i,j] + a;
              //sr.WtW(j,i) += a;
            }
          }

        }


      }
    }
  }


  if(di.basis_dim<2) {
    sl.WtW[0,0] = sl.n0;
    sl.sy_vec[0] = sl.sy;
    sr.WtW[0,0] = sr.n0;
    sr.sy_vec[0] = sr.sy;
  }

  // this is faster without bounds checking, but incrementing w/ +=above is faster *with* bounds checking?
  for(size_t j=0; j<di.basis_dim; ++j) {
    //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
    for(size_t i=0; i<j; ++i) {
      sl.WtW[j,i] = sl.WtW[i,j];
      sr.WtW[j,i] = sr.WtW[i,j];
    }
  }

}

/*
void getsuff_basis(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  bool orig = false;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;


    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

      y = di.y[i];   // extract current yi.  resolves to r.

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);


      if(xx[v] < xi[v][c]) { // Update left.

        ++(sl.n);
        sl.sy_vec += y*omega_i;

        if(orig) {
          sl.WtW += omega_i*omega_i.t();
        } else {
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sl.WtW[i,j] = sl.WtW[i,j] + a;
              //sl.WtW(j,i) += a;
            }
          }
        }



      } else { //Update right.
        ++(sr.n);
        sr.sy_vec += y*omega_i;

        if(orig) {
          sr.WtW += omega_i*omega_i.t();
        } else {
          for(size_t j=0; j<di.basis_dim; ++j) {
            //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
            sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
            for(size_t i=0; i<j; ++i) {
              double a = omega_i_tmp[j]*omega_i_tmp[i];
              sr.WtW[i,j] = sr.WtW[i,j] + a;
              //sr.WtW(j,i) += a;
            }
          }
        }

      }
    }
  }

  // this is faster without bounds checking, but incrementing w/ +=above is faster *with* bounds checking
  for(size_t j=0; j<di.basis_dim; ++j) {
    //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
    for(size_t i=0; i<j; ++i) {
      sl.WtW[j,i] = sl.WtW[i,j];
      sr.WtW[j,i] = sr.WtW[i,j];
    }
  }

}
*/
void getsuff_basis(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double *xx;//current x
  double y;  //current y
  double t;  //current t

  bool orig = false;

  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  double omega;

  sl.WtW = zeros(di.basis_dim,di.basis_dim); sl.sy_vec = zeros(di.basis_dim);
  sr.WtW = zeros(di.basis_dim,di.basis_dim); sr.sy_vec = zeros(di.basis_dim);

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    tree::tree_cp bn = x.bn(xx,xi);

    y = di.y[i];   // extract current yi.

    if(bn==nl) {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      ++(sl.n);

      if(di.basis_dim==1) {
        omega = di.omega[i];
        sl.sy += omega*y;
        sl.n0 += omega*omega;

      } else {

        sl.sy_vec += y*omega_i;
        for(size_t j=0; j<di.basis_dim; ++j) {
          //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
          sl.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
          for(size_t i=0; i<j; ++i) {
            double a = omega_i_tmp[j]*omega_i_tmp[i];
            sl.WtW[i,j] = sl.WtW[i,j] + a;
            //sl.WtW(j,i) += a;
          }
        }

      }

    }

    if(bn==nr) {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      ++(sr.n);

      if(di.basis_dim==1) {
        omega = di.omega[i];
        sr.sy += omega*y;
        sr.n0 += omega*omega;

      } else {

        sr.sy_vec += y*omega_i;
        for(size_t j=0; j<di.basis_dim; ++j) {
          //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
          sr.WtW(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
          for(size_t i=0; i<j; ++i) {
            double a = omega_i_tmp[j]*omega_i_tmp[i];
            sr.WtW[i,j] = sr.WtW[i,j] + a;
            //sr.WtW(j,i) += a;
          }
        }

      }

    }
  }


  if(di.basis_dim<2) {
    sl.WtW[0,0] = sl.n0;
    sl.sy_vec[0] = sl.sy;
    sr.WtW[0,0] = sr.n0;
    sr.sy_vec[0] = sr.sy;
  }

  for(size_t j=0; j<di.basis_dim; ++j) {
    //sv[ni].sy_vec(j) += y*omega_i_tmp[j];
    for(size_t i=0; i<j; ++i) {
      sl.WtW[j,i] = sl.WtW[i,j];
      sr.WtW[j,i] = sr.WtW[i,j];
    }
  }

}

void drmu_basis(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{

  bool debug = false;
  tree::npv bnv;
  std::vector<sinfo> sv; //will be resized in allsuff
  tree::npv bnv0; std::vector<sinfo> sv0;
  //debug is broken, would need to eat node_pointers. suff stat code works anyhow
//  if(debug) allsuff_basis(t,xi,di,bnv0,sv0,node_pointers);

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  t.getbots(bnv);   // Save bottom nodes for x to bnv variable.
  tree::npv::size_type nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  double prior_prec = pi.Prec0(0,0);
  double prior_mean = pi.mu0(0);

  for(tree::npv::size_type i=0; i<nb; i++) {
    sv[i] = bnv[i]->s;
    if(debug) {
      if(abs(sv[i].sy - sv0[i].sy)>1e-8) {

        Rcout << "depth "<< bnv[i]->depth() << endl;
        Rcout << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
        Rcout << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
        stop("shit");
      }
     else {
        //Rcout << "depth "<< bnv[i]->depth() << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec <<" new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << endl;
        //Rcoutt << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
       //Rcoutt << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
      }
     sv[i] = sv0[i];
    }

  }

  if(di.basis_dim<2) {
    vec beta_draw(1);
    double tt = prior_prec*prior_mean;
    double s2 = (pi.sigma*pi.sigma);
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      double Phi = sv[i].n0/s2 + prior_prec;
      double m = tt + sv[i].sy/s2;
      beta_draw(0) = m/Phi + gen.normal(0,1)/sqrt(Phi); //rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);

      // Check for NA result.
      if(beta_draw(0) != beta_draw(0)) {
        Rcpp::stop("drmu failed");
      }
    }
  } else {
    mat Phi;
    vec m;
    vec beta_draw;

    vec tt = pi.Prec0*pi.mu0;
    double s2 = (pi.sigma*pi.sigma);
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      //Rcoutt << "phi ";
      //Rcoutt << "WtW" << endl << sv[i].WtW << endl;
      Phi = sv[i].WtW/s2 + pi.Prec0;
      //Rcoutt << "m ";
      m = tt + sv[i].sy_vec/s2;
      //Rcoutt << "draw ";
      beta_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
        Rcpp::stop("drmu failed");
      }
    }
  }

}

double lil_basis(sinfo& s, pinfo& pi){

//  Rcout << "log likelihood calc" << endl;

  double ll=0; double tt;

  double s2 = pi.sigma*pi.sigma;

  if(false) {//(s.sy_vec.n_elem == 1) {
    /*

    This is slow as hell??
     Instead: check for 1d when computing summary stats, avoid matrix ops** everywhere**

    double precpred = pi.Prec0.at(0,0) + s.n0/s2;

    //double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    double dotp = precpred*s.sy/(s2*s2);
    double tt = 0;

    tt = -0.5*log(precpred) - 0.5*pi.logdetSigma0;

    ll = -0.5*n*log(s2) + tt + 0.5*dotp;
    */


  } else {
    mat precpred = pi.Prec0 + s.WtW/s2;

    double dotp = dot(s.sy_vec/s2, solve(precpred, s.sy_vec/s2));
    double tt = 0;

    tt = -0.5*log(det(precpred)) - 0.5*pi.logdetSigma0;

    ll = -0.5*s.n*log(s2) + tt + 0.5*dotp;
  }

  //Rcpp::Rcout<< ll << endl;

  return(ll);
}

double lil_basis(sinfo& sl, sinfo& sr, pinfo& pi){

  sinfo st = sl;
  st.n += sr.n;
  st.n0 += sr.n0;
  st.sy += sr.sy;
  st.sy_vec += sr.sy_vec;
  st.WtW += sr.WtW;

  double ll = lil_basis(st, pi);

  //Rcpp::Rcout<< ll << endl;

  return(ll);
}

void fit_basis(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla)
{
  double *xx;
  double *omega_i_tmp;
  tree::tree_cp bn;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    if(populate) {
      bn = t.bn(xx,xi);
      node_pointers[i] = bn;
    } else {
      bn = node_pointers[i];
    }

    omega_i_tmp = di.omega + i*di.basis_dim;
    //arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);
    //fv[i] = arma::dot(bn->getm(), omega_i);

    //ok
    if(vanilla) {
      fv[i] = bn->mu(0);
    } else {
      fv[i] = 0.0;
      for(size_t j=0; j<di.basis_dim; ++j) {
        fv[i] += omega_i_tmp[j]*bn->mu(j);
      }
    }


  }
}

double fit_i_basis(size_t& i, std::vector<tree>& t, xinfo& xi, dinfo& di, bool vanilla)
{
  double *xx;
  double fv = 0.0;
  double *omega_i_tmp;
  tree::tree_cp bn;

  for(size_t j=0;j<t.size();j++) {
    xx = di.x + i*di.p;
    bn = t[j].bn(xx,xi); //instead of dropping, using the node_pointer object

    omega_i_tmp = di.omega + i*di.basis_dim;

    //ok
    if(vanilla) {
      fv += bn->mu(0);
    } else {
      for(size_t k=0; k<di.basis_dim; ++k) {
        fv += omega_i_tmp[k]*bn->mu(k);
      }
    }

  }

  return fv;
}


//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
	double dif = x-m;
	return exp(-.5*dif*dif/v)/sqrt(2*M_PI*v);
}

//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p, RNG& gen)
{

	double sum;
	double u = gen.uniform();

    int i=0;
    sum=p[0];
    while(sum<u) {
		i += 1;
		sum += p[i];
    }
    return i;
}

//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os)
{
	size_t p = xi.size();
	if(p!=2) {
		Rcout << "error in grm, p !=2\n";
		return;
	}
	size_t n1 = xi[0].size();
	size_t n2 = xi[1].size();
	tree::tree_cp bp; //pointer to bottom node
	double *x = new double[2];
	for(size_t i=0;i!=n1;i++) {
		for(size_t j=0;j!=n2;j++) {
			x[0] = xi[0][i];
			x[1] = xi[1][j];
			bp = tr.bn(x,xi);
			os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << endl;
		}
	}
	delete[] x;
}

//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
	int L,U;
	bool v_found = false; //have you found a variable you can split on
	size_t v=0;
	while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) v_found=true;
		v++;
	}
	return v_found;
}

//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
	double pb;  //prob of birth to be returned
	tree::npv bnv; //all the bottom nodes
	t.getbots(bnv);
	for(size_t i=0;i!=bnv.size();i++)
		if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
	if(goodbots.size()==0) { //are there any bottom nodes you can split on?
		pb=0.0;
	} else {
		if(t.treesize()==1) pb=1.0; //is there just one node?
		else pb=pi.pb;
	}
	return pb;
}

//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
	int L,U;
	for(size_t v=0;v!=xi.size();v++) {//try each variable
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) goodvars.push_back(v);
	}
}

//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
	if(cansplit(n,xi)) {
		return pi.alpha/pow(1.0+n->depth(),pi.beta);
	} else {
		return 0.0;
	}
}

//--------------------------------------------------

//get counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
	size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	bnv.clear();
	x.getbots(bnv);

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

  std::vector<int> cts(bnv.size(), 0);

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		y=di.y[i];

		tbn = x.bn(xx,xi);
		ni = bnmap[tbn];

    cts[ni] += 1;
	}
  return(cts);
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   tree::npv& bnv, //vector of pointers to bottom nodes
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index

	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
	double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}


void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign,
                   tree::tree_cp &tbn
                   )
{
  //tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di) {
  bool good = true;
  tree::npv bnv;
  std::vector<int> cts;
  int m = 0;
  for (size_t tt=0; tt<t.size(); ++tt) {
    cts = counts(t[tt], xi, di, bnv);
    m = std::min(m, *std::min_element(cts.begin(), cts.end()));
    if(m<minct) {
      good = false;
      break;
    }
  }
  return good;
}

mat coef_basis(tree& t, xinfo& xi, dinfo& di)
{
  double *xx;

  tree::tree_cp bn;

  mat out(di.basis_dim,di.n);

  for(size_t i=0; i<di.n; i++) {
    xx = di.x + i*di.p;
    bn = t.bn(xx,xi);
    out.col(i) = bn->getm();
    /*
    for(size_t j=0; j<di.basis_dim; ++j) {
      out.at(i,j) = bn->mu.at(j);
    }
     */
  }
  return(out);
}

//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
	double *xx;
	tree::tree_cp bn;
	pv.resize(di.n);
	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		bn = t.bn(xx,xi);
		pv[i] = bn->nid();
	}
}

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
	Rcout << "xinfo: \n";
	for(size_t v=0;v!=xi.size();v++) {
		Rcout << "v: " << v << endl;
		for(size_t j=0;j!=xi[v].size();j++) Rcout << "j,xi[v][j]: " << j << ", " << xi[v][j] << endl;
	}
	Rcout << "\n\n";
}

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
	double xinc;

	//compute min and max for each x
	std::vector<double> minx(p,INFINITY);
	std::vector<double> maxx(p,-INFINITY);
	double xx;
	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}
// get min/max needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xx;

	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
}
//make xinfo = cutpoints give the minx and maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xinc;
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}


double lil_loglinear(sinfo& s, pinfo& pi){

  //  Rcout << "log likelihood calc" << endl;
  
  #define LN2 0.693147180559945309417232121458

  // calculate log-likelihood using gig_norm 

  // define a and b (should be available in prior...)

  double ll=0; 
  ll = -gig_norm(pi.c, 0, 2*pi.d)-LN2;
  double loga = gig_norm(s.sy + pi.c, 0, 2*pi.d + 2*s.n0);
  double logb = gig_norm(s.sy - pi.c, 2*pi.d, 2*s.n0);
  ll += logsumexp(loga, logb);

  //Rcpp::Rcout<< ll << endl;

  return(ll);
}

double lil_loglinear(sinfo& sl, sinfo& sr, pinfo& pi){

  // update!
  sinfo st = sl;
  st.n += sr.n;
  st.n0 += sr.n0;
  st.sy += sr.sy;
  st.sy_vec += sr.sy_vec;
  st.WtW += sr.WtW;

  // double ll = lil_basis(st, pi);
  double ll = lil_loglinear(st, pi);

  //Rcpp::Rcout<< ll << endl;

  return(ll);
}


double ll_loglinear(
  std::vector<double>& y_vals, 
  std::vector<double>& log_mean, 
  double kappa, 
  int model, 
  bool log_scale, 
  std::vector<double>& log_w, 
  int n_zero
){
  // evaluate likelihood for the 4 count models
  double ll_val = 0;

  if (model == 1){
    // poisson
    for(std::size_t i = 0; i < y_vals.size(); ++i){
      ll_val += R::dpois(y_vals[i], exp(log_mean[i]), true);
    }
  }
  else if (model == 2){
    // negative binomial
    for(std::size_t i = 0; i < y_vals.size(); ++i){
      ll_val += R::dnbinom_mu(y_vals[i], kappa, exp(log_mean[i]), true);
    }
  }
  else if (model == 3){
    // zero-inflated poisson

    // if y = 0
    for (std::size_t i = 0; i < n_zero; ++i){
      double log_wp = log_w[i] + R::dpois(y_vals[i], exp(log_mean[i]), true);
      double log_1mw = log1p(-exp(log_w[i]));
      ll_val += logsumexp(log_1mw, log_wp);
    }

    // if y > 0
    for (std::size_t i = n_zero; i < y_vals.size(); ++i){
      ll_val += log_w[i] + R::dpois(y_vals[i], exp(log_mean[i]), true);
    }
  } 
  else if (model == 4){
    // zero-inflated negative binomial

    // if y = 0
    for (std::size_t i = 0; i < n_zero; ++i){
      double log_wp = log_w[i] + R::dnbinom_mu(y_vals[i], kappa, exp(log_mean[i]), true);
      double log_1mw = log1p(-exp(log_w[i]));
      ll_val += logsumexp(log_1mw, log_wp);
    }

    // if y > 0
    for (std::size_t i = n_zero; i < y_vals.size(); ++i){
      ll_val += log_w[i] + R::dnbinom_mu(y_vals[i], kappa, exp(log_mean[i]), true);
    }
  }

  // log-likelihood or likelihood
  double final_vrsn = 0;
  if (log_scale){
    final_vrsn = ll_val;
  }
  else
  {
    final_vrsn = exp(ll_val);
  }

  return (final_vrsn);
}


void fit_loglinear(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla)
{
  double *xx;
  double *omega_i_tmp;
  tree::tree_cp bn;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    if(populate) {
      bn = t.bn(xx,xi);
      node_pointers[i] = bn;
    } else {
      bn = node_pointers[i];
    }

    omega_i_tmp = di.omega + i*di.basis_dim;
    //arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);
    //fv[i] = arma::dot(bn->getm(), omega_i);

    //ok
    if(vanilla) {
      fv[i] = bn->mu(0);
    } else {
      fv[i] = 0.0;
      for(size_t j=0; j<di.basis_dim; ++j) {
        fv[i] += omega_i_tmp[j]*bn->mu(j);
      }
    }
  }
}



void allsuff_loglinear(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, 
                   std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  double logoff;  // log offset
  double ui;      // u_i sufficient stat

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;  
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    //sv[i].n_vec.resize(di.tlen);
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0.0;
    sv[i].n_unique = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }
  
  std::vector<std::map<int, int> > unique_groups(sv.size());

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
    y = di.y[i];           // Resolves to r.

    // additional data for sufficient stat calculation
    logoff = di.offset[i];
    ui = di.u_i[i];

    //tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    ++(sv[ni].n);
    if(di.group) unique_groups[ni][di.groups[i]] += 1;

    if(di.basis_dim == 1) {
      // omega = di.omega[i];
      //sv[ni].sy += omega*y;
      //sv[ni].n0 += omega*omega;

      sv[ni].sy += ui;
      sv[ni].n0 += y;  // exp(logoff + y);
    }
    else
    {

      //get design vector
      double *omega_i_tmp = di.omega + i*di.basis_dim;
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        sv[ni].WtW.at(j,j) += omega_i_tmp[j]*omega_i_tmp[j];
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          //sv[ni].WtW(j,g) += a; //this is faster than below?
          //sv[ni].WtW[i,j] = sv[ni].WtW[i,j] + a;
          //sv[ni].WtW[j,i] = sv[ni].WtW[j,i] + a; //get this outside obs loop later
        }
      }
    }
  }
  
  
  for(size_t q=0; q<sv.size(); ++q) {
    if(di.group) {   
      sv[q].n_unique = unique_groups[q].size();
    } else {
      sv[q].n_unique = sv[q].n;
    }
  }

  if(di.basis_dim>1) {
    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW.at(j,g) = sv[q].WtW.at(g,j);
        }
      }
    }
  }

  if(di.basis_dim<2) {
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n0;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }


}


void allsuff_loglinear_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv,
                         std::vector<sinfo>& sv, sinfo& sl, sinfo& sr, std::vector<tree::tree_cp>& node_pointers)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  double t;          //current t
  double omega;

  double logoff;    // log(offset)
  double ui;        // u_i sufficient statistics

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  // Resize vectors within sufficient stats to have di.tlen length.
  for(size_t i = 0; i < nb; ++i){
    sv[i].sy = 0;
    sv[i].n0 = 0.0;
    sv[i].n = 0;
    sv[i].n_unique = 0;
    sv[i].sy_vec.zeros(di.basis_dim);
    sv[i].WtW.zeros(di.basis_dim, di.basis_dim);
  }
  
  std::vector<std::map<int, int> > unique_groups(sv.size());
  std::vector<std::map<int, int> > unique_groups_lr(2); //1st elt is l node, second is r node

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2

  double *omega_i_tmp;
  arma::vec omega_i;
  bool in_candidate_nog, left, right;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    y = di.y[i];

    ui = di.u_i[i];
    logoff = di.offset[i];

    // tbn = x.bn(xx,xi); // Find bottom node for this observation.
    tbn = node_pointers[i];
    ni = bnmap[tbn];   // Map bottom node to integer index

    left = false; right = false;

    in_candidate_nog = (tbn == nx);
    if(in_candidate_nog) {
      left  = (xx[v] < xi[v][c]);
      right = !(xx[v] < xi[v][c]);
    }

    ++(sv[ni].n);
    if(di.group) unique_groups[ni][di.groups[i]] += 1;
    
    if(left) {
      sl.n += 1;
      if(di.group) unique_groups_lr[0][di.groups[i]] += 1;
    }
    if(right) {
      sr.n += 1;
      if(di.group) unique_groups_lr[1][di.groups[i]] += 1;
    }

    if(di.basis_dim == 1) {
      //omega = di.omega[i];
      //sv[ni].sy += omega*y;
      //sv[ni].n0 += omega*omega;

      sv[ni].sy += ui;
      sv[ni].n0 += y; 

      if(left) {
        // sl.sy += omega*y;
        // sl.n0 += omega*omega;

        sl.sy += ui;
        sl.n0 += y; 
      }
      if(right) {
        // sr.sy += omega*y;
        // sr.n0 += omega*omega;

        sr.sy += ui;
        sr.n0 += y; 
      }
    } else {
      omega_i_tmp = di.omega + i*di.basis_dim;
      //get design vector
      arma::vec omega_i(omega_i_tmp, di.basis_dim, false, false);

      sv[ni].sy_vec += y*omega_i;
      if(left) sl.sy_vec += y*omega_i;
      if(right) sr.sy_vec += y*omega_i;

      for(size_t j=0; j<di.basis_dim; ++j) {
        double a = omega_i_tmp[j]*omega_i_tmp[j];
        sv[ni].WtW.at(j,j) += a;
        if(left) sl.WtW.at(j,j) += a;
        if(right) sr.WtW.at(j,j) += a;
        for(size_t g=0; g<j; ++g) {
          double a = omega_i_tmp[j]*omega_i_tmp[g];
          sv[ni].WtW.at(g,j) += a;
          if(left) sl.WtW.at(g,j) += a;
          if(right) sr.WtW.at(g,j) += a;
        }
      }
    }
  }
  
  for(size_t q=0; q<sv.size(); ++q) {
    if(di.group) {   
      sv[q].n_unique = unique_groups[q].size();
    } else {
      sv[q].n_unique = sv[q].n;
    }
  }

  if(di.group) {   
    sl.n_unique = unique_groups_lr[0].size();
    sr.n_unique = unique_groups_lr[1].size();
  } else {
    sl.n_unique = sl.n;
    sr.n_unique = sr.n;
  }
  
  if(di.basis_dim>1) {

    //Rcout << "L" << endl << sl.WtW << endl << sl.sy_vec << endl;
    //Rcout << "R" << endl << sr.WtW << endl << sr.sy_vec << endl;

    for(size_t q=0; q<sv.size(); ++q) {
      for(size_t j=0; j<di.basis_dim; ++j) {
        for(size_t g=0; g<j; ++g) {
          sv[q].WtW(j,g) = sv[q].WtW(g,j);
          sl.WtW(j,g) = sl.WtW(g,j);
          sr.WtW(j,g) = sr.WtW(g,j);
        }
      }
    }

    //Rcoutt << "L" << endl << sl.WtW << endl;
    //Rcoutt << "R" << endl << sr.WtW << endl;
  }

  if(di.basis_dim<2) {
    sl.WtW.at(0,0) = sl.n0;
    sl.sy_vec.at(0) = sl.sy;
    sr.WtW.at(0,0) = sr.n0;
    sr.sy_vec.at(0) = sr.sy;
    for(size_t j=0; j<sv.size(); ++j) {
      sv[j].WtW.at(0,0) = sv[j].n0;
      sv[j].sy_vec.at(0) = sv[j].sy;
    }
  }
}


void drmu_loglinear(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{

  bool debug = false;
  tree::npv bnv;
  std::vector<sinfo> sv; //will be resized in allsuff
  tree::npv bnv0; std::vector<sinfo> sv0;
  //debug is broken, would need to eat node_pointers. suff stat code works anyhow
  //  if(debug) allsuff_basis(t,xi,di,bnv0,sv0,node_pointers);

  bnv.clear();      // Clear the bnv variable if any value is already saved there.
  t.getbots(bnv);   // Save bottom nodes for x to bnv variable.
  tree::npv::size_type nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

  double prior_prec = pi.Prec0(0,0);
  double prior_mean = pi.mu0(0);

  double leaf_c = pi.c;
  double leaf_d = pi.d;

  for(tree::npv::size_type i=0; i<nb; i++) {
    sv[i] = bnv[i]->s;
    if(debug) {
      if(abs(sv[i].sy - sv0[i].sy)>1e-8) {

        Rcout << "depth "<< bnv[i]->depth() << endl;
        Rcout << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
        Rcout << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
        stop("shit");
      }
     else {
        //Rcout << "depth "<< bnv[i]->depth() << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec <<" new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << endl;
        //Rcoutt << " good " << sv0[i].sy << " " << sv0[i].n << " " << sv0[i].sy_vec << sv0[i].WtW << endl;
       //Rcoutt << " new " << sv[i].sy << " " << sv[i].n << " " << sv[i].sy_vec << sv[i].WtW << endl;
      }
     sv[i] = sv0[i];
    }

  }

  if(di.basis_dim<2) {
    vec beta_draw(1);
    double tt = prior_prec*prior_mean;
    double s2 = (pi.sigma*pi.sigma);
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      // double Phi = sv[i].n0/s2 + prior_prec;
      // double m = tt + sv[i].sy/s2;
      // beta_draw(0) = m/Phi + gen.normal(0,1)/sqrt(Phi); //rmvnorm_post(m, Phi);

      // compute weights
      double loga = gig_norm(sv[i].sy + leaf_c, 0, 2*leaf_d + 2*sv[i].n0);
      double logb = gig_norm(sv[i].sy - leaf_c, 2*leaf_d, 2*sv[i].n0);
      double logapb = logsumexp(loga, logb);

      double mu;

      if(log(gen.uniform())<loga-logapb) {
        // gamma
        mu = gen.gamma(sv[i].sy + leaf_c, 1 / (sv[i].n0 + leaf_d));
      } else {
        // gig
        mu = gen.gig(sv[i].sy - leaf_c, 2 * leaf_d, 2 * sv[i].n0);
      }

      // save log(mu)
      beta_draw(0) = log(mu);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);

      // Check for NA result.
      if(beta_draw(0) != beta_draw(0)) {
        Rcpp::stop("drmu failed");
      }
    }
  } else {
    mat Phi;
    vec m;
    vec beta_draw;

    vec tt = pi.Prec0*pi.mu0;
    double s2 = (pi.sigma*pi.sigma);
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      //Rcoutt << "phi ";
      //Rcoutt << "WtW" << endl << sv[i].WtW << endl;
      Phi = sv[i].WtW/s2 + pi.Prec0;
      //Rcoutt << "m ";
      m = tt + sv[i].sy_vec/s2;
      //Rcoutt << "draw ";
      beta_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(beta_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
        Rcpp::stop("drmu failed");
      }
    }
  }
}


void drxi_loglinear(std::vector<double>& log_xi, double kappa, std::vector<double>& log_mu, std::vector<double>& y, std::vector<int>& z, RNG& gen)
{
  // sample new xi_i from Gamma(kappa + y_i, kappa + mu_{0i}f(x_i))
  int num_y = y.size();

  for (size_t i = 0; i < num_y; ++i) {
    // sample new xi if Z_i = 1
    if (z[i] == 1){
      double alpha = kappa + y[i];
      double beta = kappa + exp(log_mu[i]);

      double xi_i = gen.gamma(alpha, 1 / beta);
      log_xi[i] = log(xi_i);
    }
  }
}


void drphi_loglinear(std::vector<double>& log_phi, std::vector<double>& log_f_sum, RNG& gen)
{
  // sample new phi_i from Exp(f_0(x_i) + f_1(x_i))
  int num_phi = log_phi.size();

  for (size_t i = 0; i < num_phi; ++i) {
    double phi_i = gen.exponential(exp(-log_f_sum[i]));
    log_phi[i] = log(phi_i);
  }
}

void drz_loglinear(std::vector<int>& z, std::vector<double>& log_w, double kappa, std::vector<double>& log_mu, int n_zero, int model, RNG& gen)
{
  // sample new z_i from Bern(p_i), with
  //    p_i = w_i * p_y(y = 0| x_i, kappa, f) / (1 - w_i + w_i * p_y(y = 0| x_i, kappa, f))

  for (size_t i = 0; i < n_zero; ++i) {

    // log-likelihood of y = 0 under count model
    double ll_zero = 0.0;

    if (model == 3){
      // poisson
      ll_zero = R::dpois(0, exp(log_mu[i]), true);
    }
    else if (model == 4){
      // negative binomial
      ll_zero = R::dnbinom_mu(0, kappa, exp(log_mu[i]), true);
    }

    double log_num = log_w[i] + ll_zero;
    double log_denom = logsumexp(log1p(-exp(log_w[i])), log_num);
    double prob_i = exp(log_num - log_denom);

    // sample new z values
    z[i] = gen.binom(1, prob_i);
  }
}
