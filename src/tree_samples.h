#ifndef GUARD_ts_h
#define GUARD_ts_h


//#include <Rcpp.h>
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using namespace Rcpp;

class tree_samples {
  public:
  bool init;
  size_t m,p,ndraws,basis_dim;
	xinfo xi;
	std::vector<std::vector<tree> > t;
	mutable std::vector<std::vector<flat_tree> > t_flat;
  
  void load_string(CharacterVector samples_);
  CharacterVector save_string();
  void load_list(List samples_, List x_info_list);
  arma::cube coefs(NumericMatrix &x_);
  arma::mat fits(NumericMatrix &x_, NumericMatrix &Omega);
  void scale(double s);
  
  tree_samples() {init=false;}
  tree_samples(size_t m_, size_t p_, size_t ndraws_, 
               size_t basis_dim_, xinfo& xi_) {
    m = m_;
    p = p_;
    ndraws = ndraws_;
    basis_dim = basis_dim_;
    xi = xi_;
    
    t.resize(ndraws,std::vector<tree>(m));
    
    init = true;
    
  }
  
  template<class Archive>
  void save(Archive & archive) const
  {
    t_flat.resize(t.size());
    for(size_t i=0; i<t_flat.size(); ++i) {
      t_flat[i].resize(t[i].size());
      for(size_t j=0; j<t_flat[i].size(); ++j) {
        t_flat[i][j] = t[i][j].flatten();
      }
    }
    archive( m, p, ndraws, basis_dim, xi, t_flat, init ); 
  }
  
  template<class Archive>
  void load(Archive & archive)
  {
    archive( m, p, ndraws, basis_dim, xi, t_flat, init ); 
    t.resize(t_flat.size());
    for(size_t i=0; i<t_flat.size(); ++i) {
      t[i].resize(t_flat[i].size());
      for(size_t j=0; j<t_flat[i].size(); ++j) {
        t[i][j].make_from_flat(t_flat[i][j]);
      }
    }
  }
// 
//   static void load_and_construct( Archive & ar, cereal::construct<tree_samples> & construct )
//   {
//     size_t m,p,ndraws,basis_dim;
//     xinfo xi;
//     std::vector<std::vector<flat_tree> > t_flat;
//     ar( m, p, ndraws, basis_dim, xi, t_flat );
//     construct(m, p, ndraws, basis_dim, xi, t_flat );
//   }

  
};

RCPP_EXPOSED_CLASS(tree_samples)
  
#endif
