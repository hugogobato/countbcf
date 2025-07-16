#ifndef slice_h
#define slice_h

#include <vector>
#include "funs.h"
#include "info.h"

class logdensity {
  public:
  virtual double val(double y) = 0;
};

class ld_norm: public logdensity {
  public:
  double mu;
  double sigma;
  double val(double y) { return(R::dnorm(y, mu, sigma, 1)); }
  ld_norm(double mu_, double sigma_) { mu=mu_; sigma=sigma_; }
};

class ld_dart_alpha: public logdensity {
  public:
  std::vector<double> counts;
  std::vector<double> var_sizes;
  std::vector<double> probs;
  double dart_a, dart_b, dart_rho;
  ld_dart_alpha(std::vector<double> counts_, std::vector<double> var_sizes_, std::vector<double> probs_) { 
    ld_dart_alpha(counts_, var_sizes_, probs_, dart_a = 0.5, dart_b = 1, dart_rho = counts_.size());
  }
  ld_dart_alpha(std::vector<double> counts_, std::vector<double> var_sizes_, std::vector<double> probs_,
                double dart_a_, double dart_b_, double dart_rho_) { 
    counts = counts_;
    var_sizes = var_sizes_;
    probs  = probs_;
    dart_a = dart_a_; 
    dart_b = dart_b_;
    dart_rho = dart_rho_;
  }
  double val(double y) { 
    double ret = R::dbeta(y/(y+dart_rho), dart_a, dart_b, 1) + log(dart_rho) - 2*log(y+dart_rho); //log prior
    double csum = 0.0;
    for(size_t i=0; i<counts.size(); ++i) {
      double totct = counts[i] + y/(counts.size()*var_sizes[i]);
      csum += totct;
      ret += (totct-1)*log(probs[i]) - lgamma(totct);
    }
    ret += lgamma(csum);
    return(ret);
  }
};


double slice(double x0, logdensity* g, double w=1., double m=INFINITY, 
             double lower=-INFINITY, double upper=INFINITY);
#endif slice_h
