#include "arma_config.h"
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "tree_samples.h"

using namespace Rcpp;


// [[Rcpp::export]]
List countbcfgemini(arma::vec y_,
               arma::vec offset_,
               List bart_specs,
               List bart_designs, // X's and Z's (in Omega) are found here
               arma::mat random_des,
               arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
               double random_var_df, arma::vec randeff_scales,
               int burn, int nd, int thin,      // Draw nd*thin + burn samples, saving nd draws after burn-in
               int count_model,                 // type of count model (1 = poisson, 2 = nb, 3 = zipoisson, 4 = zinb)
               double lambda, double nu,        // prior pars for sigma^2_y
               double kappa_a, double kappa_b,  // prior pars for kappa (shape parameters of beta prime distribution)
               double leaf_c, double leaf_d,    // leaf hyperparameters for f
               double z_c, double z_d,          // leaf hyperparameters for f0, f1
               double kappa_prop_sd = 0.2,      // standard deviation of kappa MH proposal distribution
               bool return_trees = true,
               bool save_trees = false,
               bool est_mod_fits = false, bool est_con_fits = false,
               bool prior_sample = false,
               int status_interval = 100,
               NumericVector lower_bd = NumericVector::create(0.0),
               NumericVector upper_bd = NumericVector::create(0.0),
               bool probit = false,
               bool text_trace = true,
               bool R_trace = false)
{

  // check if random effects have been specified
  bool randeff = true;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  }
  if(randeff) Rcout << "Using random effects." << std::endl;

  // set up random number generator, used for all draws
  RNGScope scope;
  RNG gen; 

  //-------------------------------------------------------------------------//
  // Read / format 'y'                                                       //
  //-------------------------------------------------------------------------//
  Rcout << "Reading in y\n\n";

  std::vector<double> y; // storage for y
  std::vector<double> offset; // storage for log(offset) vector

  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys; // sufficient stats for all of y, use to initialize the bart trees.

  int n_y0 = 0; // BUG FIX: Correctly initialized to 0.
  for(NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2

    if (*it == 0) n_y0++;
  }

  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n; // sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); // sample standard deviation
  if(probit) shat = 1.0;
  
  double sigma = shat;

  // log offset vector
  for(NumericVector::iterator it = offset_.begin(); it != offset_.end(); ++it){
    offset.push_back(*it);
  }

  // kappa (dispersion parameter for negative binomial models; ignored for all other models)
  double kappa = 1;
  double kappa_acpt_rate = 0;

  // initialize latent parameters
  std::vector<int> z_latent;        // latent vector z (for zero-inflated models)
  std::vector<double> log_w;        // log(w) for ZI model
  std::vector<double> log_w_denom;  // log(f_0 + f_1) for ZI model
  std::vector<double> log_phi;      // latent variable (for ZI models)
  std::vector<double> log_xi;       // latent log(xi) values (for NB models)
  std::vector<double> u_vec = y;    // vector used to store sufficient statistics for tree updates 

  for(std::size_t it = 0; it < y.size(); ++it){
    z_latent.push_back(1);
    log_w.push_back(log(0.5));
    log_w_denom.push_back(log(2));
    log_phi.push_back(0);
    log_xi.push_back(0);
  }


  //-------------------------------------------------------------------------//
  // Read, format design info                                                //
  //-------------------------------------------------------------------------//
  Rcout << "Setting up designs\n\n";
  size_t num_designs = bart_designs.size();
  std::vector<std::vector<double> > x(num_designs);
  std::vector<std::vector<int> > groups(num_designs);
  std::vector<bool> group(num_designs);
  std::vector<xinfo> x_info(num_designs);
  std::vector<arma::mat> Omega(num_designs);
  std::vector<size_t> covariate_dim(num_designs);

  for(size_t i=0; i<num_designs; i++) {
    Rcout << "design " << i << endl;
    List dtemp = bart_designs[i];
    bool g = dtemp["group"];
    group[i] = g;
    IntegerVector gt_ = dtemp["groups"];
    for(IntegerVector::iterator it=gt_.begin(); it!= gt_.end(); ++it) {
      groups[i].push_back(*it);
    }
    NumericVector xt_ = dtemp["X"];
    for(NumericVector::iterator it=xt_.begin(); it!= xt_.end(); ++it) {
      x[i].push_back(*it);
    }
    size_t p = x[i].size()/n;
    covariate_dim[i] = p;
    Rcout << "Instantiated covariate matrix " << i+1 << " with " << p << " columns" << endl;
    xinfo xi;
    xi.resize(p);
    List x_info_list = dtemp["info"];
    for(int j=0; j<p; ++j) {
      NumericVector tmp = x_info_list[j];
      std::vector<double> tmp2;
      for(size_t s=0; s<tmp.size(); ++s) {
        tmp2.push_back(tmp[s]);
      }
      xi[j] = tmp2;
    }
    x_info[i] = xi;
    Omega[i] = as<arma::mat>(dtemp["Omega"]);
  }

  //-------------------------------------------------------------------------//
  // Set up forests                                                          //
  //-------------------------------------------------------------------------//
  Rcout << "Setting up forests\n\n";

  size_t num_forests = bart_specs.size();
  std::vector<std::vector<tree> > trees(num_forests);
  std::vector<pinfo> prior_info(num_forests);
  std::vector<std::vector<double> > allfits(num_forests);
  // ROBUSTNESS IMPROVEMENT: Explicitly store component and type for each forest
  std::vector<std::string> forest_components(num_forests);
  std::vector<std::string> forest_types(num_forests);

  std::vector<double> r_tree(n, 0.0);
  std::vector<double> log_mean(n, 0.0);

  Rcout << "Number of forests: " << num_forests <<  "\n\n" << endl;
  
  std::vector<dinfo> di(num_forests);
  std::vector<std::vector<std::vector<tree::tree_cp> > > node_pointers(num_forests);
  std::vector<double> sample_eta(num_forests);
  
  std::vector<tree_samples> final_tree_trace(num_forests);
  double* ftemp  = new double[n]; 
  
  for(size_t i=0; i<num_forests; ++i) {
    List spec = bart_specs[i];
    
    // Store component and type for clarity and robustness
    forest_components[i] = as<std::string>(spec["component"]);
    forest_types[i] = as<std::string>(spec["type"]);
    Rcout << "Forest " << i << ": component='" << forest_components[i] << "', type='" << forest_types[i] << "'" << endl;
        
    sample_eta[i] = spec["sample_eta"];
    int desi = spec["design_index"];
    size_t ntree = spec["ntree"];
    trees[i].resize(ntree);
    prior_info[i].vanilla = spec["vanilla"];
    
    for(size_t j=0; j<ntree; ++j) trees[i][j].setm(zeros(Omega[desi].n_rows));

    prior_info[i].pbd = 1.0;
    prior_info[i].pb = .5;

    // Assign leaf hyperparameters based on the explicit component name
    if (forest_components[i] == "f") {
      prior_info[i].c = leaf_c;
      prior_info[i].d = leaf_d;
    } else if (forest_components[i] == "f0" || forest_components[i] == "f1") {
      prior_info[i].c = z_c;
      prior_info[i].d = z_d;
    } else {
      stop("Unknown forest component specified: " + forest_components[i]);
    }

    prior_info[i].alpha = spec["alpha"];
    prior_info[i].beta  = spec["beta"];
    prior_info[i].sigma = shat;
    prior_info[i].mu0 = as<arma::vec>(spec["mu0"]);
    prior_info[i].Sigma0 = as<arma::mat>(spec["Sigma0"]);
    prior_info[i].Prec0 = prior_info[i].Sigma0.i();
    prior_info[i].logdetSigma0 = log(det(prior_info[i].Sigma0));
    prior_info[i].eta = 1;
    prior_info[i].gamma = 1;
    prior_info[i].scale_df = spec["scale_df"];

    dinfo dtemp;
    dtemp.n=n;
    dtemp.p = covariate_dim[desi];
    dtemp.x = &(x[desi])[0];
    dtemp.y = &r_tree[0]; 
    dtemp.u_i = &u_vec[0];
    dtemp.offset = &offset[0];
    dtemp.basis_dim = Omega[desi].n_rows;
    dtemp.omega = &(Omega[desi])[0];
    dtemp.groups = &(groups[desi])[0];
    dtemp.group = group[desi];
    
    node_pointers[i].resize(ntree);
    allfits[i].assign(n, 0.0);
    for(size_t j=0; j<ntree; ++j) {
      node_pointers[i][j].resize(n);
      fit_loglinear(trees[i][j], x_info[desi], dtemp, ftemp, node_pointers[i][j], true, prior_info[i].vanilla);
      for(size_t k=0; k<n; ++k) allfits[i][k] += ftemp[k];
    }

    prior_info[i].dart = spec["dart"];
    if(prior_info[i].dart) prior_info[i].dart_alpha = 1.0;
    std::vector<double> vp_mod(covariate_dim[desi], 1.0/covariate_dim[desi]);
    prior_info[i].var_probs = vp_mod;
    
    di[i] = dtemp;
    
    tree_samples ts(ntree, di[i].p, nd, di[i].basis_dim, x_info[desi]);
    final_tree_trace[i] = ts;
  }

  Rcout << "Setup done." << endl;

  //-------------------------------------------------------------------------//
  // initialize fits                                                    //
  //-------------------------------------------------------------------------//
  bool is_zi_model = (count_model == 3) || (count_model == 4);
  
  // Calculate initial combined fits
  for (size_t k = 0; k < n; k++){
    log_mean[k] = offset[k];
    for(size_t s = 0; s < num_forests; ++s) {
        if(forest_components[s] == "f") {
            log_mean[k] += allfits[s][k];
        }
    }
  }

  // output storage setup
  NumericVector kappa_post(nd);
  NumericMatrix yhat_post(nd, n);
  NumericMatrix etas_post(nd, num_forests);
  std::vector<NumericMatrix> forest_fits_post(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    forest_fits_post[j] = NumericMatrix(nd,n);
  }
  std::vector<arma::cube> post_coefs(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    post_coefs[j] = arma::cube(di[j].basis_dim, di[j].n, nd, arma::fill::zeros);
  }
   
  //-------------------------------------------------------------------------//
  // MCMC                                                                    //
  //-------------------------------------------------------------------------//
  Rcout << "\nBeginning MCMC:\n";
  time_t tp;
  int time1 = time(&tp);      
  size_t save_ctr = 0;

  for(size_t i = 0; i < (nd*thin+burn); i++) {
    Rcpp::checkUserInterrupt();
    if(i%status_interval == 0) {
      Rcout << "iteration: " << i << " of " << (nd*thin+burn) << endl;
      if (is_zi_model && (count_model==2 || count_model==4)) Rcout << "  kappa acceptance rate: " << kappa_acpt_rate << endl;
    }

    if ((count_model == 2) || (count_model == 4)){
      double kappa_star = exp(gen.normal(log(kappa), kappa_prop_sd));
      double log_a_num = ll_loglinear(y, log_mean, kappa_star, count_model, true, log_w, n_y0) + R::dbeta(kappa_star / (1 + kappa_star), kappa_a, kappa_b, true) - 2 * log1p(kappa_star) + log(kappa_star);
      double log_a_denom = ll_loglinear(y, log_mean, kappa, count_model, true, log_w, n_y0) + R::dbeta(kappa / (1 + kappa), kappa_a, kappa_b, true) - 2 * log1p(kappa) + log(kappa);
      if(log(gen.uniform()) < log_a_num - log_a_denom){
        kappa = kappa_star;
        kappa_acpt_rate = (kappa_acpt_rate * i + 1.0) / (i + 1.0);
      } else {
        kappa_acpt_rate = (kappa_acpt_rate * i) / (i + 1.0);
      }
    }

    if (is_zi_model){
      drz_loglinear(z_latent, log_w, kappa, log_mean, n_y0, count_model, gen);
    }
 
    if ((count_model == 2) || (count_model == 4)){
      drxi_loglinear(log_xi, kappa, log_mean, y, z_latent, gen);
    }

    if (is_zi_model){
      drphi_loglinear(log_phi, log_w_denom, gen);
    }

    // Update trees
    for(size_t s = 0; s < num_forests; ++s) {
      for(size_t j = 0; j < trees[s].size(); ++j) {
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        // Calculate residuals and sufficient statistics
        for (size_t k = 0; k < n; k++){
          if (ftemp[k] != ftemp[k]) stop("nan in ftemp");

          allfits[s][k] -= ftemp[k]; // No eta scaling for log-linear model

          if (forest_components[s] == "f") {
              log_mean[k] -= ftemp[k];
              u_vec[k] = z_latent[k] * y[k];
              r_tree[k] = z_latent[k] * exp(log_xi[k] + log_mean[k]);
          } else { // ZI components (f0, f1)
              int zi_tree_role = (forest_components[s] == "f0") ? 0 : 1;
              u_vec[k] = 1 - std::abs(z_latent[k] - zi_tree_role);
              // BUG FIX: Correct sufficient stat for ZI parts. It does not depend on a pseudo-residual.
              r_tree[k] = exp(log_phi[k]);
          }
          if (r_tree[k] != r_tree[k]) stop("NaN in resid");
        }

        bd_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen, node_pointers[s][j]);
        drmu_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen);
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        // Add tree fit back
        for(size_t k=0; k<n; k++) {
          allfits[s][k] += ftemp[k];
          if (forest_components[s] == "f") {
              log_mean[k] += ftemp[k];
          }
        }
      }
    }

    // Update w(x_i) for ZI models after all trees are updated
    if (is_zi_model){
      for (size_t k = 0; k < n; k++){
        double f0_total = 0.0;
        double f1_total = 0.0;
        for(size_t s=0; s<num_forests; ++s) {
            if(forest_components[s] == "f0") f0_total += allfits[s][k];
            if(forest_components[s] == "f1") f1_total += allfits[s][k];
        }
        log_w_denom[k] = logsumexp(f0_total, f1_total);
        log_w[k] = f1_total - log_w_denom[k];
      }
    }

    // Save results
    if( (i>=burn) && ( (i-burn) % thin==0) ) {
      kappa_post(save_ctr) = kappa;
        
      for(size_t k=0; k<n; k++) {
        yhat_post(save_ctr, k) = log_mean[k];
      }
      
      for(size_t s=0; s<num_forests; ++s) {
        etas_post(save_ctr,s) = prior_info[s].eta; // eta is always 1, but saved for consistency
        for(size_t k=0; k<n; ++k) {
            forest_fits_post[s](save_ctr, k) = allfits[s][k];
        }
        for(size_t j=0; j< trees[s].size(); ++j) { 
          post_coefs[s].slice(save_ctr) += coef_basis(trees[s][j], x_info[s], di[s]);
          final_tree_trace[s].t[save_ctr][j] = trees[s][j];
          final_tree_trace[s].t[save_ctr][j].compress();
        } 
      }
  
      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;
  
  delete[] ftemp;

  std::vector<Rcpp::RawVector> Rtree_serial_streams(num_forests);
  std::vector<Rcpp::CharacterVector> Rtree_streams(num_forests);
  
  for(size_t s=0; s<num_forests; ++s) {
    Rtree_streams[s] = final_tree_trace[s].save_string();
    std::stringstream serial_stream;
    {
      cereal::BinaryOutputArchive oarchive(serial_stream);
      oarchive(final_tree_trace[s]);
    }
    serial_stream.seekg(0, serial_stream.end);
    RawVector retval(serial_stream.tellg());
    serial_stream.seekg(0, serial_stream.beg);
    serial_stream.read(reinterpret_cast<char*>(&retval[0]), retval.size());
    Rtree_serial_streams[s] = retval;
  }
  
  return(List::create(_["yhat_post"] = yhat_post,
                      _["forest_fits_post"] = forest_fits_post,
                      _["coefs"] = post_coefs,
                      _["etas"] = etas_post,
                      _["kappa"] = kappa_post,
                      _["kappa_acceptance"] = kappa_acpt_rate,
                      _["tree_streams"] = Rtree_streams,
                      _["tree_serials"] = Rtree_serial_streams,
                      _["tree_trace"] = final_tree_trace,
                      _["y_last"] = y
  ));
}