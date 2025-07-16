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
List countbart(arma::vec y_,
               arma::vec offset_,
               List bart_specs,
               List bart_designs,
               arma::mat random_des,
               arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
               double random_var_df, arma::vec randeff_scales,
               int burn, int nd, int thin,      // Draw nd*thin + burn samples, saving nd draws after burn-in
               int count_model,                 // type of count model (1 = poisson, 2 = nb, 3 = zipoisson, 4 = zinb)
               double lambda, double nu,        // prior pars for sigma^2_y
               double kappa_a, double kappa_b,  // prior pars for kappa (shape parameters of beta prime distribution)
               double leaf_c, double leaf_d,    // leaf hyperparameters
               double z_c, double z_d,          // leaf hyperparameters for zero-inflated model
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

  //  std::string treef_name = as<std::string>(treef_name_);
  //  std::ofstream treef(treef_name.c_str());
  //  std::string treef_name_serial = as<std::string>(treef_name_serial_);
  //  std::ofstream treef_serial(treef_name_serial.c_str(), std::ios::out | std::ios::binary);

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

  int n_y0 = -1; // size of 0-block for zero-inflated models (assumes y_ is sorted).

  for(NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2

    if (*it == 0) n_y0 += 1; // increase size of y=0 block
  }

  // Rcout << "n_y0: " << n_y0 << endl;
  
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
  std::vector<int> z ;              // latent vector z (for zero-inflated models)
  std::vector<double> log_w;        // f_1 / (f_0 + f_1) (for zero-inflated models)
  std::vector<double> log_w_denom;  // log(f_0 + f_1) (for zero-inflated models)
  std::vector<double> log_phi;      // latent variable (for zero-inflated models)
  std::vector<double> log_xi;       // latent log(xi) values (for neg binom models)
  std::vector<double> u_vec = y;    // vector used to store sufficient statistics for tree updates 

  for(std::size_t it = 0; it < y.size(); ++it){
    z.push_back(1);
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
    
    // Rcout << "desi " << i << endl;
    
    // the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
    // std::vector<double> x_con;
    
    bool g = dtemp["group"];
    group[i] = g;
    
    Rcout << "group setting: " << g << endl;
    
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
    
    // Rcout << "a " << i << endl;
    
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
    
    // Rcout << "b " << i << endl;
    
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

  // store 'sufficient' statistic (v_i in Murray's paper)
  std::vector<double> r_tree(n);
  std::fill(r_tree.begin(), r_tree.end(), 0.0);

  // log(mean)
  std::vector<double> log_mean(n);
  std::fill(log_mean.begin(), log_mean.end(), 0.0);

  Rcout << "Number of forests: " << num_forests <<  "\n\n" << endl;
  
  std::vector<dinfo> di(num_forests);
  std::vector<std::vector<std::vector<tree::tree_cp> > > node_pointers(num_forests);
  std::vector<double> sample_eta(num_forests);
  
  std::vector<std::stringstream> tree_streams(num_forests);
  std::vector<std::stringstream> serial_streams(num_forests);
  std::vector<Rcpp::CharacterVector> Rtree_streams(num_forests);
  
  std::vector<tree_samples> final_tree_trace(num_forests);

  // contains the fit from a single tree
  double* ftemp  = new double[n]; 
  
  for(size_t i=0; i<num_forests; ++i) {
  
    //Rcout << i << endl;
    List spec = bart_specs[i];
        
    double es = spec["sample_eta"];
    sample_eta[i] = es;
        
    int desi = spec["design_index"];
    size_t ntree = spec["ntree"];
    trees[i].resize(ntree);
    prior_info[i].vanilla = spec["vanilla"];
        
    //tree_streams[i].precision(10);
    /*
    tree_streams[i] << x_info[desi] << endl;
    tree_streams[i] << ntree << endl;
    tree_streams[i] << covariate_dim[desi] << endl;
    tree_streams[i] << Omega[desi].n_rows << endl;
    tree_streams[i] << nd << endl;
    */
        
    /*
    //save stuff to tree file
    treef << xi << endl; //cutpoints
    treef << m << endl;  //number of trees
    treef << p << endl;  //dimension of x's
    treef << (int)(nd/thin) << endl;
    */
    //Rcout << "a" << endl;
    
    for(size_t j=0; j<ntree; ++j) trees[i][j].setm(zeros(Omega[desi].n_rows));

    prior_info[i].pbd = 1.0;  // prob of birth/death move
    prior_info[i].pb = .5;    // prob of birth given birth/death

    if ((count_model == 3) || (count_model == 4)){
      if ((num_forests - i) < 3){
        // zero-inflated tree
        prior_info[i].c = z_c;
        prior_info[i].d = z_d;
      } else {
        // log-linear tree
        prior_info[i].c = leaf_c;
        prior_info[i].d = leaf_d;
      }
    } else {
      // log-linear tree
      prior_info[i].c = leaf_c;
      prior_info[i].d = leaf_d;
    }

    prior_info[i].alpha = spec["alpha"]; // prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
    prior_info[i].beta  = spec["beta"];
    prior_info[i].sigma = shat;

    prior_info[i].mu0 = as<arma::vec>(spec["mu0"]);

    // Rcout << "b" << endl;
    prior_info[i].Sigma0 = as<arma::mat>(spec["Sigma0"]);
    prior_info[i].Prec0 = prior_info[i].Sigma0.i();
    prior_info[i].logdetSigma0 = log(det(prior_info[i].Sigma0));
    prior_info[i].eta = 1;
    prior_info[i].gamma = 1;
    prior_info[i].scale_df = spec["scale_df"];

    // Rcout << "c" << endl;
    // Rcout << "d" << endl;
    // Rcout << desi << endl;
    
    // data info
    dinfo dtemp;
    dtemp.n=n;
    // Rcout << "d1" << endl;
    dtemp.p = covariate_dim[desi];
    // Rcout << "d11" << endl;
    dtemp.x = &(x[desi])[0];
    
    // Rcout << "ptr test " << &(x[desi])[0] << " " <<  &x[desi][0] << endl;
    // Rcout << "ptr test " << *(&x[desi][0]) << " " <<  *(&x[desi][3]) << endl;
    // Rcout << "d12" << endl;
    dtemp.y = &r_tree[0]; // the y for each draw will be the residual

    // adding observations (y)
    dtemp.u_i = &u_vec[0];

    // log(offset) vector
    dtemp.offset = &offset[0];

    // Rcout << "d2" << endl;
    dtemp.basis_dim = Omega[desi].n_rows;
    dtemp.omega = &(Omega[desi])[0];
    
    // Rcout<< "groups" << endl;
    dtemp.groups = &(groups[desi])[0];
    // Rcout<< "group" << endl;
    dtemp.group = group[desi];
    
    // Rcout << "e" << endl;
    
    // Initialize node pointers & allfits
    node_pointers[i].resize(ntree);
    allfits[i].resize(n);
    std::fill(allfits[i].begin(), allfits[i].end(), 0.0);


    for(size_t j=0; j<ntree; ++j) {
      
      node_pointers[i][j].resize(n);
      
      // update ftemp with j'th tree in i'th forest
      fit_loglinear(trees[i][j], x_info[desi], dtemp, ftemp, node_pointers[i][j], true, prior_info[i].vanilla);
      
      // evaluates \sum_{j = 1}^{m} g(x_k; T_j, L_j) for all trees i, units k
      for(size_t k=0; k<n; ++k) allfits[i][k] += ftemp[k];
      
      // Rcout << "allfits test" << allfits[i][1] << " " << allfits[i][2] << endl;
    }

    // DART
    prior_info[i].dart = spec["dart"];
    if(prior_info[i].dart) prior_info[i].dart_alpha = 1.0;
    std::vector<double> vp_mod(covariate_dim[desi], 1.0/covariate_dim[desi]);
    prior_info[i].var_probs = vp_mod;

    // todo: var sizes adjustment
    // //DART
    // if(dart) {
    //   pi_con.dart_alpha = 1;
    //   pi_mod.dart_alpha = 1;
    //   if(var_sizes_con.size() < di_con.p) {
    //     pi_con.var_sizes.resize(di_con.p);
    //     std::fill(pi_con.var_sizes.begin(),pi_con.var_sizes.end(), 1.0/di_con.p);
    //   }
    //   if(var_sizes_mod.size() < di_mod.p) {
    //     pi_mod.var_sizes.resize(di_mod.p);
    //     std::fill(pi_mod.var_sizes.begin(),pi_mod.var_sizes.end(), 1.0/di_mod.p);
    //   }
    // }
    
    di[i] = dtemp;
    
    tree_samples ts(ntree, di[i].p, nd, di[i].basis_dim, x_info[desi]);
    final_tree_trace[i] = ts;
  }

  Rcout << "Done." << endl;

  //-------------------------------------------------------------------------//
  // setup random effects (not used in count models)                         //
  //-------------------------------------------------------------------------//
  size_t random_dim = random_des.n_cols;
  int nr = 1;
  if(randeff) nr = n;

  arma::vec r(nr); // working residuals
  arma::vec Wtr(random_dim); // W'r
        
  arma::mat WtW = random_des.t()*random_des; //W'W
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  // PX parameters
  arma::vec eta(random_var_ix.n_cols); // random_var_ix is num random effects by num variance components
  eta.fill(1.0);

  for(size_t k=0; k<nr; ++k) {
    r(k) = y[k];
    for(size_t j=0; j<num_forests; ++j) {
      r(k) -= allfits[j][k];
    }
  }

  Wtr = random_des.t()*r;
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
  arma::vec allfit_random = random_des*gamma;
  if(!randeff) allfit_random.fill(0);

  // allfit_random.fill(0);


  //-------------------------------------------------------------------------//
  // initialize tree fits                                                    //
  //-------------------------------------------------------------------------//
  double* allfit = new double[n]; // yhat

  for(size_t i=0;i<n;i++) {
    allfit[i] = 0;
    for(size_t j=0; j<num_forests; ++j) {
      allfit[i] += allfits[j][i];
    }
    if(randeff) allfit[i] += allfit_random[i];
  }

  // initialize log_mean 
  for (size_t k = 0; k < n; k++){
    log_mean[k] = offset[k] + allfit[k];
  }

  // output storage
  NumericVector sigma_post(nd);
  NumericVector kappa_post(nd);

  // if(R_trace) forest_trace_R[s][save_ctr][j] = trees[s][j].flatten();
  std::vector<std::vector<List> > forest_trace_R(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    forest_trace_R[j].resize(nd);
    for(size_t i=0; i<nd; ++i) {
      List init_list(trees[j].size());
      forest_trace_R[j][i] = init_list;
    }
  }

  std::vector<NumericMatrix> forest_fits(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    NumericMatrix postfits(nd,n);
    forest_fits[j] = postfits;
  }
  
  // NumericMatrix m_post(nd,n);
  NumericMatrix yhat_post(nd, n);
  // NumericMatrix b_post(nd,n);

  NumericMatrix etas_post(nd, num_forests);
  
  // TODO: return DART stuff
  //  NumericMatrix var_prob_con(nd,pi_con.var_probs.size());
  //  NumericMatrix var_prob_mod(nd,pi_mod.var_probs.size());

  //  NumericMatrix m_est_post(nd,n_con_est);
  //  NumericMatrix b_est_post(nd,n_mod_est);

  arma::mat gamma_post(nd, gamma.n_elem);
  arma::mat random_sd_post(nd, random_var.n_elem);

  std::vector<arma::cube> post_coefs(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    arma::cube tt(di[j].basis_dim, di[j].n, nd);
    tt.fill(0);
    post_coefs[j] = tt;
  }
   
  //   arma::cube scoefs_mod(di_mod.basis_dim, di_mod.n, nd);
  //   arma::mat coefs_mod(di_mod.basis_dim, di_mod.n);
  // 
  //   arma::cube scoefs_con(di_con.basis_dim, di_con.n, nd);
  //   arma::mat coefs_con(di_con.basis_dim, di_con.n);

  //  NumericMatrix spred2(nd,dip.n);

  /*
  //save stuff to tree file
  treef << xi << endl; //cutpoints
  treef << m << endl;  //number of trees
  treef << p << endl;  //dimension of x's
  treef << (int)(nd/thin) << endl;
  */

  //-------------------------------------------------------------------------//
  // MCMC                                                                    //
  //-------------------------------------------------------------------------//
          
  Rcout << "\nBeginning MCMC:\n";
  time_t tp;
  int time1 = time(&tp);      
  size_t save_ctr = 0;

  for(size_t i = 0; i < (nd*thin+burn); i++) {

    if(prior_sample) {
      for(int k=0; k<n; k++) y[k] = gen.normal(allfit[k], sigma);
    }
          
    if(lower_bd.size()>1) {
      for(int k=0; k<n; k++) {
        // Rcout << y[k] << " " << allfit[k] << " " << sigma << " "<< lower_bd[k] << endl;
        if(lower_bd[k]!=-INFINITY) y[k] = rtnormlo(allfit[k], sigma, lower_bd[k]);
        // Rcout << y[k] << endl;
      }
    }
          
    if(upper_bd.size()>1) {
      for(int k=0; k<n; k++) {
        if(upper_bd[k]!=INFINITY) y[k] = -rtnormlo(-allfit[k], sigma, -upper_bd[k]);
      }
    }


    //Rcout << "a" << endl;
    Rcpp::checkUserInterrupt();
    if(i%status_interval == 0) {
      Rcout << "iteration: " << i << endl;

      if ((count_model == 2) || (count_model == 4)){
        Rcout << "dispersion parameter acceptance rate: " << kappa_acpt_rate << endl;
      }
    }

    //-----------------------------------------------------------------------//
    //  update kappa (dispersion parameter in NB model)                      //
    //-----------------------------------------------------------------------//
    if ((count_model == 2) || (count_model == 4)){

      // propose new kappa: random walk on log scale w/ Gaussian proposal        
      double kappa_star = exp(gen.normal(log(kappa), kappa_prop_sd));

      // MH acceptance ratio: 
      //      log-lik + beta' prior log-density + log-scale proposal correction
      double log_a_num = ll_loglinear(y, log_mean, kappa_star, count_model, true, log_w, n_y0)  
        + R::dbeta(kappa_star / (1 + kappa_star), kappa_a, kappa_b, true) - 2 * log1p(kappa_star)
        + log(kappa_star);
      double log_a_denom = ll_loglinear(y, log_mean, kappa, count_model, true, log_w, n_y0)
        + R::dbeta(kappa / (1 + kappa), kappa_a, kappa_b, true) - 2 * log1p(kappa)
        + log(kappa);

      if(log(gen.uniform()) < log_a_num - log_a_denom){
        // accept
        kappa = kappa_star;
        kappa_acpt_rate = ((kappa_acpt_rate * i) + 1) / (i + 1);
      }
      else
      {
        // reject
        kappa_acpt_rate = ((kappa_acpt_rate * i) + 0) / (i + 1);
      }

    }

    //-----------------------------------------------------------------------//
    //  update Z_i (zero-inflated models)                                    //
    //-----------------------------------------------------------------------//
    if ((count_model == 3) || (count_model == 4)){
      drz_loglinear(z, log_w, kappa, log_mean, n_y0, count_model, gen);
    }

    //-----------------------------------------------------------------------//
    //  update xi (negative binomial)                                        //
    //-----------------------------------------------------------------------// 
    if ((count_model == 2) || (count_model == 4)){
      drxi_loglinear(log_xi, kappa, log_mean, y, z, gen);
    }


    //-----------------------------------------------------------------------//
    //  update phi (zero-inflated)                                           //
    //-----------------------------------------------------------------------// 
    if ((count_model == 3) || (count_model == 4)){
      drphi_loglinear(log_phi, log_w_denom, gen);
    }

    //-----------------------------------------------------------------------//
    // update trees                                                          //
    //-----------------------------------------------------------------------//
    for(size_t s = 0; s < num_forests; ++s) {           // s: index over forests
      for(size_t j = 0; j < trees[s].size(); ++j) {     // j: index over trees in f_s

        // gives the fit (updated in 'ftemp') for tree j in forest s
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        for (size_t k = 0; k < n; k++){

          // checks that fit value is valid
          if (ftemp[k] != ftemp[k]){
            // Rcout << t_con[j] << endl;
            stop("nan in ftemp");
          }

          // remove fit of tree considered for update
          allfit[k] = allfit[k] - prior_info[s].eta * ftemp[k];
          allfits[s][k] = allfits[s][k] - prior_info[s].eta * ftemp[k];

          if (((count_model == 3) || (count_model == 4)) & ((num_forests - s) < 3)){

            // zero inflated tree updates

            // determine which tree is being updated
            int tree_ind = num_forests - s - 2;
            int zi_tree = std::abs(tree_ind); // 0 if updating f_0, 1 if updating tree 1

            // sufficient stats
            int zero_type = z[k] - zi_tree;
            u_vec[k] = 1 - std::abs(zero_type);           // u_i in paper
            r_tree[k] = exp(log_phi[k] + allfits[s][k]);  // f_h(x_i) * v_i in paper

          } else {
            
            // update mean structure
            log_mean[k] = log_mean[k] - prior_info[s].eta * ftemp[k];
            // update sufficient statistic vectors
            u_vec[k] = z[k] * y[k];                           // u_i in paper
            r_tree[k] = z[k] * exp(log_xi[k] + log_mean[k]);  // f_h(x_i) * v_i in paper

          }

          if (r_tree[k] != r_tree[k]){
            // Rcout << (y[k]-allfit[k]) << endl;
            // Rcout << pi_con.eta << endl;
            // Rcout << r_con[k] << endl;
            stop("NaN in resid");
          }
        }

        // Rcout << " bd " << endl;
        bd_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen, node_pointers[s][j]);
        
        // Rcout << " drmu" << endl;
        drmu_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen);

        // Rcout << " second fit" << endl;
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        // add tree info back into allfits object
        for(size_t k=0;k<n;k++) {
          if (((count_model == 3) || (count_model == 4)) & ((num_forests - s) < 3)){
            allfit[k] += prior_info[s].eta*ftemp[k];
            allfits[s][k] += prior_info[s].eta*ftemp[k];
          } else {
            allfit[k] += prior_info[s].eta*ftemp[k];
            allfits[s][k] += prior_info[s].eta*ftemp[k];
            log_mean[k] += prior_info[s].eta * ftemp[k];
          }
        }
      }
    }

    //-----------------------------------------------------------------------//
    // update  w(x_i)                                                        //
    //-----------------------------------------------------------------------//
    if ((count_model == 3) || (count_model == 4)){
      for (size_t k = 0; k < n; k++){
        log_w_denom[k] = logsumexp(allfits[num_forests - 2][k], allfits[num_forests - 1][k]);
        log_w[k] = allfits[num_forests - 1][k] - log_w_denom[k];
      }
    }

    /*
      //  update PX parameters (NOT USED FOR COUNT MODELS, set to 'false' below)
      double eta_old;
      if(false) {
        for(size_t s=0; s<num_forests; ++s) {
          if(sample_eta[s]>0) {
            for(size_t k=0;k<n;k++) {
              ftemp[k] = y[k] - (allfit[k] - allfits[s][k]);
            }
            eta_old = prior_info[s].eta;
            //update_scale(ftemp, &(allfits[s])[0], n, sigma, prior_info[s], gen); <- seems to work
            update_scale(ftemp, allfits[s], n, sigma, prior_info[s], gen);

            //Rcout << "s = " << s << " gamma = " << prior_info[s].gamma << " eta = " << prior_info[s].eta << endl;

            for(size_t k=0; k<n; ++k) {
              allfit[k] -= allfits[s][k];
              allfits[s][k] = allfits[s][k] * prior_info[s].eta / eta_old;
              allfit[k] += allfits[s][k];
            }

            prior_info[s].sigma = sigma/fabs(prior_info[s].eta);
          }
        }
      }
    */

    //-----------------------------------------------------------------------//
    // update  random effects (not used for count models)                    //
    //-----------------------------------------------------------------------//
    if(randeff) {
      
      // update random effects
      for(size_t k=0; k<n; ++k) {
        r(k) = y[k] - allfit[k] + allfit_random[k];
        allfit[k] -= allfit_random[k];
      }

      Wtr = random_des.t()*r;

      arma::mat adj = diagmat(random_var_ix*eta);
      // Rcout << adj << endl << endl;
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
      Phi = 0.5*(Phi + Phi.t());
      arma::vec m = adj*Wtr/(sigma*sigma);
      // Rcout << m << Phi << endl << Sigma_inv_random;
      gamma = rmvnorm_post(m, Phi);

      // Rcout << "updated gamma";

      // Update px parameters eta
      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
      Phi2 = 0.5*(Phi2 + Phi2.t());
      eta = rmvnorm_post(m2, Phi2);

      // Rcout << "updated eta";

      // Update variance parameters
      arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
      // Rcout << "A";
      arma::rowvec counts = sum(random_var_ix, 0);
      // Rcout << "B";
      for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df/randeff_scales(ii)*randeff_scales(ii) + ssqs(ii));
      }
      // Rcout << "updated vars" << endl;
      Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

      // Rcout << random_var_ix*random_var;
      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;

      // is rebuilding allfits still necessary?
      for(size_t k=0; k<n; ++k) {
        allfit[k] = allfit_random(k);
        for(size_t s=0; s<num_forests; ++s) {
          allfit[k] += allfits[s][k];
        }
        // allfit[k] = allfit_con[k] + allfit_mod[k] + ; //+= allfit_random[k];
      }
    }
  
    /*
    //draw sigma
    double rss = 0.0;
    double restemp = 0.0;
    for(size_t k=0;k<n;k++) {
      restemp = y[k]-allfit[k];
      rss += restemp*restemp;
    }
    //Rcout << y[0] << " " << y[5] << endl;
    //Rcout << allfit[0] << " " << allfit[5] << endl;
    //Rcout << "rss " << rss << endl;
    if(!probit) sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    //pi_con.sigma = sigma/fabs(pi_con.eta);
    //pi_mod.sigma = sigma/fabs(pi_mod.eta);
      
    for(size_t s=0; s<num_forests; ++s) {
      // Rcout << "sigma " << sigma << " eta " <<prior_info[s].eta << endl;
      prior_info[s].sigma = sigma/fabs(prior_info[s].eta);
    }
    */ 


    //-----------------------------------------------------------------------//
    // save results                                                          //
    //-----------------------------------------------------------------------//
    if( ((i>=burn) & (i % thin==0))) {
      
      // for(size_t j=0;j<m;j++) treef << t[j] << endl;
        
      //      msd_post(save_ctr) = fabs(pi_con.eta)*con_sd;
      //      bsd_post(save_ctr) = fabs(pi_mod.eta)*mod_sd;
        
      //pi_mod.var_probs
      
      // for(size_t j=0; j<pi_con.var_probs.size(); ++j) {
      //   var_prob_con(save_ctr, j) = pi_con.var_probs[j];
      // }
      // for(size_t j=0; j<pi_mod.var_probs.size(); ++j) {
      //   var_prob_mod(save_ctr, j) = pi_mod.var_probs[j];
      // }

      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_sd_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();
        
      sigma_post(save_ctr) = sigma;
      kappa_post(save_ctr) = kappa;
        
      // eta_con_post(save_ctr) = pi_con.eta;
      // eta_mod_post(save_ctr) = pi_mod.eta;
        
      for(size_t k=0; k<n; k++) {
        //        m_post(save_ctr, k) = allfit_con[k];
        //        b_post(save_ctr, k) = allfit_mod[k];
        yhat_post(save_ctr, k) = log_mean[k];
      }
      
      for(size_t s=0; s<num_forests; ++s) {
        etas_post(save_ctr,s) = prior_info[s].eta;
        for(size_t j=0; j< trees[s].size(); ++j) { 
          post_coefs[s].slice(save_ctr) += prior_info[s].eta*coef_basis(trees[s][j], x_info[s], di[s]);
          // if(text_trace) tree_streams[s] << trees[s][j];
          final_tree_trace[s].t[save_ctr][j] = trees[s][j];
          final_tree_trace[s].t[save_ctr][j].compress();
          final_tree_trace[s].t[save_ctr][j].scale(prior_info[s].eta);
          // if(R_trace) forest_trace_R[s][save_ctr][j] = trees[s][j].flatten(prior_info[s].eta);
        } 
      }
  
      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;
  
  delete[] allfit;
  delete[] ftemp;

  std::vector<Rcpp::RawVector> Rtree_serial_streams(num_forests);
      
  //   std::stringstream ss;
  //   {
  //     cereal::BinaryOutputArchive oarchive(ss); // Create an output archive
  //     oarchive(my_instance);
  //   }
  //   ss.seekg(0, ss.end);
  //   RawVector retval(ss.tellg());
  //   ss.seekg(0, ss.beg);
  //   ss.read(reinterpret_cast<char*>(&retval[0]), retval.size());
  //   return retval;
  // }
  
  for(size_t s=0; s<num_forests; ++s) {
    Rtree_streams[s] = final_tree_trace[s].save_string(); // tree_streams[s].str();
    {
      cereal::BinaryOutputArchive oarchive(serial_streams[s]); // Create an output archive
      oarchive(final_tree_trace[s]); // Write the data to the archive
    }
    // Rtree_serial_streams[s] = serial_streams[s].str();
    serial_streams[s].seekg(0, serial_streams[s].end);
    RawVector retval(serial_streams[s].tellg());
    serial_streams[s].seekg(0, serial_streams[s].beg);
    serial_streams[s].read(reinterpret_cast<char*>(&retval[0]), retval.size());
    Rtree_serial_streams[s] = retval;
  }
  
  return(List::create(_["yhat_post"] = yhat_post,
                      _["coefs"] = post_coefs,
                      _["etas"] = etas_post,
                      _["sigma"] = sigma_post, 
                      _["kappa"] = kappa_post,
                      _["kappa_acceptance"] = kappa_acpt_rate,
                      _["gamma"] = gamma_post, 
                      _["random_sd_post"] = random_sd_post,
                      _["tree_streams"] = Rtree_streams,
                      _["tree_serials"] = Rtree_serial_streams,
                      _["tree_trace"] = final_tree_trace,
                      _["y_last"] = y
  ));
}
