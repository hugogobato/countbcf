#include "arma_config.h"
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


// [[Rcpp::export]]
List bcfCore(arma::vec y_, arma::mat Omega_con, arma::mat Omega_mod, arma::mat Omega_con_est, arma::mat Omega_mod_est,
                  NumericVector x_con_, NumericVector x_mod_, NumericVector x_con_est_, NumericVector x_mod_est_,
                  List x_con_info_list, List x_mod_info_list,
                  arma::mat random_des, //needs to come in with n rows no matter what?
                  arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
                  double random_var_df, arma::vec randeff_scales,
                  int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
                  int ntree_mod, int ntree_con,
                  double lambda, double nu, //prior pars for sigma^2_y
                  arma::mat Sigma0_con, arma::mat Sigma0_mod,
                  double con_alpha, double con_beta,
                  double mod_alpha, double mod_beta,
                  CharacterVector treef_name_, 
                  bool est_mod_fits = false, bool est_con_fits = false,
                  bool prior_sample = false,
                  bool use_con_scale = true, bool use_mod_scale = true,
                  double con_scale_df = 1, double mod_scale_df = -1,
                  int status_interval=100,
                  bool vanilla = false, 
                  bool dart = false,
                  NumericVector var_sizes_con = NumericVector::create(0.0),
                  NumericVector var_sizes_mod = NumericVector::create(0.0),
                  NumericVector lower_bd = NumericVector::create(0.0),
                  NumericVector upper_bd = NumericVector::create(0.0), bool probit=false)
{

  bool randeff = true;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  }

  if(randeff) Rcout << "Using random effects." << std::endl;

  std::string treef_name = as<std::string>(treef_name_);
  std::ofstream treef(treef_name.c_str());

  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws

  //Rcout << "\n*****Into bart main\n";

  /*****************************************************************************
  /* Read, format y
  *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.

  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
  if(probit) shat = 1.0;

  /*****************************************************************************
  /* Read, format X_con
  *****************************************************************************/
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> x_con;
  for(NumericVector::iterator it=x_con_.begin(); it!= x_con_.end(); ++it) {
    x_con.push_back(*it);
  }
  size_t p_con = x_con.size()/n;

  Rcout << "Using " << p_con << " control variables." << std::endl;

  //x cutpoints
  xinfo xi_con;

  xi_con.resize(p_con);
  for(int i=0; i<p_con; ++i) {
    NumericVector tmp = x_con_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_con[i] = tmp2;
  }

  /*****************************************************************************
  /* Read, format X_mod
  *****************************************************************************/
  std::vector<double> x_mod;
  for(NumericVector::iterator it=x_mod_.begin(); it!= x_mod_.end(); ++it) {
    x_mod.push_back(*it);
  }
  size_t p_mod = x_mod.size()/n;

  Rcout << "Using " << p_mod << " effect moderators." << std::endl;

  //x cutpoints
  xinfo xi_mod;

  xi_mod.resize(p_mod);
  for(int i=0; i<p_mod; ++i) {
    NumericVector tmp = x_mod_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_mod[i] = tmp2;
  }

  /*****************************************************************************
  /* Setup the model
  *****************************************************************************/
  //--------------------------------------------------
  //trees

  
  //Initialize to all zeros for now
  std::vector<tree> t_mod(ntree_mod);
  //arma::vec betahat = solve(Omega_mod*Omega_mod.t() + 0.05*eye(Omega_mod.n_rows,Omega_mod.n_rows), Omega_mod*y_);
  for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(zeros(Omega_mod.n_rows));

  //betahat = solve(Omega_con*Omega_con.t()+ 0.05*eye(Omega_con.n_rows,Omega_con.n_rows), Omega_con*y_);
  std::vector<tree> t_con(ntree_con);
  for(size_t i=0;i<ntree_con;i++) t_con[i].setm(zeros(Omega_con.n_rows));

  //--------------------------------------------------
  //prior parameters


  pinfo pi_mod;
  pi_mod.pbd = 1.0; //prob of birth/death move
  pi_mod.pb = .5; //prob of birth given  birth/death

  pi_mod.alpha = mod_alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi_mod.beta  = mod_beta;
  pi_mod.sigma = shat;

  pi_mod.mu0 = zeros(Omega_mod.n_rows);
  
  pi_mod.Sigma0 = Sigma0_mod;
  pi_mod.Prec0 = pi_mod.Sigma0.i();
  pi_mod.logdetSigma0 = log(det(pi_mod.Sigma0));
  pi_mod.eta = 1;
  pi_mod.gamma = 1;
  pi_mod.scale_df = mod_scale_df;

  // DART
  pi_mod.dart = dart;
  std::vector<double> vp_mod(p_mod, 1.0/p_mod);
  pi_mod.var_probs = vp_mod;

  pinfo pi_con;
  pi_con.pbd = 1.0; //prob of birth/death move
  pi_con.pb = .5; //prob of birth given  birth/death

  pi_con.alpha = con_alpha;
  pi_con.beta  = con_beta;

  pi_con.mu0 = zeros(Omega_con.n_rows);
  pi_con.Sigma0 = Sigma0_con;
  pi_con.Prec0 = pi_con.Sigma0.i();
  pi_con.logdetSigma0 = log(det(pi_con.Sigma0));
  pi_con.eta = 1;
  pi_con.gamma = 1;
  pi_con.scale_df = con_scale_df;

  //DART
  pi_con.dart = dart;
  std::vector<double> vp_con(p_con, 1.0/p_con);
  pi_con.var_probs = vp_con;

  pi_con.sigma = shat;

  double sigma = shat;


  //--------------------------------------------------
  //dinfo for control function m(x)
//  Rcout << "ybar " << ybar << endl;
  double* allfit_con = new double[n]; //sum of fit of all trees
  //TESTING
  for(size_t i=0;i<n;i++) allfit_con[i] = 0.0;// ybar;
  //END TESTING
  double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_con;
  di_con.n=n;
  di_con.p=p_con;
  di_con.x = &x_con[0];
  di_con.y=r_con; //the y for each draw will be the residual
  di_con.basis_dim = Omega_con.n_rows;
  di_con.omega = &Omega_con[0];

  //--------------------------------------------------
  //dinfo for trt effect function b(x)
  double* allfit_mod = new double[n]; //sum of fit of all trees
  //TESTING
  for(size_t i=0;i<n;i++) allfit_mod[i] = 0.0; //(z_[i]*bscale1 + (1-z_[i])*bscale0)*trt_init;
  //END TESTING
  double* r_mod = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_mod;
  di_mod.n=n; di_mod.p=p_mod; di_mod.x = &x_mod[0]; di_mod.y=r_mod; //the y for each draw will be the residual
  di_mod.basis_dim = Omega_mod.n_rows; di_mod.omega = &Omega_mod[0];


  //DART
  if(dart) {
    pi_con.dart_alpha = 1;
    pi_mod.dart_alpha = 1;
    if(var_sizes_con.size() < di_con.p) {
      pi_con.var_sizes.resize(di_con.p);
      std::fill(pi_con.var_sizes.begin(),pi_con.var_sizes.end(), 1.0/di_con.p);
    }
    if(var_sizes_mod.size() < di_mod.p) {
      pi_mod.var_sizes.resize(di_mod.p);
      std::fill(pi_mod.var_sizes.begin(),pi_mod.var_sizes.end(), 1.0/di_mod.p);
    }
  }

  //--------------------------------------------------
  //dinfo and design for trt effect function out of sample
  //x for predictions
  dinfo di_mod_est; //data information for prediction
  std::vector<double> x_mod_est;     //stored like x
  size_t n_mod_est;
  //  if(x_mod_est_.size()) {
  for(NumericVector::iterator it=x_mod_est_.begin(); it!=x_mod_est_.end(); ++it) {
    x_mod_est.push_back(*it);
  }
  n_mod_est = x_mod_est.size()/p_mod;
//  Rcout << "n_mod_est " << n_mod_est << std::endl;
  if(x_mod_est.size() != n_mod_est*p_mod) stop("error, wrong number of elements in effect estimate data set\n");
  //if(n_mod_est)
  di_mod_est.n=n_mod_est;
  di_mod_est.p=p_mod;
  di_mod_est.x = &x_mod_est[0];
  di_mod_est.y=0; //there are no y's!
  di_mod_est.basis_dim = Omega_mod_est.n_rows;
  di_mod_est.omega = &Omega_mod_est[0]; //placeholder, wopn't be touched

  dinfo di_con_est; //data information for prediction
  std::vector<double> x_con_est;     //stored like x
  size_t n_con_est;
  //  if(x_con_est_.size()) {
  for(NumericVector::iterator it=x_con_est_.begin(); it!=x_con_est_.end(); ++it) {
    x_con_est.push_back(*it);
  }
  n_con_est = x_con_est.size()/p_con;
  //Rcout << "n_con_est " << n_con_est << std::endl;
  if(x_con_est.size() != n_con_est*p_con) stop("error, wrong number of elements in effect estimate data set\n");
  //if(n_con_est)
  di_con_est.n=n_con_est;
  di_con_est.p=p_con;
  di_con_est.x = &x_con_est[0];
  di_con_est.y=0; //there are no y's!
  di_con_est.basis_dim = Omega_con_est.n_rows;
  di_con_est.omega = &Omega_con_est[0]; //placeholder, wopn't be touched
  //  }

  //  }
  //--------------------------------------------------
  //storage for ouput

  //--------------------------------------------------
  //setup for random effects
  size_t random_dim = random_des.n_cols;
  int nr=1;
  if(randeff) nr = n;

  arma::vec r(nr); //working residuals
  arma::vec Wtr(random_dim); // W'r

  arma::mat WtW = random_des.t()*random_des; //W'W
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  // PX parameters
  arma::vec eta(random_var_ix.n_cols); //random_var_ix is num random effects by num variance components
  eta.fill(1.0);

  for(size_t k=0; k<nr; ++k) {
    r(k) = y[k] - allfit_con[k] - allfit_mod[k];
  }

  Wtr = random_des.t()*r;
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
  arma::vec allfit_random = random_des*gamma;
  if(!randeff) allfit_random.fill(0);

  //allfit_random.fill(0);

  //--------------------------------------------------
  //storage for the fits
  double* allfit = new double[n]; //yhat
  for(size_t i=0;i<n;i++) {
    allfit[i] = allfit_mod[i] + allfit_con[i];
    if(randeff) allfit[i] += allfit_random[i];
  }
  double* ftemp  = new double[n]; //fit of current tree

  // init node pointers
  std::vector<std::vector<tree::tree_cp> > node_pointers_mod(ntree_mod);
  std::vector<std::vector<tree::tree_cp> > node_pointers_con(ntree_con);
  for(size_t j=0; j<ntree_mod; ++j) {
    node_pointers_mod[j].resize(n);
    fit_basis(t_mod[j],xi_mod,di_mod,ftemp,node_pointers_mod[j],true,false);
  }
  for(size_t j=0; j<ntree_con; ++j) {
    node_pointers_con[j].resize(n);
    fit_basis(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],true,vanilla);
  }

  NumericVector sigma_post(nd);
//  NumericVector msd_post(nd);
//  NumericVector bsd_post(nd);
  NumericMatrix m_post(nd,n);
  NumericMatrix yhat_post(nd,n);
  NumericMatrix b_post(nd,n);
  NumericVector eta_con_post(nd);
  NumericVector eta_mod_post(nd);

  NumericMatrix var_prob_con(nd,pi_con.var_probs.size());
  NumericMatrix var_prob_mod(nd,pi_mod.var_probs.size());

  NumericMatrix m_est_post(nd,n_con_est);
  NumericMatrix b_est_post(nd,n_mod_est);

  arma::mat gamma_post(nd,gamma.n_elem);
  arma::mat random_sd_post(nd,random_var.n_elem);

  arma::cube scoefs_mod(di_mod.basis_dim, di_mod.n, nd);
  arma::mat coefs_mod(di_mod.basis_dim, di_mod.n);
  
  arma::cube scoefs_con(di_con.basis_dim, di_con.n, nd);
  arma::mat coefs_con(di_con.basis_dim, di_con.n);

  arma::cube scoefs_con_est(di_con_est.basis_dim, di_con_est.n, nd);
  arma::mat coefs_con_est(di_con_est.basis_dim, di_con_est.n);

  arma::cube scoefs_mod_est(di_mod_est.basis_dim, di_mod_est.n, nd);
  arma::mat coefs_mod_est(di_mod_est.basis_dim, di_mod_est.n);

  //  NumericMatrix spred2(nd,dip.n);

  /*
  //save stuff to tree file
  treef << xi << endl; //cutpoints
  treef << m << endl;  //number of trees
  treef << p << endl;  //dimension of x's
  treef << (int)(nd/thin) << endl;
  */

  //*****************************************************************************
  /* MCMC
   * note: the allfit objects are all carrying the appropriate scales
   */
  //*****************************************************************************
  
  Rcout << "\nBeginning MCMC:\n";
  time_t tp;
  int time1 = time(&tp);

  size_t save_ctr = 0;
  for(size_t i=0;i<(nd*thin+burn);i++) {

    //Rcout << allfit_con[0] << endl;

    if(prior_sample) {
      for(int k=0; k<n; k++) y[k] = gen.normal(allfit[k], sigma);
    }

    if(lower_bd.size()>1) {
      for(int k=0; k<n; k++) {
        //Rcout << y[k] << " " << allfit[k] << " " << sigma << " "<< lower_bd[k] << endl;
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
    if(i%status_interval==0) {
      Rcout << "iteration: " << i << " sigma/SD(Y): "<< sigma << endl;
    }

    //draw trees for m(x)

    //Rcout << "b" << endl;
    for(size_t j=0;j<ntree_con;j++) {

      //Rcout << "tree " << j;
      //Rcout << " first fit" << endl;
      fit_basis(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],false,vanilla);

      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcout << "control tree " << j <<" obs "<< k<<" "<< endl;
          Rcout << t_con[j] << endl;
          stop("nan in ftemp");
        }
        allfit[k] = allfit[k]-pi_con.eta*ftemp[k];
        allfit_con[k] = allfit_con[k]-pi_con.eta*ftemp[k];
        r_con[k] = (y[k]-allfit[k])/pi_con.eta;
        if(r_con[k] != r_con[k]) {
          Rcout << (y[k]-allfit[k]) << endl;
          Rcout << pi_con.eta << endl;
          Rcout << r_con[k] << endl;
          stop("NaN in resid");
        }
      }

      //Rcout << " bd " << endl;
      double aa = bd_basis(t_con[j],xi_con,di_con,pi_con,gen,node_pointers_con[j]);
      //Rcout << " aa " << aa << endl;
      //Rcout << " drmu" << endl;
      drmu_basis(t_con[j],xi_con,di_con,pi_con,gen);
      //Rcout << " second fit" << endl;
      fit_basis(t_con[j],xi_con,di_con,ftemp,node_pointers_con[j],false,vanilla);
      for(size_t k=0;k<n;k++) {
        allfit[k] += pi_con.eta*ftemp[k];
        allfit_con[k] += pi_con.eta*ftemp[k];
      }
    }

    //Rcout << " c" << endl;
    // moderator trees
    for(size_t j=0;j<ntree_mod;j++) {
      fit_basis(t_mod[j],xi_mod,di_mod,ftemp,node_pointers_mod[j], false, false);
      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcout << "moderator tree " << j <<" obs "<< k<<" "<< endl;
          Rcout << t_mod[j] << endl;
          stop("nan in ftemp");
        }

        allfit[k] = allfit[k]-pi_mod.eta*ftemp[k];
        allfit_mod[k] = allfit_mod[k]-pi_mod.eta*ftemp[k];
        r_mod[k] = (y[k]-allfit[k])/pi_mod.eta;
      }

      double aa = bd_basis(t_mod[j],xi_mod,di_mod,pi_mod,gen,node_pointers_mod[j]);
      drmu_basis(t_mod[j],xi_mod,di_mod,pi_mod,gen);
      fit_basis(t_mod[j],xi_mod,di_mod,ftemp,node_pointers_mod[j],false,false);

      for(size_t k=0;k<n;k++) {
        allfit[k] += pi_mod.eta*ftemp[k];
        allfit_mod[k] += pi_mod.eta*ftemp[k];
      }
    }




    //DART
    if(dart & (i>(0.25*burn))) {
      //Rcout << "DART updates" << endl;
      update_dart(t_mod, pi_mod, di_mod, xi_mod, gen);
      update_dart(t_con, pi_con, di_con, xi_con, gen);
      //Rcout << "DART updates complete" << endl;
    }


    //update the PX parameters for control function
    double eta_old;

    for(size_t k=0;k<n;k++) {
      ftemp[k] = y[k] - allfit_mod[k];
      if(randeff) ftemp[k] -= allfit_random[k];
    }
    eta_old = pi_con.eta; // Save previous eta before drawing new one, for adjusting scaling.
    update_scale(ftemp, allfit_con, n, sigma, pi_con, gen);

    for(size_t k=0; k<n; ++k) {
      allfit[k] -= allfit_con[k];
      allfit_con[k] = allfit_con[k] * pi_con.eta / eta_old;
      allfit[k] += allfit_con[k];
    }

    pi_con.sigma = sigma/fabs(pi_con.eta);
    
    //update the PX parameters for moderation function
    
    for(size_t k=0;k<n;k++) {
      ftemp[k] = y[k] - allfit_con[k];
      if(randeff) ftemp[k] -= allfit_random[k];
    }
    eta_old = pi_mod.eta; // Save previous eta before drawing new one, for adjusting scaling.
    
    update_scale(ftemp, allfit_mod, n, sigma, pi_mod, gen);
    
    // Update fits to have new pi.eta scaling.
    for(size_t k=0; k<n; ++k) {
      allfit[k] -= allfit_mod[k];
      allfit_mod[k] = allfit_mod[k] * pi_mod.eta / eta_old;
      allfit[k] += allfit_mod[k];
    }
    //
    //Rcout << pi_mod.eta << endl;
    
    pi_mod.sigma = sigma/fabs(pi_mod.eta);

    //Rcout << "e" << endl;

    if(randeff) {
      //update random effects
      for(size_t k=0; k<n; ++k) {
        r(k) = y[k] - allfit_con[k] - allfit_mod[k];
        allfit[k] -= allfit_random[k];
      }

      Wtr = random_des.t()*r;

      arma::mat adj = diagmat(random_var_ix*eta);
      //    Rcout << adj << endl << endl;
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
      Phi = 0.5*(Phi + Phi.t());
      arma::vec m = adj*Wtr/(sigma*sigma);
      //Rcout << m << Phi << endl << Sigma_inv_random;
      gamma = rmvnorm_post(m, Phi);

      //Rcout << "updated gamma";

      // Update px parameters eta

      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
      Phi2 = 0.5*(Phi2 + Phi2.t());
      eta = rmvnorm_post(m2, Phi2);

      //Rcout << "updated eta";

      // Update variance parameters

      arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
      //Rcout << "A";
      arma::rowvec counts = sum(random_var_ix, 0);
      //Rcout << "B";
      for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df/randeff_scales(ii)*randeff_scales(ii) + ssqs(ii));
      }
      //Rcout << "updated vars" << endl;
      Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

      //Rcout << random_var_ix*random_var;

      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;

      //Rcout << "recom allfit vars" << endl;

      for(size_t k=0; k<n; ++k) {
        allfit[k] = allfit_con[k] + allfit_mod[k] + allfit_random(k); //+= allfit_random[k];
      }
    }

    //draw sigma
    double rss = 0.0;
    double restemp = 0.0;
    for(size_t k=0;k<n;k++) {
      restemp = y[k]-allfit[k];
      rss += restemp*restemp;
    }
    if(!probit) sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    pi_con.sigma = sigma/fabs(pi_con.eta);
    pi_mod.sigma = sigma/fabs(pi_mod.eta);

    if( ((i>=burn) & (i % thin==0)) )  {
      //for(size_t j=0;j<m;j++) treef << t[j] << endl;

//      msd_post(save_ctr) = fabs(pi_con.eta)*con_sd;
//      bsd_post(save_ctr) = fabs(pi_mod.eta)*mod_sd;

      //pi_mod.var_probs

      for(size_t j=0; j<pi_con.var_probs.size(); ++j) {
        var_prob_con(save_ctr, j) = pi_con.var_probs[j];
      }
      for(size_t j=0; j<pi_mod.var_probs.size(); ++j) {
        var_prob_mod(save_ctr, j) = pi_mod.var_probs[j];
      }

      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_sd_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();

      sigma_post(save_ctr) = sigma;
      eta_con_post(save_ctr) = pi_con.eta;
      eta_mod_post(save_ctr) = pi_mod.eta;

      for(size_t k=0;k<n;k++) {
        m_post(save_ctr, k) = allfit_con[k];
        b_post(save_ctr, k) = allfit_mod[k];
        yhat_post(save_ctr, k) = allfit[k];
      }
      
      if(est_con_fits) {
        for(size_t k=0;k<di_con_est.n;k++) {
          m_est_post(save_ctr, k) = pi_con.eta*fit_i_basis(k, t_con, xi_con, di_con_est, vanilla);
        }
      }
      
      if(est_mod_fits) {
        for(size_t k=0;k<di_mod_est.n;k++) {
          b_est_post(save_ctr, k) = pi_mod.eta*fit_i_basis(k, t_mod, xi_mod, di_mod_est, false);
        }
      }

      coefs_mod.zeros();
      for(size_t j=0; j<ntree_mod; ++j) {
        coefs_mod += pi_mod.eta*coef_basis(t_mod[j], xi_mod, di_mod);
      }
      scoefs_mod.slice(save_ctr) = coefs_mod;
      
      coefs_con.zeros();
      for(size_t j=0; j<ntree_con; ++j) {
        coefs_con += pi_con.eta*coef_basis(t_con[j], xi_con, di_con);
      }
      scoefs_con.slice(save_ctr) = coefs_con;

      coefs_mod_est.zeros();
      //if(di_mod_est.n) {
      for(size_t j=0; j<ntree_mod; ++j) {
        coefs_mod_est += pi_mod.eta*coef_basis(t_mod[j], xi_mod, di_mod_est);
      }
      scoefs_mod_est.slice(save_ctr) = coefs_mod_est;

      coefs_con_est.zeros();
      //if(di_con_est.n) {
      for(size_t j=0; j<ntree_con; ++j) {
        coefs_con_est += pi_con.eta*coef_basis(t_con[j], xi_con, di_con_est);
      }
      scoefs_con_est.slice(save_ctr) = coefs_con_est;

      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;

  t_mod.clear(); t_con.clear();
  delete[] allfit;
  delete[] allfit_mod;
  delete[] allfit_con;
  delete[] r_mod;
  delete[] r_con;
  delete[] ftemp;

  treef.close();

  return(List::create(_["yhat_post"] = yhat_post,
                      _["coefs_mod"] = scoefs_mod,
                      _["coefs_con"] = scoefs_con,
                      _["coefs_mod_est"] = scoefs_mod_est,
                      _["coefs_con_est"] = scoefs_con_est,
                      _["var_prob_con"] = var_prob_con,
                      _["var_prob_mod"] = var_prob_mod,
                      _["m_post"] = m_post,
                      _["m_est_post"] = m_est_post,
                      _["b_post"] = b_post, _["b_est_post"] = b_est_post,
                      _["eta_con_post"] = eta_con_post,
                      _["eta_mod_post"] = eta_mod_post,
                      _["sigma"] = sigma_post, //_["msd"] = msd_post, _["bsd"] = bsd_post,
                      _["gamma"] = gamma_post, _["random_sd_post"] = random_sd_post,
                      _["y_last"] = y
  ));
}
