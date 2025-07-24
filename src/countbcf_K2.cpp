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
List countbcf(arma::vec y_,
              arma::vec offset_,
              List bart_specs,
              List bart_designs,
              arma::mat random_des,
              arma::mat random_var, arma::mat random_var_ix,
              double random_var_df, arma::vec randeff_scales,
              int burn, int nd, int thin,
              int count_model,                 // 1=poisson, 2=nb, 3=zip, 4=zinb
              double lambda, double nu,
              double kappa_a, double kappa_b,
              double leaf_c, double leaf_d,
              double z_c, double z_d,
              double kappa_prop_sd = 0.2,
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
  /* ---------- 0.  House-keeping ---------- */
  bool randeff = (random_var_ix.n_elem > 1);
  if (randeff) Rcout << "Using random effects." << std::endl;

  RNGScope scope;
  RNG gen;

  /* ---------- 1.  Read y and offset ---------- */
  std::vector<double> y, offset;
  double miny =  INFINITY, maxy = -INFINITY;
  int n_y0 = -1;                     // # zeros for ZIP / ZINB
  for (auto it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
    if (*it < miny) miny = *it;
    if (*it > maxy) maxy = *it;
    if (*it == 0) ++n_y0;
  }
  size_t n = y.size();
  for (auto it = offset_.begin(); it != offset_.end(); ++it) offset.push_back(*it);

  double sigma = (probit) ? 1.0 : sqrt(((arma::accu(y_%y_) - n*std::pow(arma::mean(y_),2)))/(n-1));

  /* ---------- 2.  Design matrices ---------- */
  size_t num_designs = bart_designs.size();
  std::vector<std::vector<double> > x(num_designs);
  std::vector<std::vector<int> >   groups(num_designs);
  std::vector<bool>                group(num_designs);
  std::vector<xinfo>               x_info(num_designs);
  std::vector<arma::mat>           Omega(num_designs);
  std::vector<size_t>              covariate_dim(num_designs);

  for (size_t i = 0; i < num_designs; ++i) {
    List d = bart_designs[i];
    group[i] = d["group"];
    IntegerVector gt = d["groups"];
    for (auto gi : gt) groups[i].push_back(gi);
    NumericVector xt = d["X"];
    for (auto xi : xt) x[i].push_back(xi);
    covariate_dim[i] = x[i].size() / n;

    List x_info_list = d["info"];
    xinfo xi; xi.resize(covariate_dim[i]);
    for (size_t j = 0; j < covariate_dim[i]; ++j) {
      NumericVector tmp = x_info_list[j];
      std::vector<double> tmp2(tmp.begin(), tmp.end());
      xi[j] = tmp2;
    }
    x_info[i] = xi;
    Omega[i] = as<arma::mat>(d["Omega"]);
  }

  /* ---------- 3.  Forests: 6 = 3 sub-models × (mu, tau) ---------- */
  size_t num_forests = 6;      // f_mu, f_tau, f0_mu, f0_tau, f1_mu, f1_tau
  std::vector<std::vector<tree> >   trees(num_forests);
  std::vector<pinfo>                prior_info(num_forests);
  std::vector<std::vector<double> > allfits(num_forests);
  std::vector<double>               r_tree(n);
  std::vector<dinfo>                di(num_forests);
  std::vector<std::vector<std::vector<tree::tree_cp> > > node_pointers(num_forests);
  std::vector<double>               sample_eta(num_forests);

  double *ftemp = new double[n];
  for (size_t s = 0; s < num_forests; ++s) {
    List spec = bart_specs[s];
    sample_eta[s] = spec["sample_eta"];
    size_t desi = spec["design_index"];
    size_t ntree = spec["ntree"];
    trees[s].resize(ntree);
    prior_info[s].vanilla = spec["vanilla"];

    for (size_t j = 0; j < ntree; ++j) trees[s][j].setm(zeros(Omega[desi].n_rows));

    prior_info[s].pbd = 1.0; prior_info[s].pb = .5;
    prior_info[s].alpha = spec["alpha"];
    prior_info[s].beta  = spec["beta"];
    prior_info[s].sigma = sigma;
    prior_info[s].mu0   = as<arma::vec>(spec["mu0"]);
    prior_info[s].Sigma0 = as<arma::mat>(spec["Sigma0"]);
    prior_info[s].Prec0  = prior_info[s].Sigma0.i();
    prior_info[s].logdetSigma0 = log(det(prior_info[s].Sigma0));
    prior_info[s].eta = 1; prior_info[s].gamma = 1;
    prior_info[s].scale_df = spec["scale_df"];

    // hyper-priors for leaves
    if ((count_model == 3 || count_model == 4) && ((s == 2) || (s == 3))) {
      prior_info[s].c = z_c; prior_info[s].d = z_d;
    } else {
      prior_info[s].c = leaf_c; prior_info[s].d = leaf_d;
    }

    prior_info[s].dart = spec["dart"];
    if (prior_info[s].dart) prior_info[s].dart_alpha = 1.0;
    std::vector<double> vp(covariate_dim[desi], 1.0/covariate_dim[desi]);
    prior_info[s].var_probs = vp;

    dinfo dtemp;
    dtemp.n = n; dtemp.p = covariate_dim[desi];
    dtemp.x = &(x[desi])[0];
    dtemp.y = &r_tree[0];
    dtemp.u_i = &y[0];
    dtemp.offset = &offset[0];
    dtemp.basis_dim = Omega[desi].n_rows;
    dtemp.omega     = &(Omega[desi])[0];
    dtemp.groups    = &(groups[desi])[0];
    dtemp.group     = group[desi];
    di[s] = dtemp;

    allfits[s].resize(n); std::fill(allfits[s].begin(), allfits[s].end(), 0.0);
    node_pointers[s].resize(ntree);
    for (size_t j = 0; j < ntree; ++j) {
      node_pointers[s][j].resize(n);
      fit_loglinear(trees[s][j], x_info[desi], di[s], ftemp, node_pointers[s][j], true, prior_info[s].vanilla);
      for (size_t k = 0; k < n; ++k) allfits[s][k] += ftemp[k];
    }
  }

  /* ---------- 4.  Random effects ---------- */
  size_t random_dim = random_des.n_cols;
  int nr = (randeff) ? n : 1;
  arma::vec r(nr), Wtr(random_dim);
  arma::mat WtW = random_des.t()*random_des;
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*arma::vec(random_var.n_elem, arma::fill::ones)));
  arma::vec eta(random_var_ix.n_cols, arma::fill::ones);
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random,
                          random_des.t()*arma::vec(n, arma::fill::zeros)/(sigma*sigma));
  arma::vec allfit_random = (randeff) ? random_des*gamma : arma::vec(n, arma::fill::zeros);

  /* ---------- 5.  Latent variables ---------- */
  std::vector<int>    z(n, 1);
  std::vector<double> log_w(n, log(0.5));
  std::vector<double> log_w_denom(n, log(2));
  std::vector<double> log_phi(n, 0);
  std::vector<double> log_xi(n, 0);
  std::vector<double> u_vec = y;   // sufficient stat for trees

  double kappa = 1.0, kappa_acpt_rate = 0.0;

  /* ---------- 6.  Output containers ---------- */
  NumericVector sigma_post(nd), kappa_post(nd);
  NumericMatrix etas_post(nd, num_forests);
  arma::mat gamma_post(nd, gamma.n_elem);
  arma::mat random_sd_post(nd, random_var.n_elem);
  std::vector<arma::cube> post_coefs(num_forests);
  for (size_t s = 0; s < num_forests; ++s) {
    post_coefs[s] = arma::cube(di[s].basis_dim, n, nd);
    post_coefs[s].zeros();
  }

  std::vector<tree_samples> final_tree_trace(num_forests);
  for (size_t s = 0; s < num_forests; ++s)
    final_tree_trace[s] = tree_samples(trees[s].size(), di[s].p, nd, di[s].basis_dim, x_info[0]);

  /* ---------- 7.  MCMC ---------- */
  Rcout << "\nBeginning MCMC:\n";
  time_t tp; int time1 = time(&tp);
  size_t save_ctr = 0;

  // Build Z vector (treatment indicator) from bart_designs
  // Assuming last design (index 1) contains 1/0 treatment indicator
  std::vector<double> Z_vec(n);
  for (size_t k = 0; k < n; ++k) Z_vec[k] = x[1][k];

  for (size_t iter = 0; iter < burn + nd*thin; ++iter) {

    if (prior_sample) for (size_t k = 0; k < n; ++k) y[k] = gen.normal(allfit_random[k], sigma);

    /* ---- kappa ---- */
    if (count_model == 2 || count_model == 4) {
      double kappa_star = exp(gen.normal(log(kappa), kappa_prop_sd));
      double log_a_num = ll_loglinear(y, log_mean, kappa_star, count_model, true, log_w, n_y0)
        + R::dbeta(kappa_star/(1+kappa_star), kappa_a, kappa_b, true) - 2*log1p(kappa_star) + log(kappa_star);
      double log_a_denom = ll_loglinear(y, log_mean, kappa, count_model, true, log_w, n_y0)
        + R::dbeta(kappa/(1+kappa), kappa_a, kappa_b, true) - 2*log1p(kappa) + log(kappa);
      if (log(gen.uniform()) < log_a_num - log_a_denom) { kappa = kappa_star; kappa_acpt_rate = ((kappa_acpt_rate*iter)+1)/(iter+1); }
      else { kappa_acpt_rate = ((kappa_acpt_rate*iter)+0)/(iter+1); }
    }

    /* ---- Z_i ---- */
    if (count_model == 3 || count_model == 4) drz_loglinear(z, log_w, kappa, log_mean, n_y0, count_model, gen);

    /* ---- xi (NB) ---- */
    if (count_model == 2 || count_model == 4) drxi_loglinear(log_xi, kappa, log_mean, y, z, gen);

    /* ---- phi (ZI) ---- */
    if (count_model == 3 || count_model == 4) drphi_loglinear(log_phi, log_w_denom, gen);

    /* ---- Forest updates (mu, tau pairs) ---- */
    for (size_t pair = 0; pair < 3; ++pair) {
      size_t mu_idx = 2*pair, tau_idx = 2*pair+1;

      /* --- mu forest (full sample) --- */
      for (size_t j = 0; j < trees[mu_idx].size(); ++j) {
        fit_loglinear(trees[mu_idx][j], x_info[0], di[mu_idx], ftemp, node_pointers[mu_idx][j], false, prior_info[mu_idx].vanilla);
        for (size_t k = 0; k < n; ++k) {
          allfits[mu_idx][k] -= prior_info[mu_idx].eta * ftemp[k];
        }
        for (size_t k = 0; k < n; ++k) {
          if (pair == 0) {           // count model
            u_vec[k] = z[k]*y[k];
            r_tree[k] = z[k]*exp(log_xi[k] + offset[k] + allfits[mu_idx][k] + Z_vec[k]*allfits[tau_idx][k]);
          } else if (pair == 1) {    // f0
            int zero_type = z[k] - 0;
            u_vec[k] = 1 - std::abs(zero_type);
            r_tree[k] = exp(log_phi[k] + allfits[mu_idx][k] + Z_vec[k]*allfits[tau_idx][k]);
          } else {                   // f1
            int zero_type = z[k] - 1;
            u_vec[k] = 1 - std::abs(zero_type);
            r_tree[k] = exp(log_phi[k] + allfits[mu_idx][k] + Z_vec[k]*allfits[tau_idx][k]);
          }
        }
        bd_loglinear(trees[mu_idx][j], x_info[0], di[mu_idx], prior_info[mu_idx], gen, node_pointers[mu_idx][j]);
        drmu_loglinear(trees[mu_idx][j], x_info[0], di[mu_idx], prior_info[mu_idx], gen);
        fit_loglinear(trees[mu_idx][j], x_info[0], di[mu_idx], ftemp, node_pointers[mu_idx][j], false, prior_info[mu_idx].vanilla);
        for (size_t k = 0; k < n; ++k) allfits[mu_idx][k] += prior_info[mu_idx].eta * ftemp[k];
      }

      /* --- tau forest (treated units only) --- */
      for (size_t j = 0; j < trees[tau_idx].size(); ++j) {
        fit_loglinear(trees[tau_idx][j], x_info[0], di[tau_idx], ftemp, node_pointers[tau_idx][j], false, prior_info[tau_idx].vanilla);
        for (size_t k = 0; k < n; ++k) {
          allfits[tau_idx][k] -= prior_info[tau_idx].eta * ftemp[k];
        }
        for (size_t k = 0; k < n; ++k) {
          r_tree[k] = Z_vec[k]*(y[k] - (pair==0 ? exp(log_xi[k] + offset[k] + allfits[mu_idx][k]) : 0)); // residual for tau
        }
        bd_loglinear(trees[tau_idx][j], x_info[0], di[tau_idx], prior_info[tau_idx], gen, node_pointers[tau_idx][j]);
        drmu_loglinear(trees[tau_idx][j], x_info[0], di[tau_idx], prior_info[tau_idx], gen);
        fit_loglinear(trees[tau_idx][j], x_info[0], di[tau_idx], ftemp, node_pointers[tau_idx][j], false, prior_info[tau_idx].vanilla);
        for (size_t k = 0; k < n; ++k) allfits[tau_idx][k] += prior_info[tau_idx].eta * ftemp[k];
      }
    }

    /* ---- Update w(x_i) for ZI ---- */
    if (count_model == 3 || count_model == 4) {
      for (size_t k = 0; k < n; ++k) {
        double f0_val = allfits[2][k] + Z_vec[k]*allfits[3][k];
        double f1_val = allfits[4][k] + Z_vec[k]*allfits[5][k];
        log_w_denom[k] = logsumexp(f0_val, f1_val);
        log_w[k] = f1_val - log_w_denom[k];
      }
    }

    /* ---- Random effects (unchanged) ---- */
    if (randeff) {
      for (size_t k = 0; k < n; ++k) {
        r(k) = y[k];
        for (size_t s = 0; s < num_forests; ++s) r(k) -= allfits[s][k];
        r(k) += allfit_random[k];
      }
      arma::mat adj = diagmat(random_var_ix*eta);
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + arma::diagmat(1/(random_var_ix*arma::vec(random_var.n_elem, arma::fill::ones)));
      Phi = 0.5*(Phi + Phi.t());
      arma::vec m = adj*(random_des.t()*r)/(sigma*sigma);
      gamma = rmvnorm_post(m, Phi);

      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*(random_des.t()*r)/(sigma*sigma);
      Phi2 = 0.5*(Phi2 + Phi2.t());
      eta = rmvnorm_post(m2, Phi2);

      arma::vec ssqs = random_var_ix.t()*(gamma % gamma);
      arma::rowvec counts = sum(random_var_ix, 0);
      for (size_t ii = 0; ii < random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*
          2.0/(random_var_df/randeff_scales(ii)*randeff_scales(ii) + ssqs(ii));
      }
      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;
    }

    /* ---- Update log_mean ---- */
    for (size_t k = 0; k < n; ++k) {
      // count model (f_mu + Z*f_tau)
      log_mean[k] = offset[k] + allfits[0][k] + Z_vec[k]*allfits[1][k];
    }

    /* ---- Save ---- */
    if (iter >= burn && iter % thin == 0) {
      sigma_post(save_ctr) = sigma;
      kappa_post(save_ctr) = kappa;
      for (size_t s = 0; s < num_forests; ++s) {
        etas_post(save_ctr, s) = prior_info[s].eta;
        for (size_t j = 0; j < trees[s].size(); ++j) {
          post_coefs[s].slice(save_ctr) += prior_info[s].eta*
            coef_basis(trees[s][j], x_info[0], di[s]);
          final_tree_trace[s].t[save_ctr][j] = trees[s][j];
          final_tree_trace[s].t[save_ctr][j].compress();
          final_tree_trace[s].t[save_ctr][j].scale(prior_info[s].eta);
        }
      }
      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_sd_post.row(save_ctr) = sqrt(eta % eta % random_var).t();
      ++save_ctr;
    }
  }

  Rcout << "time for loop: " << time(nullptr) - time1 << " seconds\n";
  delete[] ftemp;

  /* ---------- 8.  Serialize trees ---------- */
  std::vector<Rcpp::RawVector> serial_streams(num_forests);
  std::vector<Rcpp::CharacterVector> Rtree_streams(num_forests);
  for (size_t s = 0; s < num_forests; ++s) {
    Rtree_streams[s] = final_tree_trace[s].save_string();
    std::stringstream ss; cereal::BinaryOutputArchive oa(ss);
    oa(final_tree_trace[s]);
    ss.seekg(0, ss.end); RawVector rv(ss.tellg());
    ss.seekg(0, ss.beg); ss.read(reinterpret_cast<char*>(&rv[0]), rv.size());
    serial_streams[s] = rv;
  }

  return List::create(
    _["yhat_post"]        = NumericMatrix(),   // could expose actual fits if desired
    _["coefs"]            = post_coefs,
    _["etas"]             = etas_post,
    _["sigma"]            = sigma_post,
    _["kappa"]            = kappa_post,
    _["kappa_acceptance"] = kappa_acpt_rate,
    _["gamma"]            = gamma_post,
    _["random_sd_post"]   = random_sd_post,
    _["tree_streams"]     = Rtree_streams,
    _["tree_serials"]     = serial_streams,
    _["tree_trace"]       = final_tree_trace,
    _["y_last"]           = y
  );
}