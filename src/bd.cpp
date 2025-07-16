#include <iostream>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funs.h"

using std::cout;
using std::endl;

double bd_basis(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers)
{
  tree::npv goodbots;                    // nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots);  // prob of a birth at x
  
  //Rcpp::Rcout << "PBx: " << PBx << endl;
  
  // If statement for selecting birth or death proposal.
  if(gen.uniform() < PBx) {
    
    //Rcpp::Rcout << "birth" << endl;
    
    //--------------------------------------------------
    // BIRTH PROPOSAL
    //--------------------------------------------------
    
    //--------------------------------------------------
    //draw proposal
    
    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(gen.uniform()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at
    
    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars);
    //size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    //size_t v = goodvars[vi];
    
    // DART
    size_t vi;
    if(pi.dart) {
      std::vector<double> wt(goodvars.size());
      for(size_t j=0; j<goodvars.size(); ++j) {
        //Rcpp::Rcout << goodvars[j] << " ";
        wt[j] = log(pi.var_probs[goodvars[j]]);
        //Rcpp::Rcout << wt[j] << endl;
      }
      // vi = goodvars[rdisc_log_inplace(wt)];
      vi = rdisc_log_inplace(wt);
    } else {
      vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    }
    size_t v = goodvars[vi];
    
    //Rcpp::Rcout << "goodvars size " <<goodvars.size() << endl;
    //Rcpp::Rcout << "vi " <<vi << endl;
    //Rcpp::Rcout << "v " <<v << endl;
    
    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points
    
    //--------------------------------------------------
    //compute things needed for metropolis ratio
    
    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx
    
    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }
    
    double PDy; //prob of proposing death at y
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }
    
    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
      //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
      if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }
    
    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    sl.sy = 0;
    sl.n0 = 0.0;
    sl.sy_vec.zeros(di.basis_dim);
    sl.WtW.zeros(di.basis_dim, di.basis_dim);
    sr.sy = 0;
    sr.n0 = 0.0;
    sr.sy_vec.zeros(di.basis_dim);
    sr.WtW.zeros(di.basis_dim, di.basis_dim);
    
    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, si ze_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_basis_birth(x,nx,v,c,xi,di,bnv,sv,sl,sr,node_pointers);
    
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }
    
    //--------------------------------------------------
    //compute alpha
    
    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;
    
    if((sl.n_unique>=5) && (sr.n_unique>=5)) { //cludge?
      
      //lil_basis(double n, vec sy_vec, mat WtW, double sigma, vec mu0, mat Prec0, mat Sigma0)
      
      // lill = lil_basis(sl.n, sl.sy_vec, sl.WtW, pi);                        // USED lil(...)
      // lilr = lil_basis(sr.n, sr.sy_vec, sr.WtW, pi);                        // USED lil(...)
      // lilt = lil_basis(sl.n + sr.n, sl.sy_vec + sr.sy_vec, sl.WtW + sr.WtW, pi);  // USED lil(...)
      // 
      lill = lil_basis(sl, pi);                        // USD lil(...)
      lilr = lil_basis(sr, pi);                        // USED lil(...)
      lilt = lil_basis(sl,sr,pi);  // USED lil(...)
      
      
      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);
      
    } else {
      alpha=0.0;
    }
    
    /*
    cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
    cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
    cout << "L,U: " << L << "," << U << endl;
    cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
      PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
        ", " << Pnogy << "," << Pbotx << endl;
    cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
    cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
    cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
    cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
    */
    
    //--------------------------------------------------
    //finally ready to try metrop
    //--------------------------------------------------
    
    vec mul, mur;
    mul = zeros(di.basis_dim);
    mur = zeros(di.basis_dim);
    if(gen.uniform() < alpha) {
      
      //--------------------------------------------------
      // do birth:
      // Set mul and mur to zero vectors, since we will immediately
      // fill them in the MCMC by using drmu.  This saves computation cost.
      
      //Rcpp::Rcout << "Prebirth... " << sr.sy << " " << sr.sy_vec << endl;
      //Rcpp::Rcout << "Prebirth... " << sl.sy << " " << sl.sy_vec << endl;
      
      //Rcpp::Rcout << "Birthing";
      
      
      //Rcpp::Rcout << " before birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      x.birth(nx->nid(),v,c,mul,mur,sl,sr);
      
      double *xx;        //current x
      for(size_t i=0;i<di.n;i++) {
        xx = di.x + i*di.p;  
        bool in_candidate_nog, left, right;
        in_candidate_nog = (node_pointers[i] == nx);
        left  = in_candidate_nog & (xx[v] < xi[v][c]);
        right = in_candidate_nog & !(xx[v] < xi[v][c]);
        if(left)  node_pointers[i] = nx->getl();
        if(right) node_pointers[i] = nx->getr();
      }
      
      //Rcpp::Rcout << " after birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      //Rcpp::Rcout << " Birthed";
      return alpha;
    } else {
      /*
      sinfo st = sl;
      st.n += sr.n;
      st.n0 += sr.n0;
      st.sy += sr.sy;
      st.sy_vec += sr.sy_vec;
      st.WtW += sr.WtW;
      
      nx->s = st;
       */
      
      return alpha+10;
    }
  } else {
    //    Rcpp::Rcout << "death" << endl;
    //--------------------------------------------------
    // DEATH PROPOSAL
    //--------------------------------------------------
    
    //--------------------------------------------------
    //draw proposal
    
    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at
    
    //--------------------------------------------------
    //compute things needed for metropolis ratio
    
    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);
    
    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);
    
    double PBy;  //prob of birth move at y
    //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }
    
    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;
    
    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();
    
    //--------------------------------------------------
    //compute sufficient statistics
    //sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    // Test if problem is getl and getr.
    //getsuff_basis(x,nx->getl(),nx->getr(),xi,di,sl,sr);
    
    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_basis(x,xi,di,bnv,sv,node_pointers);
    
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }
    
    sl = nx->getl()->s;
    sr = nx->getr()->s;
    
    double lill = lil_basis(sl, pi);                        // USED lil(...)
    double lilr = lil_basis(sr, pi);                        // USED lil(...)
    double lilt = lil_basis(sl,sr, pi);  // USED lil(...)
    
    
    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);
    
    //--------------------------------------------------
    //finally ready to try metrop
    
    // All values wrong, but ok bc get reset immediately in the draw mu.
    if(gen.uniform()<alpha) {
      
      tree::tree_cp lptr = nx->getl();
      tree::tree_cp rptr = nx->getr();
      
      for(size_t i=0;i<di.n;i++) {
        if((node_pointers[i]==lptr) | (node_pointers[i]==rptr)) {
          node_pointers[i] = nx;
        }
      }
      
      //Rcpp::Rcout << "Deathing";
      
      //draw mu for nog (which will be bot)
      vec mu;
      mu = zeros(di.basis_dim);
      
      sinfo st = sl;
      st.n += sr.n;
      st.n0 += sr.n0;
      st.sy += sr.sy;
      st.sy_vec += sr.sy_vec;
      st.WtW += sr.WtW;
      
      //do death
      x.death(nx->nid(),mu, st);
      
      return -alpha;
    } else {
      
      //nx->getl()->s = sl;
      //nx->getr()->s = sr;
      
      return -alpha-10;
    }
  }
}


double bd_loglinear(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers)
{
  tree::npv goodbots;                    // nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots);  // prob of a birth at x
  
  //Rcpp::Rcout << "PBx: " << PBx << endl;
  
  // If statement for selecting birth or death proposal.
  if(gen.uniform() < PBx) {
    
    //Rcpp::Rcout << "birth" << endl;
    
    //--------------------------------------------------
    // BIRTH PROPOSAL
    //--------------------------------------------------
    
    //--------------------------------------------------
    //draw proposal
    
    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(gen.uniform()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at
    
    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars);
    //size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    //size_t v = goodvars[vi];
    
    // DART
    size_t vi;
    if(pi.dart) {
      std::vector<double> wt(goodvars.size());
      for(size_t j=0; j<goodvars.size(); ++j) {
        //Rcpp::Rcout << goodvars[j] << " ";
        wt[j] = log(pi.var_probs[goodvars[j]]);
        //Rcpp::Rcout << wt[j] << endl;
      }
      // vi = goodvars[rdisc_log_inplace(wt)];
      vi = rdisc_log_inplace(wt);
    } else {
      vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
    }
    size_t v = goodvars[vi];
    
    //Rcpp::Rcout << "goodvars size " <<goodvars.size() << endl;
    //Rcpp::Rcout << "vi " <<vi << endl;
    //Rcpp::Rcout << "v " <<v << endl;
    
    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points
    
    //--------------------------------------------------
    //compute things needed for metropolis ratio
    
    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx
    
    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }
    
    double PDy; //prob of proposing death at y
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }
    
    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
      //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
      if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }
    
    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    sl.sy = 0;
    sl.n0 = 0.0;
    sl.sy_vec.zeros(di.basis_dim);
    sl.WtW.zeros(di.basis_dim, di.basis_dim);
    sr.sy = 0;
    sr.n0 = 0.0;
    sr.sy_vec.zeros(di.basis_dim);
    sr.WtW.zeros(di.basis_dim, di.basis_dim);
    
    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, si ze_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_loglinear_birth(x,nx,v,c,xi,di,bnv,sv,sl,sr,node_pointers);
    
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }
    
    //--------------------------------------------------
    //compute alpha
    
    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;
    
    if((sl.n_unique>=5) && (sr.n_unique>=5)) { //cludge?
      
      //lil_basis(double n, vec sy_vec, mat WtW, double sigma, vec mu0, mat Prec0, mat Sigma0)
      
      // lill = lil_basis(sl.n, sl.sy_vec, sl.WtW, pi);                        // USED lil(...)
      // lilr = lil_basis(sr.n, sr.sy_vec, sr.WtW, pi);                        // USED lil(...)
      // lilt = lil_basis(sl.n + sr.n, sl.sy_vec + sr.sy_vec, sl.WtW + sr.WtW, pi);  // USED lil(...)
      // 
      lill = lil_loglinear(sl, pi);                        // USD lil(...)
      lilr = lil_loglinear(sr, pi);                        // USED lil(...)
      lilt = lil_loglinear(sl,sr,pi);  // USED lil(...)
      
      
      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);
      
    } else {
      alpha=0.0;
    }
    
    //
    //cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
    //cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
    //cout << "L,U: " << L << "," << U << endl;
    //cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
    //  PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
    //    ", " << Pnogy << "," << Pbotx << endl;
    //cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
    //cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
    //cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
    //cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
    //
    
    //--------------------------------------------------
    //finally ready to try metrop
    //--------------------------------------------------
    
    vec mul, mur;
    mul = zeros(di.basis_dim);
    mur = zeros(di.basis_dim);
    if(gen.uniform() < alpha) {
      
      //--------------------------------------------------
      // do birth:
      // Set mul and mur to zero vectors, since we will immediately
      // fill them in the MCMC by using drmu.  This saves computation cost.
      
      //Rcpp::Rcout << "Prebirth... " << sr.sy << " " << sr.sy_vec << endl;
      //Rcpp::Rcout << "Prebirth... " << sl.sy << " " << sl.sy_vec << endl;
      
      //Rcpp::Rcout << "Birthing";
      
      
      //Rcpp::Rcout << " before birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      x.birth(nx->nid(),v,c,mul,mur,sl,sr);
      
      double *xx;        //current x
      for(size_t i=0;i<di.n;i++) {
        xx = di.x + i*di.p;  
        bool in_candidate_nog, left, right;
        in_candidate_nog = (node_pointers[i] == nx);
        left  = in_candidate_nog & (xx[v] < xi[v][c]);
        right = in_candidate_nog & !(xx[v] < xi[v][c]);
        if(left)  node_pointers[i] = nx->getl();
        if(right) node_pointers[i] = nx->getr();
      }
      
      //Rcpp::Rcout << " after birth" << endl;
      //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
      //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
      //Rcpp::Rcout << " Birthed";
      return alpha;
    } else {
      

      //sinfo st = sl;
      //st.n += sr.n;
      //st.n0 += sr.n0;
      //st.sy += sr.sy;
      //st.sy_vec += sr.sy_vec;
      //st.WtW += sr.WtW;
      //nx->s = st;
      
      return alpha+10;
    }
  } else {
    //    Rcpp::Rcout << "death" << endl;
    //--------------------------------------------------
    // DEATH PROPOSAL
    //--------------------------------------------------
    
    //--------------------------------------------------
    //draw proposal
    
    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at
    
    //--------------------------------------------------
    //compute things needed for metropolis ratio
    
    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);
    
    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);
    
    double PBy;  //prob of birth move at y
    //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }
    
    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;
    
    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();
    
    //--------------------------------------------------
    //compute sufficient statistics
    //sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    // Test if problem is getl and getr.
    //getsuff_basis(x,nx->getl(),nx->getr(),xi,di,sl,sr);
    
    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl(di.basis_dim),sr(di.basis_dim); //sl for left from nx and sr for right from nx (using rule (v,c))
    
    // void allsuff_basis_birth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv, sinfo& sl, sinfo& sr)
    tree::npv bnv;
    std::vector<sinfo> sv; //will be resized in allsuff
    //allsuff_basis(t,xi,di,bnv0,sv0);
    allsuff_loglinear(x,xi,di,bnv,sv,node_pointers);
    
    for(size_t tt=0; tt<bnv.size(); ++tt) {
      bnv[tt]->s = sv[tt];
    }
    
    sl = nx->getl()->s;
    sr = nx->getr()->s;
    
    double lill = lil_loglinear(sl, pi);                        // USED lil(...)
    double lilr = lil_loglinear(sr, pi);                        // USED lil(...)
    double lilt = lil_loglinear(sl,sr, pi);  // USED lil(...)
    
    
    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);
    
    //--------------------------------------------------
    //finally ready to try metrop
    
    // All values wrong, but ok bc get reset immediately in the draw mu.
    if(gen.uniform()<alpha) {
      
      tree::tree_cp lptr = nx->getl();
      tree::tree_cp rptr = nx->getr();
      
      for(size_t i=0;i<di.n;i++) {
        if((node_pointers[i]==lptr) | (node_pointers[i]==rptr)) {
          node_pointers[i] = nx;
        }
      }
      
      //Rcpp::Rcout << "Deathing";
      
      //draw mu for nog (which will be bot)
      vec mu;
      mu = zeros(di.basis_dim);
      
      sinfo st = sl;
      st.n += sr.n;
      st.n0 += sr.n0;
      st.sy += sr.sy;
      st.sy_vec += sr.sy_vec;
      st.WtW += sr.WtW;
      
      //do death
      x.death(nx->nid(),mu, st);
      
      return -alpha;
    } else {
      
      //nx->getl()->s = sl;
      //nx->getr()->s = sr;
      
      return -alpha-10;
    }
  }
}

