#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>
#include <cereal/types/vector.hpp>

#include "info.h"
#include "rng.h"

/*
Class is a bit confusing because sometimes you
want to think about an instance as a node and
sometimes as the whole tree including the node
and all its children.
Fundamentally, we have a tree not a node.
*/

/*
three ways to access a node:
(i) node id: is the integer assigned by the node numbering system
     assuming they are all there
(ii) node ind: is the index into the array of
   node(tree) pointers returned by getnodes or getbots or getnogs
   which means you go left to right across the bottom of the tree
(iii) by its pointer (should only "give out" const pointers)
*/

//info contained in a node, used by input operator
struct node_info {
   std::size_t id; //node id
   std::size_t v;  //variable
   std::size_t c;  //cut point
   vec m;          //mu CHANGED: was double m;
};

struct node_info_serial {
   std::size_t id; //node id
   std::size_t v;  //variable
   std::size_t c;  //cut point
   std::vector<double> m;
   
   template<class Archive>
   void serialize(Archive & archive)
   {
      archive( id, v, c, m ); 
   }
};

typedef std::vector<node_info_serial> flat_tree;

//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
// left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

class tree {
public:
  
  sinfo s;

   //------------------------------
   //typedefs
   typedef tree* tree_p;
   typedef const tree* tree_cp;
   typedef std::vector<tree_p> npv; //Node Pointer Vector
   typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector
   
   //------------------------------
   //friends
   friend std::istream& operator>>(std::istream&, tree&);
#ifdef MPIBART
   friend double bd(tree& x, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves);
#else
   friend double bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
   friend double bd_basis(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers);
   friend double bd_loglinear(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen, std::vector<tree::tree_cp>& node_pointers);
   friend void fit_basis(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla);
   friend void fit_loglinear(tree& t, xinfo& xi, dinfo& di, double* fv, std::vector<tree::tree_cp>& node_pointers, bool populate, bool vanilla);
   friend double fit_i_basis(size_t& i, std::vector<tree>& t, xinfo& xi, dinfo& di, bool vanilla);
   friend mat coef_basis(tree& t, xinfo& xi, dinfo& di); 
   //friend bool bd_rj(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
   friend bool bdprec(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
   friend double bdhet(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen);
   friend void unflatten(Rcpp::List&, tree&);
#endif

   //------------------------------
   //tree constructors, destructors
   tree();
   tree(const tree&);
   tree(vec);        
   tree(flat_tree&);
   
   ~tree() {tonull();}

   //------------------------------
   //operators
   tree& operator=(const tree&);

   //------------------------------
   //node access
   //you are freely allowed to change mu, v, c
   //set----------
   void setm(vec mu) {this->mu=mu;}   // UPDATED: WAS void setm(double mu) {this->mu=mu;}
   void setm(std::vector<double> mu) {this->mu=arma::conv_to< arma::vec >::from(mu);}   // UPDATED: WAS void setm(double mu) {this->mu=mu;}
   void setm(double mu) {this->mu[0]=mu;}   // UPDATED: WAS void setm(double mu) {this->mu=mu;}
   
   void setv(size_t v) {this->v = v;}
   void setc(size_t c) {this->c = c;}
   //get----------
   vec getm() const {return mu;}    // UPDATED: WAS double getm() const {return mu;}
   vec getm(uvec idx) const {return mu(idx);}  // NEW FUNCTION: returns specific indices of mu.

   double getm(double tval, vec tref) const {
      uvec idx = find(tref==tval);
      return as_scalar(mu(idx));} // NEW FUNCTION: Returns single index of mu, for time index matching tval input.

   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() const {return p;}  //should this be tree_cp?
   tree_p getl() const {return l;}
   tree_p getr() const {return r;}

   //------------------------------
   //tree functions
   //stuff about the tree----------
   size_t treesize() const; //number of nodes in tree
   size_t nnogs() const;    //number of nog nodes (no grandchildren nodes)
   size_t nbots() const;    //number of bottom nodes
   void pr(bool pc=true) const; //to screen, pc is "print children"
   //birth death using nid----------
   bool birth(size_t nid,size_t v, size_t c, vec &ml, vec &mr, sinfo &sl, sinfo &sr);  //UPDATED: WAS ..., double ml, double mr)
   bool birth(size_t nid,size_t v, size_t c, vec &ml, vec &mr);  //UPDATED: WAS ..., double ml, double mr)
     
   bool death(size_t nid, vec &mu);
   bool death(size_t nid, vec &mu, sinfo &s);
   //vectors of node pointers----------
   void getbots(npv& bv);         //get bottom nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnobots(npv& v);        //get all nodes except leafs
   void getnodes(cnpv& v) const;  //get all nodes
   //find node from x and region for var----------
   tree_cp bn(double *x,xinfo& xi);  //find bottom node for x
   void rg(size_t v, int* L, int* U) const; //find region [L,U] for var v.
   size_t nuse(size_t v); //how many times var v is used in a rule.
   void varsplits(std::set<size_t> &splits, size_t v); //splitting values for var v.
   //------------------------------
   //node functions
   size_t depth() const; //depth of a node
   size_t nid() const;   //node id
   char ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
   tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there.
   bool isnog() const;
   void tonull(); //like a "clear", null tree has just one node
   
   void compress();
   void scale(double);
   
   //Rcpp::List flatten(double);
   Rcpp::NumericMatrix flatten_mat(double);
   flat_tree flatten() const;
   void make_from_flat(flat_tree&);
   

private:
   //------------------------------
   //parameter for node
   vec mu; // UPDATED: WAS double mu;
   //------------------------------
   //rule: left if x[v] < xinfo[v][c]
   size_t v;
   size_t c;
   //------------------------------
   //tree structure
   tree_p p; //parent
   tree_p l; //left child
   tree_p r; //right child
   //------------------------------
   //utiity functions
   void cp(tree_p n,  tree_cp o); //copy tree o to n
   //void birthp(tree_p np,size_t v, size_t c, vec ml, vec mr);
   //void deathp(tree_p nb, vec mu); //kill children of nog node nb
};

std::istream& operator>>(std::istream&, tree&);
std::ostream& operator<<(std::ostream&, const tree&);
std::istream& operator>>(std::istream&, xinfo&);
std::ostream& operator<<(std::ostream&, const xinfo&);

void unflatten(Rcpp::List&, tree&);


#endif
