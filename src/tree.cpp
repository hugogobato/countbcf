#include <string>
#include <vector>
#include <map>
#include "tree.h"

using std::string;
using std::cout;
using std::endl;

using namespace arma;
using namespace Rcpp;

//--------------------------------------------------
// constructors
tree::tree(): mu(zeros(1)),v(0),c(0),p(0),l(0),r(0) {}     // UPDATED: was mu(0.0)
tree::tree(vec m): mu(m),v(0),c(0),p(0),l(0),r(0) {}     // UPDATED: was tree::tree(double m):
tree::tree(const tree& n): mu(zeros(1)),v(0),c(0),p(0),l(0),r(0) { cp(this,&n); } // UPDATED: was mu(0.0)
tree::tree(flat_tree& n): mu(zeros(1)),v(0),c(0),p(0),l(0),r(0) { 
  this->make_from_flat(n);
} 

void tree::make_from_flat(flat_tree& nv) {
  this->tonull();
  size_t tid,pid; //tid: id of current node, pid: parent's id
  std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
  size_t nn = nv.size(); //number of nodes
  
  //t.tonull(); // obliterate old tree (if there)
  
  //std::vector<node_info> nv(nn);   // Vector of node_info objects.
  
  /*
  for(size_t i=0;i!=nn;i++) {
    //int b_size; //basis size
    //is >> nv[i].id >> nv[i].v >> nv[i].c >> b_size;
    
    List ft = flat_tree[i];
    //Rcout << "ft size " << ft.size();
    
    nv[i].id = ft["nid"];
    nv[i].v = ft["v"];
    nv[i].c = ft["c"];
    nv[i].m = as<arma::vec>(ft["m"]);
    
  }
  */
  
  //first node has to be the top one
  pts[1] = this; //careful! this is not the first pts, it is pointer of id 1.
  this->setv(nv[0].v); this->setc(nv[0].c); this->setm(nv[0].m);
  this->p=0;
  
  //now loop through the rest of the nodes knowing parent is already there.
  for(size_t i=1;i!=nv.size();i++) {
    tree::tree_p np = new tree;
    np->v = nv[i].v; np->c=nv[i].c; np->mu=nv[i].m;
    tid = nv[i].id;
    pts[tid] = np;
    pid = tid/2;
    // set pointers
    if(tid % 2 == 0) { //left child has even id
      pts[pid]->l = np;
    } else {
      pts[pid]->r = np;
    }
    np->p = pts[pid];
  }
  return;
}


//--------------------------------------------------
//operators
tree& tree::operator=(const tree& rhs)
{
   if(&rhs != this) {
      tonull(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}
//--------------------------------------------------
//public functions
// find bottom node pointer given x
//--------------------
tree::tree_cp tree::bn(double *x,xinfo& xi)
{
   if(l==0) return this; //bottom node
   if(x[v] < xi[v][c]) {
      return l->bn(x,xi);
   } else {
      return r->bn(x,xi);
   }
}
//--------------------
//find region for a given variable
void tree::rg(size_t v, int* L, int* U) const
{
   if(p==0)  { //no parent
      return;
   }
   if(p->v == v) { //does my parent use v?
      if(this == p->l) { //am I left or right child
         if((int)(p->c) <= (*U)) *U = (p->c)-1;
      } else {
         if((int)(p->c) >= *L) *L = (p->c)+1;
      }
   }
   p->rg(v,L,U);
}
//--------------------
//tree size
size_t tree::treesize() const
{
   if(l==0) return 1;  //if bottom node, tree size is 1
   else return (1+l->treesize()+r->treesize());
}
//--------------------
size_t tree::nnogs() const
{
   if(!l) return 0; //bottom node
   if(l->l || r->l) { //not a nog
      return (l->nnogs() + r->nnogs());
   } else { //is a nog
      return 1;
   }
}
size_t tree::nuse(size_t v)
{
   npv nds;
   this->getnodes(nds);
   size_t nu=0; //return value
   for(size_t i=0;i!=nds.size();i++) {
      if(nds[i]->l && nds[i]->v==v) nu+=1;
   }
   return nu;
}

void tree::varsplits(std::set<size_t> &splits, size_t v)
{
   npv nds;
   this->getnodes(nds);
   //size_t nu=0; //return value
   //std::set out;
   for(size_t i=0;i!=nds.size();i++) {
      if(nds[i]->l && nds[i]->v==v) {
        splits.insert(nds[i]->c); //c is index of split rule
      }
   }
}

//--------------------
size_t tree::nbots() const
{
   if(l==0) { //if a bottom node
      return 1;
   } else {
      return l->nbots() + r->nbots();
   }
}
//--------------------
//depth of node
size_t tree::depth() const
{
   if(!p) return 0; //no parents
   else return (1+p->depth());
}
//--------------------
// node id
size_t tree::nid() const
//recursion up the tree
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->nid()); //if you are a left child
   else return 2*(p->nid())+1; //else you are a right child
}
//--------------------
//node type
char tree::ntype() const
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}
//--------------------
//get bottom nodes
//recursion down the tree
void tree::getbots(npv& bv)
{
   if(l) { //have children
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}
//--------------------
//get nog nodes
//recursion down the tree
void tree::getnogs(npv& nv)
{
   if(l) { //have children
      if((l->l) || (r->l)) {  //have grandchildren
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}
//--------------------
//get all nodes
//recursion down the tree
void tree::getnodes(npv& v)
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
// Get nodes that are not leafs
void tree::getnobots(npv& v)
{
  if(this->l) //left node has children
  {
    v.push_back(this);
    this->l->getnobots(v);
    if(this->r->l)
      this->r->getnobots(v);
  }
}
void tree::getnodes(cnpv& v)  const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
//--------------------
//add children to  bot node nid
bool tree::birth(size_t nid,size_t v, size_t c, vec &ml, vec &mr)  //UPDATED: WAS ..., double ml, double mr)
{
   tree_p np = getptr(nid);
   if(np==0) {
      Rcout << "error in birth: bottom node not found\n";
      return false; //did not find note with that nid
   }
   if(np->l) {
      Rcout << "error in birth: found node has children\n";
      return false; //node is not a bottom node
   }

   //add children to bottom node np
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;

   return true;
}
//add children to  bot node nid
bool tree::birth(size_t nid,size_t v, size_t c, vec &ml, vec &mr, sinfo &sl, sinfo &sr)  //UPDATED: WAS ..., double ml, double mr)
{
  tree_p np = getptr(nid);
  if(np==0) {
    Rcout << "error in birth: bottom node not found\n";
    return false; //did not find note with that nid
  }
  if(np->l) {
    Rcout << "error in birth: found node has children\n";
    return false; //node is not a bottom node
  }
  
  //Rcpp::Rcout << " during birth" << endl;
  //Rcpp::Rcout << "L" << endl << sl.WtW << endl;
  //Rcpp::Rcout << "R" << endl << sr.WtW << endl;
  
  //add children to bottom node np
  tree_p l = new tree;
  l->mu = ml;
  l->s = sl;
  
  tree_p r = new tree;
  r->mu = mr;
  r->s = sr;
  
  //Rcpp::Rcout << " during birth 2" << endl;
  //Rcpp::Rcout << "L" << endl << l->s.WtW << endl;
  //Rcpp::Rcout << "R" << endl << r->s.WtW << endl;
  
  np->l=l;
  np->r=r;
  np->v = v; np->c=c;
  l->p = np;
  r->p = np;
  
  return true;
}


//--------------------
//is the node a nog node
bool tree::isnog() const
{
   bool isnog=true;
   if(l) {
      if(l->l || r->l) isnog=false; //one of the children has children.
   } else {
      isnog=false; //no children
   }
   return isnog;
}
//--------------------
//kill children of  nog node nid
bool tree::death(size_t nid, vec &mu)   // UPDATED, WAS:  bool tree::death(size_t nid, double mu)
{
   tree_p nb = getptr(nid);
   if(nb==0) {
      Rcout << "error in death, nid invalid\n";
      return false;
   }
   if(nb->isnog()) {
      delete nb->l;
      delete nb->r;
      nb->l=0;
      nb->r=0;
      nb->v=0;
      nb->c=0;
      nb->mu=mu;
      return true;
   } else {
      Rcout << "error in death, node is not a nog node\n";
      return false;
   }
}

bool tree::death(size_t nid, vec &mu, sinfo &s)   // UPDATED, WAS:  bool tree::death(size_t nid, double mu)
{
  tree_p nb = getptr(nid);
  if(nb==0) {
    Rcout << "error in death, nid invalid\n";
    return false;
  }
  if(nb->isnog()) {
    delete nb->l;
    delete nb->r;
    nb->s=s;
    nb->l=0;
    nb->r=0;
    nb->v=0;
    nb->c=0;
    nb->mu=mu;
    return true;
  } else {
    Rcout << "error in death, node is not a nog node\n";
    return false;
  }
}


//--------------------
//add children to bot node *np
/*
void tree::birthp(tree_p np,size_t v, size_t c, vec ml, vec mr)  //UPDATED, WAS: ...,  double ml, double mr)
{
  Rcout << "birthp";
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;
}

//--------------------
//kill children of  nog node *nb
void tree::deathp(tree_p nb, vec mu)   // UPDATED, WAS: void tree::deathp(tree_p nb, double mu)
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->mu=mu;
}
 
 */

//--------------------
//print out tree(pc=true) or node(pc=false) information
//uses recursion down
void tree::pr(bool pc) const
{
   size_t d = depth();
   size_t id = nid();

   size_t pid;
   if(!p) pid=0; //parent of top node
   else pid = p->nid();

   string pad(2*d,' ');
   string sp(", ");
   if(pc && (ntype()=='t'))
      Rcout << "tree size: " << treesize() << endl;
   //Rcout << pad << "(id,parent): " << id << sp << pid;
   Rcout << pad << "id: " << id;
   Rcout << sp << "(v,c): " << v << sp << c;
   Rcout << sp << "mu: " << mu;
   Rcout << sp << "type: " << ntype();
   Rcout << sp << "depth: " << depth();
   //Rcout << sp << "pointer: " << this << endl;
   Rcout << endl;

   if(pc) {
      if(l) {
         l->pr(pc);
         r->pr(pc);
      }
   }
}
//--------------------------------------------------
//private functions
//--------------------
//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
   if(n->l) {
      Rcout << "cp:error node has children\n";
      return;
   }

   n->mu = o->mu;
   n->s = o->s;
   n->v = o->v;
   n->c = o->c;

   if(o->l) { //if o has children
      n->l = new tree;
      (n->l)->p = n;
      cp(n->l,o->l);
      n->r = new tree;
      (n->r)->p = n;
      cp(n->r,o->r);
   }
}
//--------------------
//cut back to one node
void tree::tonull()
{
   size_t ts = treesize();
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = treesize();
   }
   mu=zeros(1);    // UPDATED, WAS: mu=0.0;
   v=0;c=0;
   p=0;l=0;r=0;
}
//--------------------------------------------------
//functions
//--------------------

//flatten

/*
Rcpp::List tree::flatten(double scale=1.0) 
{
  tree::cnpv nds;
  this->getnodes(nds);
  Rcpp::List out(nds.size());
  
  arma::vec mu_null;
  arma::vec mu_save;
  mu_null.resize(0);
  for(size_t i=0;i<nds.size();i++) {
    //os << nds[i]->nid() << " ";
    //os << nds[i]->getv() << " ";
    //os << nds[i]->getc() << " ";
    
    if(nds[i]->getl() == 0) {
      mu_save = nds[i]->getm()*scale;
    } else {
      mu_save = mu_null;
    }
    
    out[i] = Rcpp::List::create(Rcpp::Named("nid") = nds[i]->nid(),
                                Rcpp::Named("v") = nds[i]->getv(),
                                Rcpp::Named("c") = nds[i]->getc(),
                                Rcpp::Named("m") = mu_save
                                );
  }
  return out;
}
*/
Rcpp::NumericMatrix tree::flatten_mat(double scale=1.0) 
{ 
  tree::cnpv nds;
  this->getnodes(nds);
  Rcpp::NumericMatrix out(nds.size(), 4);
  
  for(size_t i=0;i<nds.size();i++) {
    
    out(i,1) = nds[i]->nid();
    out(i,2) = nds[i]->getv();
    out(i,3) = nds[i]->getc();
    arma::vec m = nds[i]->getm();
    out(i,4) = m[0];
    //os << nds[i]->nid() << " ";
    //os << nds[i]->getv() << " ";
    //os << nds[i]->getc() << " ";
    
  }
  return out;
}


flat_tree tree::flatten() const 
{ 
  tree::cnpv nds;
  this->getnodes(nds);
  flat_tree out(nds.size());
  
  for(size_t i=0;i<nds.size();i++) {
    
    out[i].id = nds[i]->nid();
    out[i].v = nds[i]->getv();
    out[i].c = nds[i]->getc();
    out[i].m = arma::conv_to< std::vector<double> >::from(nds[i]->getm());
    //os << nds[i]->nid() << " ";
    //os << nds[i]->getv() << " ";
    //os << nds[i]->getc() << " ";
    
  }
  return out;
}

void unflatten(Rcpp::List& flat_tree, tree& t) {
//  Rcout<< "unflattening tree with " << flat_tree.size() << " nodes "<< endl;
  size_t tid,pid; //tid: id of current node, pid: parent's id
  std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
  size_t nn = flat_tree.size(); //number of nodes
  
  t.tonull(); // obliterate old tree (if there)
  
  std::vector<node_info> nv(nn);   // Vector of node_info objects.
  
  arma::vec mu_aux;
  arma::vec mu_null;
  mu_null.resize(0);
  for(size_t i=0;i!=nn;i++) {
    //int b_size; //basis size
    //is >> nv[i].id >> nv[i].v >> nv[i].c >> b_size;
    
    List ft = flat_tree[i];
    //Rcout << "ft size " << ft.size();
    
    nv[i].id = ft["nid"];
    nv[i].v = ft["v"];
    nv[i].c = ft["c"];
    nv[i].m = as<arma::vec>(ft["m"]);

  }
  
  //first node has to be the top one
  pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
  t.setv(nv[0].v); t.setc(nv[0].c); t.setm(nv[0].m);
  t.p=0;
  
  //now loop through the rest of the nodes knowing parent is already there.
  for(size_t i=1;i!=nv.size();i++) {
    tree::tree_p np = new tree;
    np->v = nv[i].v; np->c=nv[i].c; np->mu=nv[i].m;
    tid = nv[i].id;
    pts[tid] = np;
    pid = tid/2;
    // set pointers
    if(tid % 2 == 0) { //left child has even id
      pts[pid]->l = np;
    } else {
      pts[pid]->r = np;
    }
    np->p = pts[pid];
  }
  return;
}


// [[Rcpp::export]]
void unflatten_test(Rcpp::List flat_tree) {
  tree t;
  unflatten(flat_tree, t);
  Rcout << t;
}
// [[Rcpp::export]]
void unflatten_test_predict(Rcpp::List flat_tree, List x_info_list) {
  
  
  int p = x_info_list.size();
  xinfo xi;
  xi.resize(p);
  for(int j=0; j<p; ++j) {
    NumericVector tmp = x_info_list[j];
    std::vector<double> tmp2;
    for(size_t s=0; s<tmp.size(); ++s) {
      tmp2.push_back(tmp[s]);
    }
    xi[j] = tmp2;
  }
  
  tree t;
  unflatten(flat_tree, t);
  Rcout << t;
  
  
}

//output operator
std::ostream& operator<<(std::ostream& os, const tree& t)
{
   tree::cnpv nds;
   t.getnodes(nds);
   os << nds.size() << endl;
   for(size_t i=0;i<nds.size();i++) {
      os << nds[i]->nid() << " ";
      os << nds[i]->getv() << " ";
      os << nds[i]->getc() << " ";
      
      if(nds[i]->getl() == 0) {
        os << (nds[i]->getm()).size();
        for(size_t hh=0; hh < (nds[i]->getm()).size(); hh++){
          os << " " << (nds[i]->getm())[hh]; // All the parameters in the leaves
        }
      } else {
        os << -1;
      }
      os << endl;
   }
   return os;
}
//--------------------
//input operator
std::istream& operator>>(std::istream& is, tree& t)  //UPDATE: Added double T to inputs for streaming in means.
{
   size_t tid,pid; //tid: id of current node, pid: parent's id
   std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
   size_t nn; //number of nodes

   t.tonull(); // obliterate old tree (if there)

   //read number of nodes----------
   is >> nn;
   
   //Rcout << "tree has " << nn << " nodes" << endl;
   if(!is) {
      //cout << ">> error: unable to read number of nodes" << endl;
      return is;
   }

   // //read in vector of node information----------
   // std::vector<node_info> nv(nn);   // Vector of node_info objects.
   // for(size_t i=0;i!=nn;i++) {
   //    is >> nv[i].id >> nv[i].v >> nv[i].c >> nv[i].m;
   //    if(!is) {
   //       //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
   //       return is;
   //    }
   // }

   //read in vector of node information----------
   std::vector<node_info> nv(nn);   // Vector of node_info objects.
  
   arma::vec mu_aux;
   arma::vec mu_null;
   mu_null.resize(0);
   for(size_t i=0;i!=nn;i++) {
     int b_size; //basis size
     is >> nv[i].id >> nv[i].v >> nv[i].c >> b_size;
     
     //Rcout << endl << nv[i].id << " " << nv[i].v<< " " << nv[i].c<< " " << b_size << endl;
     
     if(b_size>0) {
       mu_aux.set_size(b_size);
       for(size_t dd = 0; dd<b_size; dd++){
         is >> mu_aux[dd];
       }
       nv[i].m = mu_aux;
     } else {
       nv[i].m = mu_null;
     }

      if(!is) {
         //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
         return is;
      }
   }

   //first node has to be the top one
   pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
   t.setv(nv[0].v); t.setc(nv[0].c); t.setm(nv[0].m);
   t.p=0;

   //now loop through the rest of the nodes knowing parent is already there.
   for(size_t i=1;i!=nv.size();i++) {
      tree::tree_p np = new tree;
      np->v = nv[i].v; np->c=nv[i].c; np->mu=nv[i].m;
      tid = nv[i].id;
      pts[tid] = np;
      pid = tid/2;
      // set pointers
      if(tid % 2 == 0) { //left child has even id
         pts[pid]->l = np;
      } else {
         pts[pid]->r = np;
      }
      np->p = pts[pid];
   }
   return is;
}
std::ostream& operator<<(std::ostream& os, const xinfo& xi)
{
	os << xi.size() << endl;
	for(size_t i=0;i<xi.size();i++) {
		os << xi[i].size() << endl;
		for(size_t j=0;j<xi[i].size();j++)
			os << xi[i][j]<<endl;
		os << endl;
	}

	return os;
}
std::istream& operator>>(std::istream& is, xinfo& xi)
{
	size_t xin;
	size_t vecdn;

	xi.resize(0); // reset old xinfo (if there)

	is >> xin;
	if(!is) {
		//cout << ">> error: unable to read size of xinfo" << endl;
		return is;
	}

	std::vector<double> vec_d;
	double vecdelem;

	for(size_t i=0;i<xin;i++) {
		is >> vecdn;
		for(size_t j=0;j<vecdn;j++) {
			is >> vecdelem;
			vec_d.push_back(vecdelem);
		}
		xi.push_back(vec_d);
		vec_d.resize(0);
	}

	return is;
}


//--------------------
// get pointer for node from its nid
tree::tree_p tree::getptr(size_t nid)
{
   if(this->nid() == nid) return this; //found it
   if(l==0) return 0; //no children, did not find it
   tree_p lp = l->getptr(nid);
   if(lp) return lp; //found on left
   tree_p rp = r->getptr(nid);
   if(rp) return rp; //found on right
   return 0; //never found it
}

//--------------------
// remove parameter vectors from interior nodes
void tree::compress()
{
  npv interior_nodes;
  this->getnobots(interior_nodes);
  arma::vec m; m.resize(0);
  for(size_t b=0; b<interior_nodes.size(); ++b) {
    interior_nodes[b]->setm(m);
  }
}

//--------------------
// scale parameter vectors in end nodes
void tree::scale(double s)
{
  npv end_nodes;
  this->getbots(end_nodes);
  for(size_t b=0; b<end_nodes.size(); ++b) {
    end_nodes[b]->setm(end_nodes[b]->getm()*s);
  }
}
