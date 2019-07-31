#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace std;

arma::Col<int> ff_C(int x);
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par);

// Two state spinor function
arma::Col<int> ff_C(int x){
  arma::Col<int> aspinor(2);
  
  aspinor[0] = (x == 1);
  aspinor[1] = (x == 2);
  
  return(aspinor);
  
}

// CRF default node_par and edge_par are 3D "Cubes" with 1 slice.
// This makes them arma integer matrices to head of annoying type 
// conflicts and code bloat.
// [[Rcpp::export]]
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par){
  
  // Remove third index from node_par
  arma::Mat<int> node_par_new(node_par.n_rows,2);
  node_par_new = node_par.slice(0);
  //Rcout << node_par_new << endl;
  
  // Remove third index from each of the edge_par
  arma::Mat<int> aedge_par_mat(2,2);
  IntegerVector tmpe(4);
  List edge_par_new(edge_par.size());
  for(int i=0; i<edge_par.size(); ++i) {
    
    tmpe = as<IntegerVector>(edge_par(i));
    aedge_par_mat(0,0) = tmpe(0);
    aedge_par_mat(0,1) = tmpe(1);
    aedge_par_mat(1,0) = tmpe(2);
    aedge_par_mat(1,1) = tmpe(3);
    
    edge_par_new(i) = aedge_par_mat;
    
    //Rcout << aedge_par_mat << endl;
    
  }
  
  List theta_pars;
  theta_pars["node_par"] = node_par_new;
  theta_pars["edge_par"] = edge_par_new;
  
  return theta_pars;
  
}
