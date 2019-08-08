#include "RcppArmadillo.h"

#ifndef LOCUTIL_H
#define LOCUTIL_H

using namespace Rcpp;
using namespace std;

arma::Mat<int> ff_C(int x);
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par);
List fix_node_and_edge_par2(arma::Cube<int> node_par, List edge_par);

// Two state spinor function for internal use. We Don't want to pass in R functions so use the C version
arma::Mat<int> ff_C(int x){
  arma::Mat<int> aspinor(2,1);
  
  aspinor[0] = (x == 1);
  aspinor[1] = (x == 2);
  
  return(aspinor);
  
}

// CRF default node_par and edge_par are 3D "Cubes" with 1 slice (3 index R arrays).
// This function makes them armadillo integer matrices (2D) to head of annoying type 
// conflicts and code bloat from re-use. For now, we export to R as well for testing.
// Ultimately intended for internal use.
// Armadillo version
// [[Rcpp::export]]
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par){
  
  arma::Mat<int> node_par_new(node_par.n_rows,2);
  arma::Mat<int> aedge_par_mat(2,2);
  
  node_par_new = node_par.slice(0);

  // Remove third index from each of the edge_par
  IntegerVector tmpe(4);
  List edge_par_new(edge_par.size());
  for(int i=0; i<edge_par.size(); ++i) {
    
    tmpe = as<IntegerVector>(edge_par(i));
    aedge_par_mat(0,0) = tmpe(0);
    aedge_par_mat(0,1) = tmpe(1);
    aedge_par_mat(1,0) = tmpe(2);
    aedge_par_mat(1,1) = tmpe(3);
    
    edge_par_new(i) = aedge_par_mat;
    
  }
  
  List theta_pars;
  theta_pars["node_par"] = node_par_new;
  theta_pars["edge_par"] = edge_par_new;
  
  return theta_pars;
  
}

// Rcpp version. **********GET RID OF THIS IF WE NEVER USE IT
// [[Rcpp::export]]
List fix_node_and_edge_par2(arma::Cube<int> node_par, List edge_par){
  
  IntegerMatrix node_par_new(node_par.n_rows,2);
  IntegerMatrix aedge_par_mat(2,2);
  
  node_par_new = wrap(node_par.slice(0));
  
  // Remove third index from each of the edge_par
  IntegerVector tmpe(4);
  List edge_par_new(edge_par.size());
  for(int i=0; i<edge_par.size(); ++i) {
    
    tmpe = as<IntegerVector>(edge_par(i));
    aedge_par_mat(0,0) = tmpe(0);
    aedge_par_mat(0,1) = tmpe(1);
    aedge_par_mat(1,0) = tmpe(2);
    aedge_par_mat(1,1) = tmpe(3);
    
    edge_par_new(i) = aedge_par_mat;
    
  }
  
  List theta_pars;
  theta_pars["node_par"] = node_par_new;
  theta_pars["edge_par"] = edge_par_new;
  
  return theta_pars;
  
}

// [[Rcpp::export]]
arma::uvec row_match(arma::Mat<int> x, arma::Mat<int> table){
  
  arma::uvec matching_row_idxs;
  
  // Note: Below assume x is a effectively a vector with one row or one column
  // We pass it in as a matrix though. Noticed Armadillo routenes often ran faster that way
  
  // Get common first set of possible row match indices   
  matching_row_idxs = intersect(find(table.col(0) == x(0)), find(table.col(1) == x(1)));
    
  // Look over the rest and update (pairdown) the possible matching rows 
  for(int i=2; i<x.size(); ++i){
    matching_row_idxs = intersect(matching_row_idxs, find(table.col(i) == x(i)));
  }

  // Note: a size 0 return means no matching row was found

  return matching_row_idxs;
}

#endif