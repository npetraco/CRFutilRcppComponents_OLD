#include "RcppArmadillo.h"
#include "loc_util.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
arma::Mat<int> phi_features_C(arma::Mat<int> config, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default=0) {
  
  int num_nodes = config.size();
  int num_edges = edge_mat.n_rows;

  // The original R function determines this everytime it is called. This is stupid and will
  // slow things down. To keep parity with the R function signature, by default the number of
  // parameters will be calculated. However the user can also prespcify it and change the default
  // value from 0 so the loop below isn't excecuted over an over when calling the function
  // multiple times.
  int num_params;
  if(num_params_default == 0) {
    
    num_params = node_par.max();
    int amax       = 0;
    for(int i=0; i<edge_par.size(); ++i) {
      amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
      if(amax > num_params){
        num_params = amax;
      } else {
        num_params = num_params_default;
      }
    }
    
  } 
  
  // Initialize a phi (row) vector. Choose row vector. Thats what the R code in compute.model.matrix assumes.
  // We do this here to keep parity with the R code.
  arma::Mat<int> phi_vec(1,num_params);
  phi_vec.zeros();
  
  // Nodes: \phi_i({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \tau}_i, 0}
  int phi_off = -1;                                        // -1 bec we want offsets NOT indices
  for(int i=0; i<num_nodes; ++i) {
    
    // Usually reversed but these are rowvec . colvec. Doing it this way avoids a transpose 
    phi_off = dot( node_par.row(i), ff_C(config(i)) ) - 1; // -1 bec we want offsets NOT indices

    // Kronecker delta part:
    if(phi_off != -1) {
      phi_vec(phi_off) = 1; 
    }
    
  }

  // Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  arma::Mat<int> aepm;
  for(int i=0; i<num_edges; ++i) {

    // Check and see if we reached the end of phi. No point in doing the rest of the edges if we did:
    // NOTE: assumes parameters have unique consecutive indices, so no wierd parameterizations.
    // Sticking to "standard" or "flexible" should be safe.
    if(phi_off == num_params) {
      break;
    }

    aepm = as<arma::Mat<int>>(edge_par(i));
    int left_off  = edge_mat(i,0) - 1;                                    // offset NOT index, do -1
    int right_off = edge_mat(i,1) - 1;                                    // offset NOT index, do -1
    
    arma::Mat<int> tmp = ff_C(config(left_off)).t() * aepm * ff_C(config(right_off));
    
    phi_off = tmp(0,0) - 1;                                               // offset NOT index, do -1

    // Kronecker delta part:
    if(phi_off != -1) {
      phi_vec(phi_off) = 1; 
    }
    
  }

  return phi_vec;
  
}

// [[Rcpp::export]]
arma::Mat<int> compute_model_matrix(arma::Mat<int> configs, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default = 0) {
  
  int num_configs = configs.n_rows;
  int num_params;

  if(num_params_default == 0) {

    num_params = node_par.max();
    int amax       = 0;
    for(int i=0; i<edge_par.size(); ++i) {
      amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
      if(amax > num_params){
        num_params = amax;
      } else {
        num_params = num_params_default;
      }
    }

  }
  
  //Rcout << num_configs << endl;
  //Rcout << num_params << endl;
  
  arma::Mat<int> model_mat(num_configs, num_params);
  
  //arma::Mat<int> aphi(num_params,1);
  for(int i=0; i<num_configs; ++i){
    //Rcout << i << endl;
    //Rcout << configs.row(i) << endl;
    //aphi = phi_features_C(configs.row(i), edge_mat, node_par, edge_par);
    //model_mat.row(i) = aphi;
    //Rcout << aphi << endl;
    model_mat.row(i) = phi_features_C(configs.row(i), edge_mat, node_par, edge_par);
  }
  
  return(model_mat);
  
  
}
