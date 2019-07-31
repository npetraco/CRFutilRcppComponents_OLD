#include "RcppArmadillo.h"
#include "loc_util.hpp"

using namespace Rcpp;
using namespace std;

// //[[Rcpp::export]]
// arma::Col<int> ff_C(int x){
//   arma::Col<int> aspinor(2);
//   
//   //Rcout << x << endl;
//   aspinor[0] = (x == 1);
//   aspinor[1] = (x == 2);
//   
//   return(aspinor);
//   
// }
//Place dependent functions above call or//declare the function definition with:
//arma::Col<int> ff_C(int x);


// [[Rcpp::export]]
arma::Row<int> phi_features_C(arma::Col<int> config, arma::mat edge_mat, arma::Mat<int> node_par, List edge_par) {
  
  int num_nodes = config.size();
  int num_edges = edge_mat.n_rows;

  // Find num_params buy digging out max in par matrices. 
  // Stupid not to just send this in, but done to keep 
  // parity with the function signature of the original R code.....  
  // MOVE THIS OUT LATER SO WE DON'T HAVE TO DO IT OVER AND OVER AGAIN IF WE USE THIS FUNCTION IN A LOOP
  int num_params = node_par.max();
  int amax       = 0;
  for(int i=0; i<edge_par.size(); ++i) {
    amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
    if(amax > num_params){
      num_params = amax;
    }
  }
  
  // Initialize a phi vector:
  arma::Row<int> phi_vec(num_params);
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

  return phi_vec;
  
}

