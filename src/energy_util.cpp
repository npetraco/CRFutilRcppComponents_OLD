#include "RcppArmadillo.h"
#include "CRFutilRcppComponents.h"

using namespace Rcpp;
using namespace CRFutilRcppComponents;
using namespace std;

// [[Rcpp::export]]
int symbolic_conditional_energy_C(
    arma::Mat<int>                config, 
    Rcpp::Nullable<int>           condition_element_number, 
    Rcpp::Nullable<IntegerMatrix> node_par,
    Rcpp::Nullable<List>          edge_par,
    Rcpp::Nullable<IntegerMatrix> edge_mat) {
  
  //arma::Mat<int> aphi = phi_features_C(config, edge_mat, node_par, edge_par);

  from_head_test();
  
  int y = 1;
  return y;
  //config, condition.element.number, crf, ff, format="tex", printQ=FALSE
  
}