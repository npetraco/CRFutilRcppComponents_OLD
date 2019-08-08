#include "RcppArmadillo.h"

// [[Rcpp::interfaces(r,cpp)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int from_head_test() {
  
  Rcout << "This work??" << endl;
  
  return 0;
}
  