// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// phi_features_C
arma::Mat<int> phi_features_C(arma::Mat<int> config, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default);
RcppExport SEXP _CRFutilRcppComponents_phi_features_C(SEXP configSEXP, SEXP edge_matSEXP, SEXP node_parSEXP, SEXP edge_parSEXP, SEXP num_params_defaultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    Rcpp::traits::input_parameter< int >::type num_params_default(num_params_defaultSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_features_C(config, edge_mat, node_par, edge_par, num_params_default));
    return rcpp_result_gen;
END_RCPP
}
// compute_model_matrix
arma::Mat<int> compute_model_matrix(arma::Mat<int> configs, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default);
RcppExport SEXP _CRFutilRcppComponents_compute_model_matrix(SEXP configsSEXP, SEXP edge_matSEXP, SEXP node_parSEXP, SEXP edge_parSEXP, SEXP num_params_defaultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type configs(configsSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    Rcpp::traits::input_parameter< int >::type num_params_default(num_params_defaultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_model_matrix(configs, edge_mat, node_par, edge_par, num_params_default));
    return rcpp_result_gen;
END_RCPP
}
// get_par_idx
int get_par_idx(arma::Mat<int> config, Rcpp::Nullable<int> i, Rcpp::Nullable<int> j, Rcpp::Nullable<IntegerMatrix> node_par_in);
RcppExport SEXP _CRFutilRcppComponents_get_par_idx(SEXP configSEXP, SEXP iSEXP, SEXP jSEXP, SEXP node_par_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type j(jSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<IntegerMatrix> >::type node_par_in(node_par_inSEXP);
    rcpp_result_gen = Rcpp::wrap(get_par_idx(config, i, j, node_par_in));
    return rcpp_result_gen;
END_RCPP
}
// fix_node_and_edge_par
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par);
RcppExport SEXP _CRFutilRcppComponents_fix_node_and_edge_par(SEXP node_parSEXP, SEXP edge_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Cube<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    rcpp_result_gen = Rcpp::wrap(fix_node_and_edge_par(node_par, edge_par));
    return rcpp_result_gen;
END_RCPP
}
// fix_node_and_edge_par2
List fix_node_and_edge_par2(arma::Cube<int> node_par, List edge_par);
RcppExport SEXP _CRFutilRcppComponents_fix_node_and_edge_par2(SEXP node_parSEXP, SEXP edge_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Cube<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    rcpp_result_gen = Rcpp::wrap(fix_node_and_edge_par2(node_par, edge_par));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _CRFutilRcppComponents_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _CRFutilRcppComponents_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _CRFutilRcppComponents_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _CRFutilRcppComponents_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CRFutilRcppComponents_phi_features_C", (DL_FUNC) &_CRFutilRcppComponents_phi_features_C, 5},
    {"_CRFutilRcppComponents_compute_model_matrix", (DL_FUNC) &_CRFutilRcppComponents_compute_model_matrix, 5},
    {"_CRFutilRcppComponents_get_par_idx", (DL_FUNC) &_CRFutilRcppComponents_get_par_idx, 4},
    {"_CRFutilRcppComponents_fix_node_and_edge_par", (DL_FUNC) &_CRFutilRcppComponents_fix_node_and_edge_par, 2},
    {"_CRFutilRcppComponents_fix_node_and_edge_par2", (DL_FUNC) &_CRFutilRcppComponents_fix_node_and_edge_par2, 2},
    {"_CRFutilRcppComponents_rcpparma_hello_world", (DL_FUNC) &_CRFutilRcppComponents_rcpparma_hello_world, 0},
    {"_CRFutilRcppComponents_rcpparma_outerproduct", (DL_FUNC) &_CRFutilRcppComponents_rcpparma_outerproduct, 1},
    {"_CRFutilRcppComponents_rcpparma_innerproduct", (DL_FUNC) &_CRFutilRcppComponents_rcpparma_innerproduct, 1},
    {"_CRFutilRcppComponents_rcpparma_bothproducts", (DL_FUNC) &_CRFutilRcppComponents_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_CRFutilRcppComponents(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
