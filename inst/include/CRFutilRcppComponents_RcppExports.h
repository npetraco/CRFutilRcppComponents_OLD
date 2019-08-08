// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_CRFutilRcppComponents_RCPPEXPORTS_H_GEN_
#define RCPP_CRFutilRcppComponents_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace CRFutilRcppComponents {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("CRFutilRcppComponents", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("CRFutilRcppComponents", "_CRFutilRcppComponents_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in CRFutilRcppComponents");
            }
        }
    }

    inline int from_head_test() {
        typedef SEXP(*Ptr_from_head_test)();
        static Ptr_from_head_test p_from_head_test = NULL;
        if (p_from_head_test == NULL) {
            validateSignature("int(*from_head_test)()");
            p_from_head_test = (Ptr_from_head_test)R_GetCCallable("CRFutilRcppComponents", "_CRFutilRcppComponents_from_head_test");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_from_head_test();
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<int >(rcpp_result_gen);
    }

    inline arma::Mat<int> ff_C(int x) {
        typedef SEXP(*Ptr_ff_C)(SEXP);
        static Ptr_ff_C p_ff_C = NULL;
        if (p_ff_C == NULL) {
            validateSignature("arma::Mat<int>(*ff_C)(int)");
            p_ff_C = (Ptr_ff_C)R_GetCCallable("CRFutilRcppComponents", "_CRFutilRcppComponents_ff_C");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ff_C(Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::Mat<int> >(rcpp_result_gen);
    }

    inline List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par) {
        typedef SEXP(*Ptr_fix_node_and_edge_par)(SEXP,SEXP);
        static Ptr_fix_node_and_edge_par p_fix_node_and_edge_par = NULL;
        if (p_fix_node_and_edge_par == NULL) {
            validateSignature("List(*fix_node_and_edge_par)(arma::Cube<int>,List)");
            p_fix_node_and_edge_par = (Ptr_fix_node_and_edge_par)R_GetCCallable("CRFutilRcppComponents", "_CRFutilRcppComponents_fix_node_and_edge_par");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_fix_node_and_edge_par(Shield<SEXP>(Rcpp::wrap(node_par)), Shield<SEXP>(Rcpp::wrap(edge_par)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline arma::uvec row_match(arma::Mat<int> x, arma::Mat<int> table) {
        typedef SEXP(*Ptr_row_match)(SEXP,SEXP);
        static Ptr_row_match p_row_match = NULL;
        if (p_row_match == NULL) {
            validateSignature("arma::uvec(*row_match)(arma::Mat<int>,arma::Mat<int>)");
            p_row_match = (Ptr_row_match)R_GetCCallable("CRFutilRcppComponents", "_CRFutilRcppComponents_row_match");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_row_match(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(table)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::uvec >(rcpp_result_gen);
    }

}

#endif // RCPP_CRFutilRcppComponents_RCPPEXPORTS_H_GEN_
