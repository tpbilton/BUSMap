// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MH_Bayes_Hir
NumericMatrix MH_Bayes_Hir(IntegerMatrix ref, IntegerMatrix alt, IntegerVector OPGP, int nInd, int nSnps, NumericVector startVal, IntegerVector simPar, unsigned seed);
RcppExport SEXP _BUSMap_MH_Bayes_Hir(SEXP refSEXP, SEXP altSEXP, SEXP OPGPSEXP, SEXP nIndSEXP, SEXP nSnpsSEXP, SEXP startValSEXP, SEXP simParSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type alt(altSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type OPGP(OPGPSEXP);
    Rcpp::traits::input_parameter< int >::type nInd(nIndSEXP);
    Rcpp::traits::input_parameter< int >::type nSnps(nSnpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startVal(startValSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type simPar(simParSEXP);
    Rcpp::traits::input_parameter< unsigned >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_Bayes_Hir(ref, alt, OPGP, nInd, nSnps, startVal, simPar, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BUSMap_MH_Bayes_Hir", (DL_FUNC) &_BUSMap_MH_Bayes_Hir, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_BUSMap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}