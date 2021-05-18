#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>      // std::iota

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Recursively forecasting a VAR
// Warning: returns in companion form. e.g. only first k columns are of interest
// [[Rcpp::export()]]
arma::mat recursiveforecast_cpp(arma::vec &inity, arma::mat &FF, arma::vec &c, int h){
  arma::mat fcst(h+1, FF.n_cols);
  fcst.row(0) = inity.t();

  for (int hh=1; hh<=h; ++hh){
    arma::vec y_hh = c + FF*fcst.row(hh-1).t();
    fcst.row(hh) = y_hh.t();
  }

  return(fcst);
}
