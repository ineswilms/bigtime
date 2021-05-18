#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>      // std::iota

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// simulating a VAR
// to be called from simVAR.R
// [[Rcpp::export]]
arma::mat simVAR_cpp(int periods, int k, int p, arma::mat &coef_mat, 
                     arma::vec &c, arma::mat &e, arma::vec &init_y, int burnin){
  
  
  // Setting up the simulation
  int tot_length = periods + burnin;
  arma::mat Y(tot_length+1, k*p);
  Y.row(0) = init_y.t();
  
  // Running the recursive loop
  // Note that we use the companion form throughout
  for (int t=1; t<tot_length; ++t){
    arma::vec y = c + coef_mat*Y.row(t-1).t() + e.col(t);
    Y.row(t) = y.t();
  }
  
  return(Y.rows(tot_length-periods, Y.n_rows-2).cols(0, k-1));
}







