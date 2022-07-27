#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rho_update(arma::mat X, arma::vec E_logp, arma::cube E_residuals, double E_tau)
{
  
  // Initialization
  int n = X.n_rows;
  int T = X.n_cols;
  int H = E_logp.n_elem;
  
  arma::mat rho(n,H); 
  arma::vec lprob(H); // Log-probabilities
  arma::vec prob(H); // Probabilities
  
  // Auxiliary quantities
  
  for(int i=0; i < n; i++) {
    for(int h =0; h < H; h++) {
      lprob[h] = E_logp(h);
      for(int t=0; t < T; t++) {
        if(!NumericVector::is_na(X(i,t))) {
          lprob(h) = lprob(h) - 0.5*E_tau*E_residuals(h, i, t);
        }
      }
    }
    lprob = lprob - max(lprob);
    prob  = exp(lprob);
    prob  = prob/sum(prob);
    rho.row(i) = prob.t();
  }
  return(rho);
}