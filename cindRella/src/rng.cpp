#include "rng.h"

mat rmvnorm_prec(int n, int d, mat &Phi);

double rgamma(const double &a, const double &b) {
  return R::rgamma(a, b);
}
double rbeta(const double &a, const double &b) {
  return R::rbeta(a, b);
}
double runif() {
  return R::runif(0.0, 1.0);
}
double rchisq(const double &df) {
  return Rf_rchisq(df);
}
double znorm() {
  return R::rnorm(0.0, 1.0);
}
double rnorm(const double &mu, const double &sd) {
  return R::rnorm(mu, sd);
}
double rlogstick(double &alpha) {
  return log(1-pow(runif(), 1.0/alpha));
}

double rlogbeta(const double &alpha, const double &beta) {
  return log(rbeta(alpha, beta));
}

double rloggamma(const double &alpha) {
  if(alpha<0.5) {
    return log(rgamma(alpha+1.0, 1.0)) + log(runif())/alpha;
    //return(log(2) + log(runif())/alpha + log(rgamma(alpha+1, 1.0)));
  } else {
    return log(rgamma(alpha, 1.0));
    //return log(rgamma(alpha, 2.0));
  }
}

// written by Neil based on Jared's code (0 mean)
mat rmvnorm_prec(int n, int d, mat &Phi) {
  RNGScope tmp;
  mat z(d, n);
  std::generate(z.begin(), z.end(), znorm);
  mat res = solve(trimatu(chol(Phi)), z);
  return(res);
}


// // // [[Rcpp::export]]
mat rwish_I_root(int p, double df) {
  mat T(p, p);
  for(int i=0; i<p; i++) {
    T(i,i) = sqrt(rchisq(df-i));//no +1 because of 0 indexing!
      for(int j=0; j<i; j++) {
        T(i,j) = znorm();//Rf_rnorm(0., 1.);
      }
  }
  return(trimatl(T));
}

// // [[Rcpp::export]]
mat rwish_root(mat L, double df) {
  mat T = rwish_I_root(L.n_rows, df);
  return(trimatl(L)*trimatl(T));
}

// // [[Rcpp::export]]
mat rwish(mat S, double df) {
  mat R = chol(S);
  mat out = rwish_root(R.t(), df);
  out = out*out.t();
  return(out);
}

// // [[Rcpp::export]]
mat cpp_rmvnorm_prec(int n, vec &mu, mat &Phi) {
  RNGScope tmp;
  //NumericVector z_ = rnorm(n*mu.size());
  //mat z(z_.begin(), mu.size(), n, false, false);
  mat z(mu.size(), n);
  std::generate(z.begin(), z.end(), znorm);
  mat res = solve(trimatu(chol(Phi)), z);
  res.each_col() += mu;
  return(res);
}

//sample from N(Phi^(-1)m, Phi^(-1))
// // [[Rcpp::export]]
mat rmvnorm_post(int n, vec &m, mat &Phi) {
  RNGScope tmp;
  mat R = chol(Phi);
  vec mu = solve(trimatu(R), solve(trimatl(R.t()), m));
  //NumericVector z_ = rnorm(n*mu.size());
  //mat z(z_.begin(), mu.size(), n, false, false);
  mat z(mu.size(), n);
  std::generate(z.begin(), z.end(), znorm);
  mat res = solve(trimatu(chol(Phi)), z);
  res.each_col() += mu;
  return(res);
}


