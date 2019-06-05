#include <RcppArmadillo.h>
#include "probabilitytracking.h"
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)  = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  return(out);
}



NumericVector logprob(NumericVector deviationsx, NumericVector deviationsy, NumericVector px, NumericVector py, NumericVector kernelx, NumericVector kernely, double sigma, int Kx, int Ky, int ngrid, NumericVector w, NumericVector w_indices, arma::vec logw, NumericVector U, double q, NumericMatrix Z, NumericMatrix contacts, NumericVector naccidentals, int n, int nshoes, NumericVector logPhi, NumericVector Phi, double qa, double qb, arma::mat wprecision){
  NumericVector loglike(5);
  //likelihood of accidentals given locations

  for(int k = 0; k < (Kx + 1); ++k){
    loglike(0) += deviationsx[k] * log(kernelx[k]);
  }

  for(int k = 0; k < (Ky + 1); ++k){
    loglike(0) += deviationsy[k] * log(kernely[k]);
  }

  // probability of locations and u
  for(int d = 0; d < nshoes; ++d){
    loglike(1) = loglike(1) + log(U[d]) * (naccidentals[d] - 1) - ngrid * lgamma(q) - lgamma(naccidentals[d]);
    for(int a = 0; a < ngrid; ++a){
      loglike(1) = loglike(1) + Z(d, a) *(logw[w_indices[a]] + logPhi[contacts(d, a)]) - (Z(d, a) + q) * log(U[d] * w[w_indices[a]] * Phi[contacts(d, a)] + 1) + lgamma(Z(d, a) + q);
    }
  }
  // probability of ws
  loglike(2) = -as_scalar(logw.t() * wprecision * logw);


  // probability of q
  loglike(3) = loglike(3) + (qa - 1) * log(q) - (qb) * q;

  // probability of kernels;
  for(int i = 0; i < Kx; ++i){
    loglike(4) += px[i] * px[i];
  }
  for(int i = 0; i < Ky; ++i){
    loglike(4) += py[i] * py[i];
  }
  // pi is 2 *acos(0.0). brought that 2 into the other to omit the square root
  loglike(4) = -loglike(4) * 0.5 * 1.0/(sigma * sigma) - (Kx + Ky) * (log(acos(0.0))/2.0 - log(sigma * 2));

  return(loglike);

}

NumericVector logprob_noZ(int ngrid, NumericVector w, NumericVector w_indices, arma::vec logw, NumericVector U, double q, NumericMatrix Z, NumericMatrix contacts, NumericVector naccidentals, int n, int nshoes, NumericVector logPhi, NumericVector Phi, double qa, double qb, arma::mat wprecision){
  NumericVector loglike(3);
  // probability of locations and u
  for(int d = 0; d < nshoes; ++d){
    loglike(0) = loglike(0) + log(U[d]) * (naccidentals[d] - 1) - ngrid * lgamma(q) - lgamma(naccidentals[d]);
    for(int a = 0; a < ngrid; ++a){
      loglike(0) = loglike(0) + Z(d, a) *(logw[w_indices[a]] + logPhi[contacts(d, a)]) - (Z(d, a) + q) * log(U[d] * w[w_indices[a]] * Phi[contacts(d, a)] + 1) + lgamma(Z(d, a) + q);
    }
  }
  // probability of ws
  loglike(1) = -as_scalar(logw.t() * wprecision * logw);


  // probability of q
  loglike(2) = loglike(2) + (qa - 1) * log(q) - (qb) * q;

  return(loglike);
}

NumericVector logprob_now(NumericVector deviationsx, NumericVector deviationsy, NumericVector px, NumericVector py, NumericVector kernelx, NumericVector kernely, double sigma, int Kx, int Ky, int ngrid, NumericVector U, double q, NumericMatrix Z, NumericMatrix contacts, NumericVector naccidentals, int n, int nshoes, NumericVector logPhi, NumericVector Phi, double qa, double qb){
  NumericVector loglike(4);
  //likelihood of accidentals given locations

  for(int k = 0; k < (Kx + 1); ++k){
    loglike(0) += deviationsx[k] * log(kernelx[k]);
  }

  for(int k = 0; k < (Ky + 1); ++k){
    loglike(0) += deviationsy[k] * log(kernely[k]);
  }

  // probability of locations and u
  for(int d = 0; d < nshoes; ++d){
    loglike(1) = loglike(1) + log(U[d]) * (naccidentals[d] - 1) - ngrid * lgamma(q) - lgamma(naccidentals[d]);
    for(int a = 0; a < ngrid; ++a){
      loglike(1) = loglike(1) + Z(d, a) *(logPhi[contacts(d, a)]) - (Z(d, a) + q) * log(U[d] * Phi[contacts(d, a)] + 1) + lgamma(Z(d, a) + q);
    }
  }


  // probability of q
  loglike(2) = loglike(2) + (qa - 1) * log(q) - (qb) * q;

  // probability of kernels;
  for(int i = 0; i < Kx; ++i){
    loglike(3) += px[i] * px[i];
  }
  for(int i = 0; i < Ky; ++i){
    loglike(3) += py[i] * py[i];
  }
  // pi is 2 *acos(0.0). brought that 2 into the other to omit the square root
  loglike(3) = -loglike(3) * 0.5 * 1.0/(sigma * sigma) - (Kx + Ky) * (log(acos(0.0))/2.0 - log(sigma * 2));

  return(loglike);

}


NumericVector logprob_noepsilon(NumericVector deviationsx, NumericVector deviationsy, NumericVector px, NumericVector py, NumericVector kernelx, NumericVector kernely, double sigma, int Kx, int Ky, int ngrid, NumericVector w, NumericVector w_indices, arma::vec logw, NumericVector zhome, NumericVector z, NumericMatrix zoptions, NumericVector shoenums, NumericVector denoms, NumericMatrix contacts, NumericVector naccidentals, int n, int nshoes, NumericVector logPhi, NumericVector Phi, arma::mat wprecision){
  NumericVector loglike(4);
  //likelihood of accidentals given locations

  for(int k = 0; k < (Kx + 1); ++k){
    loglike(0) += deviationsx[k] * log(kernelx[k]);
  }

  for(int k = 0; k < (Ky + 1); ++k){
    loglike(0) += deviationsy[k] * log(kernely[k]);
  }

  // probability of locations
  for(int i = 0; i < n; ++i){
    loglike(1) = loglike(1) + logw[w_indices[zoptions(zhome[i], z[i])]] + logPhi[contacts(shoenums[i],zoptions(zhome[i], z[i]))] - log(denoms[shoenums[i]]);
  }
  // probability of ws
  loglike(2) = -as_scalar(logw.t() * wprecision * logw);

  // probability of kernels;
  for(int i = 0; i < Kx; ++i){
    loglike(3) += px[i] * px[i];
  }
  for(int i = 0; i < Ky; ++i){
    loglike(3) += py[i] * py[i];
  }
  // pi is 2 *acos(0.0). brought that 2 into the other to omit the square root
  loglike(3) = -loglike(3) * 0.5 * 1.0/(sigma * sigma) - (Kx + Ky) * (log(acos(0.0))/2.0 - log(sigma * 2));

  return(loglike);

}


NumericVector logprob_noepsilonnoZ(NumericVector w_indices, arma::vec logw, NumericVector zhome, NumericVector shoenums, NumericVector denoms, NumericMatrix contacts, NumericVector naccidentals, int n, NumericVector logPhi, NumericVector Phi, arma::mat wprecision){
  NumericVector loglike(2);
  //likelihood of accidentals given locations
  // probability of locations
  for(int i = 0; i < n; ++i){
    loglike(0) = loglike(1) + logw[w_indices[zhome[i]]] + logPhi[contacts(shoenums[i],zhome[i])] - log(denoms[shoenums[i]]);
  }
  // probability of ws
  loglike(1) = -as_scalar(logw.t() * wprecision * logw);
  return(loglike);

}

NumericVector logprob_noPhi(NumericVector deviationsx, NumericVector deviationsy, NumericVector px, NumericVector py, NumericVector kernelx, NumericVector kernely, double sigma, int Kx, int Ky, int ngrid, NumericVector w, NumericVector w_indices, arma::vec logw, NumericVector U, double q, NumericMatrix Z, NumericVector naccidentals, int n, int nshoes, double qa, double qb, arma::mat wprecision){
  NumericVector loglike(5);
  //likelihood of accidentals given locations

  for(int k = 0; k < (Kx + 1); ++k){
    loglike(0) += deviationsx[k] * log(kernelx[k]);
  }

  for(int k = 0; k < (Ky + 1); ++k){
    loglike(0) += deviationsy[k] * log(kernely[k]);
  }

  // probability of locations and u
  for(int d = 0; d < nshoes; ++d){
    loglike(1) = loglike(1) + log(U[d]) * (naccidentals[d] - 1) - ngrid * lgamma(q) - lgamma(naccidentals[d]);
    for(int a = 0; a < ngrid; ++a){
      loglike(1) = loglike(1) + Z(d, a) *(logw[w_indices[a]]) - (Z(d, a) + q) * log(U[d] * w[w_indices[a]] + 1) + lgamma(Z(d, a) + q);
    }
  }
  // probability of ws
  loglike(2) = -as_scalar(logw.t() * wprecision * logw);


  // probability of q
  loglike(3) = loglike(3) + (qa - 1) * log(q) - (qb) * q;

  // probability of kernels;
  for(int i = 0; i < Kx; ++i){
    loglike(4) += px[i] * px[i];
  }
  for(int i = 0; i < Ky; ++i){
    loglike(4) += py[i] * py[i];
  }
  // pi is 2 *acos(0.0). brought that 2 into the other to omit the square root
  loglike(4) = -loglike(4) * 0.5 * 1.0/(sigma * sigma) - (Kx + Ky) * (log(acos(0.0))/2.0 - log(sigma * 2));

  return(loglike);

}



