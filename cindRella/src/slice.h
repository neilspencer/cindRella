#ifndef slice_h
#define slice_h
#include <RcppArmadillo.h>
#include "common.h"
#include "rng.h"

class logdensity {
public:
  virtual double val(double y) = 0;
};

class logdensityellipse{
public:
  virtual double val(arma::vec y) = 0;
};

// define elliptical slice sampler for w
class ld_logw: public logdensityellipse {
public:
  NumericVector Phi;
  double q;
  NumericVector U;
  NumericMatrix Z;
  NumericMatrix contacts;
  int nshoes;
  int ngrid;
  NumericVector windices;

  double val(arma::vec logw) {
    double res = 0;
    for(int a = 0; a < ngrid; ++a){
      double coeff = 0;
      double result = 0;
      for(int d = 0; d < nshoes; ++d){
        coeff = coeff + Z(d, a);
        result = result + (q + Z(d, a))* (-log(U[d] * exp(logw[windices[a]]) * Phi[contacts(d, a)] + 1));
      }
      res =  res + result + coeff * logw[windices[a]];
    }

    return(res);
  }
  ld_logw(NumericVector Phi_, double q_, NumericVector U_, NumericMatrix Z_, NumericMatrix contacts_, int nshoes_, int ngrid_, NumericVector windices_){Phi = Phi_; q=q_; U=U_; Z = Z_; contacts = contacts_; nshoes = nshoes_; ngrid = ngrid_; windices = windices_;}
};

// define elliptical slice sampler for w when Phi is not used
class ld_logw_noPhi: public logdensityellipse {
public:
  double q;
  NumericVector U;
  NumericMatrix Z;
  int nshoes;
  int ngrid;
  NumericVector windices;

  double val(arma::vec logw) {
    double res = 0;
    for(int a = 0; a < ngrid; ++a){
      double coeff = 0;
      double result = 0;
      for(int d = 0; d < nshoes; ++d){
        coeff = coeff + Z(d, a);
        result = result + (q + Z(d, a))* (-log(U[d] * exp(logw[windices[a]]) + 1));
      }
      res =  res + result + coeff * logw[windices[a]];
    }

    return(res);
  }
  ld_logw_noPhi(double q_, NumericVector U_, NumericMatrix Z_, int nshoes_, int ngrid_, NumericVector windices_){q=q_; U=U_; Z = Z_; nshoes = nshoes_; ngrid = ngrid_; windices = windices_;}
};

// define elliptical slice sampler for w when epsilon terms are not used
class ld_logw_noepsilon: public logdensityellipse {
public:
  NumericVector Phi;
  NumericMatrix Z;
  NumericMatrix contacts;
  int nshoes;
  int ngrid;
  NumericVector naccidentals;
  NumericVector windices;

  double val(arma::vec logw) {
    double res = 0;
    arma::vec w = exp(logw);
    for(int d = 0; d < nshoes; ++d){
      double denom = 0;
      for(int a = 0; a < ngrid; ++a){
        res = res + Z(d, a) * logw[windices[a]];
        denom = denom + w[windices[a]] * Phi[contacts(d, a)];
      }
      res =  res - log(denom) * naccidentals[d];
    }

    return(res);
  }
  ld_logw_noepsilon(NumericVector Phi_, NumericMatrix Z_, NumericMatrix contacts_, int nshoes_, int ngrid_, NumericVector naccidentals_, NumericVector windices_){Phi = Phi_; Z = Z_; contacts = contacts_; nshoes = nshoes_; ngrid = ngrid_; naccidentals = naccidentals_; windices = windices_;}
};

// define slice sampler for u
class ld_u: public logdensity {
public:
  double q;
  NumericVector naccidentals;
  int d;
  NumericVector w;
  NumericMatrix Z;
  NumericMatrix contacts;
  NumericVector Phi;
  int ngrid;
  NumericVector w_indices;

  double val(double u) {
    double result = (naccidentals[d] - 1) * log(u);
    for(int a = 0; a < ngrid; ++a){
      result = result - (q + Z(d,a)) * log(u * w[w_indices[a]] * Phi[contacts(d,a)] + 1.0);
    }
    return(result);
  }
  ld_u(NumericVector Phi_, double q_, NumericVector naccidentals_, int d_, NumericVector w_, NumericMatrix Z_, NumericMatrix contacts_, int ngrid_, NumericVector w_indices_){ Phi = Phi_; q=q_; naccidentals = naccidentals_; d = d_; w = w_; Z = Z_; contacts = contacts_; ngrid = ngrid_; w_indices = w_indices_; }
};

// define slice sampler for u when phi is not used
class ld_u_noPhi: public logdensity {
public:
  double q;
  NumericVector naccidentals;
  int d;
  NumericVector w;
  NumericMatrix Z;
  int ngrid;
  NumericVector w_indices;

  double val(double u) {
    double result = (naccidentals[d] - 1) * log(u);
    for(int a = 0; a < ngrid; ++a){
      result = result - (q + Z(d,a)) * log(u * w[w_indices[a]] + 1.0);
    }
    return(result);
  }
  ld_u_noPhi(double q_, NumericVector naccidentals_, int d_, NumericVector w_, NumericMatrix Z_, int ngrid_, NumericVector w_indices_){ q=q_; naccidentals = naccidentals_; d = d_; w = w_; Z = Z_; ngrid = ngrid_; w_indices = w_indices_; }
};

// slice sampler for phi
class ld_phi: public logdensity {
public:
  int numb;
  double q;
  std::vector<NumericMatrix> philocs;
  NumericMatrix Z;
  NumericVector U;
  NumericVector w;
  NumericVector w_indices;
  double val(double phi){
    NumericMatrix phinumb = philocs[numb];
    double result = 0.0;
    for(int i=0; i < phinumb.nrow(); ++i){
      result = result + Z(phinumb(i, 0), phinumb(i, 1)) * log(phi) + (q + Z(phinumb(i,0), phinumb(i,1))) * (-log(U[phinumb(i,0)] * w[w_indices[phinumb(i,1)]] * phi + 1));
    }
    return(result);
  }
  ld_phi(int numb_, double q_, std::vector<NumericMatrix> philocs_, NumericVector U_, NumericMatrix Z_, NumericVector w_, NumericVector w_indices_){numb = numb_; q = q_; philocs = philocs_; U = U_; Z = Z_; w = w_; w_indices = w_indices_; }
};

//slice sampler for phi when epsilon is not included
class ld_phi_noepsilon: public logdensity {
public:
  int numb;
  std::vector<NumericMatrix> philocs;
  NumericMatrix Z;
  NumericVector w;
  NumericVector w_indices;
  NumericVector naccidentals;
  NumericVector denoms;
  double oldphi;

  double val(double phi){
    NumericMatrix phinumb = philocs[numb];
    double result = 0.0;
    NumericVector denomsnew = clone(denoms);
    for(int i=0; i < phinumb.nrow(); ++i){
      result = result + Z(phinumb(i, 0), phinumb(i, 1)) * log(phi);
      denomsnew(phinumb(i, 0)) = denomsnew(phinumb(i, 0)) + (phi - oldphi) * w[w_indices[phinumb(i, 1)]];
    }
    for(int d = 0; d < Z.nrow(); ++d){
      result = result - naccidentals[d] * log(denomsnew(d));
    }
    return(result);
  }
  ld_phi_noepsilon(int numb_, std::vector<NumericMatrix> philocs_, NumericMatrix Z_, NumericVector w_, NumericVector w_indices_, NumericVector naccidentals_, NumericVector denoms_,  double oldphi_){numb = numb_; philocs = philocs_; Z = Z_; w = w_; w_indices = w_indices_; naccidentals = naccidentals_; denoms = denoms_; oldphi = oldphi_; }
};

// slice sampler for q
class ld_q: public logdensity {
public:
  double sumlogdenom;
  double qa;
  double qb;
  int highestcount;
  NumericVector thecounts;

  double val(double qnew) {
    // from prior
    double result = - qb * qnew + (qa - 1) * log(qnew);

    //from the denominator
    result = result - qnew * sumlogdenom;

    // the gammas
    double lgammaminus = 0;

    // this is a trick for evaluating the gammas
    for(int i = 0; i < (highestcount + 1); ++i){
      result = result + lgammaminus * thecounts[i];
      lgammaminus = lgammaminus + log(i + qnew);
    }
    return(result);
  }
  ld_q(double sumlogdenom_, double qa_, double qb_, int highestcount_, NumericVector thecounts_){
    sumlogdenom = sumlogdenom_; qa = qa_; qb = qb_; highestcount = highestcount_; thecounts = thecounts_;}
};

//slice sampler for p (kernel variables)
class ld_p: public logdensity {
public:
  NumericVector deviations;
  int i;
  double sumexpp;
  int K;
  double sigma;
  NumericVector p;
  int n;
  NumericVector kernely_u;
  double val(double pnew) {
    double diff =  exp(pnew) - exp(p[i]);
    double sumexpp_new = sumexpp + diff;
    double result = 0;
    for(int j = 0; j <= i; ++j){
      result += deviations[j] * log(kernely_u[j] + diff/(1 + 2 * i));
    }
    //normalize
    result -= n * log(sumexpp_new);

    //prior for px
    result -= (pnew * pnew)/(2.0 * sigma * sigma);
    return(result);
  }
  ld_p(NumericVector deviations_, int i_, double sumexpp_, int K_, double sigma_, NumericVector p_, int n_, NumericVector kernely_u_){
    deviations = deviations_; i = i_; sumexpp = sumexpp_; K = K_; sigma = sigma_, p = p_; n = n_; kernely_u = kernely_u_;}
};

double slice(double x0,              // the initial point; use the current value of the mcmc chain
             logdensity* g,          // The log of the possibly unnormalized density
             double w=1.,            // The stepping out window (see Neal 2003)
             double m=INFINITY,      // The max number of times to step out, inf is only safe if bounded
             double lower=-INFINITY, // lower bound of the density
             double upper=INFINITY   // upper bound of the density
);

arma::vec ellipticalslice(arma::vec f,              // the initial point; use the current value of the mcmc chain
                          arma::mat Sigma,    // the precision matrix
                          int n,    // the dimension of f
                          logdensityellipse* logL
);

#endif
