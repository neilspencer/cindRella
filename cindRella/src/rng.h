#ifndef rng_h
#define rng_h

#include "common.h"

double rgamma(const double &a, const double &b);
double rbeta(const double &a, const double &b);
double runif();
double rchisq(const double &df);
double znorm();
double rnorm(const double &mu, const double &sd);
double rlogstick(double &alpha);
double rlogbeta(const double &alpha, const double &beta);

double rloggamma(const double &alpha);
mat rmvnorm_prec(int n, int d, mat &Phi);
int rdisc(NumericVector weights);

template<class T>
int rdisc(T &weights, bool normalize=true)
{
    int n = weights.size();
    double s = 1.0;
    double u = Rf_runif(0.0, 1.0);
    if(normalize) {
      s = std::accumulate(weights.begin(), weights.end(), 0.0);
      u = s * Rf_runif(0.0, 1.0);
    }
    double cs = weights[0];
    int i = 0;
    while ((u > cs) & (i < n)) {
      i++;
      cs += weights[i];
    }
    return i;
  }

template<class T>
  int rdisc_log(T &logweights, int n=-1) {
    if(n==-1) n = logweights.size();
    typename T::iterator itb = logweights.begin();
    typename T::iterator ite = logweights.begin() + n;

    double m = *std::max_element(itb, ite);
    double s = 0;
    std::vector<double> weights(n, 0.0);
    for(int i=0; i<n; ++i) {
      weights[i] = exp(logweights[i] - m);
      s += weights[i];
    }
    double u = s*runif();
    double cs = weights[0];
    int i=0;
    while((u>cs) & (i<n)) {
      ++i;
      cs += weights[i];
    }
    return i;
  }


template<class T>
  int rdisc_log_inplace(T &logweights, int n=-1, double u=-1.0) {
    if(n==-1) n = logweights.size();
    if(u==-1) u = runif();
    typename T::iterator itb = logweights.begin();
    typename T::iterator ite = logweights.begin() + n;

    double m = *std::max_element(itb, ite);
    double s = 0;
    //vector<double> weights(n, 0.0);
    for(int i=0; i<n; ++i) {
      logweights[i] = exp(logweights[i] - m);
      s += logweights[i];
    }
    u = s*u;
    double cs = logweights[0];
    int i=0;
    while((u>cs) & (i<n)) {
      ++i;
      cs += logweights[i];
    }
    return i;
  }

#endif

