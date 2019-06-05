#ifndef common_h
#define common_h

//#define NDEBUG

//#define BOOST_DISABLE_ASSERTS
//#include <boost/assert.hpp>

#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

inline double logsumexp(const double &a, const double &b){
  return a < b ? b + log(1.0 + exp(a - b)) : a + log(1.0 + exp(b - a));
}

template <class T>
double logsumexp(T &x) {
  double m = *std::max_element(x.begin(), x.end());
  double s = 0.0;
  typename T::iterator it;
  for(it=x.begin(); it!=x.end(); ++it) {
    s += exp(*it-m);
  }
  return(m+log(s));
}

#endif
