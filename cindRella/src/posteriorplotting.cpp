#include "rng.h"
#include "common.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List plot_posterior(NumericVector contact, // data
                    NumericMatrix w,  NumericMatrix Phi, NumericMatrix kernelx, NumericMatrix kernely, // variables to be sampled
                    NumericVector q,
                    NumericVector w_indices, // the w index to which each pixel is assigned
                    NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                    NumericVector nzoptions, // the number of possible grid for each accidental assignment
                    NumericMatrix optiondiffsx, // the differences in x for each of the options
                    NumericMatrix optiondiffsy// the differences in y for each of the options
) {
  int chain_length = q.size();
  int ngrid = w_indices.size();
  int Kx = kernelx.nrow() - 1;
  int Ky = kernely.nrow() - 1;
  NumericMatrix predprobs(ngrid, chain_length);
  NumericVector logweights(chain_length);

  NumericVector phi(Phi.nrow());
  NumericVector wi(w.nrow());
  double qi;
  NumericVector kernelxi(Kx + 1);
  NumericVector kernelyi(Ky + 1);

  int ind = 0;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    wi = w(_, i);
    qi = q[i];
    kernelxi = kernelx(_, i);
    kernelyi = kernely(_, i);

    double sumphiw = 0;
    for(int a = 0; a < ngrid; ++a){
      sumphiw += phi[contact[a]] * wi[w_indices[a]];
    }

    // importance sample U
    double U = R::rgamma(1.0, 1.0/(sumphiw * qi));
    NumericVector zdenoms(ngrid);
    NumericVector zprobs(ngrid);
    // weight includes terms that are constant across locations but depend on u,w,q,phi
    double logweight  = U * qi * sumphiw - log(sumphiw);
    double toaddtologweight = 0;
    for(int a = 0; a < ngrid; ++a){
      zdenoms[a] = log(1.0 + U * wi[w_indices[a]] * phi[contact[a]]);
      zprobs[a] = wi[w_indices[a]] * phi[contact[a]] * exp(-zdenoms[a]);
      toaddtologweight += zdenoms[a];
    }
    logweights[i] = logweight - qi * toaddtologweight;
    for(int a = 0; a < ngrid; ++a){
      for(int rr = 0; rr < nzoptions[a]; ++rr){
        predprobs(zoptions(a, rr), i) += zprobs[a] * kernelxi[optiondiffsx(a, rr)] * kernelyi[optiondiffsy(a, rr)];
      }
    }
  }

  //aggregate results with weights
  NumericVector results(ngrid);
  for(int a = 0; a < ngrid; ++a){
    NumericVector coord(chain_length);
    for(int j = 0; j < chain_length; ++j){
      coord(j) = log(predprobs(a, j)) + logweights[j];
    }
    results[a] = logsumexp(coord) - log(chain_length);
  }
  //return results;
  return List::create(_["predprobs"]= predprobs, _["logweights"]= logweights, _["results"] = results);
}
