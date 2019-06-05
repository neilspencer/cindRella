#include "rng.h"
#include "common.h"
#include <Rcpp.h>
using namespace Rcpp;


// the integrand in the expectation (the weight times the likelihood)
double integrand(int ngrid, double U, double qi, double sumwqphi, NumericVector wi, NumericVector contact, NumericVector phi, NumericVector denom, NumericVector Z, int naccidentals, NumericVector w_indices){

  double logresult = - ngrid * lgamma(qi);
  // include info about U importance weight
  logresult = logresult - naccidentals * log(sumwqphi) + U * sumwqphi;

  for(int n = 0; n < naccidentals; ++n){
    logresult += denom[n];
  }
  for(int a = 0; a < ngrid; ++a){
    logresult = logresult  + lgamma(Z[a] + qi) - (Z[a] + qi) * log(U * wi[w_indices[a]] * phi[contact[a]] + 1);
  }
  // adjust for it being uniform over 1/natoms chunk
  logresult += (naccidentals * log(ngrid));
  return(logresult);
}

// the integrand in the expectation (the weight times the likelihood)
double integrand_noPhi(int ngrid, double U, double qi, double sumwq, NumericVector wi, NumericVector contact, NumericVector denom, NumericVector Z, int naccidentals, NumericVector w_indices){

  double logresult = - ngrid * lgamma(qi);
  // include info about U importance weight
  logresult = logresult - naccidentals * log(sumwq) + U * sumwq;

  for(int n = 0; n < naccidentals; ++n){
    logresult += denom[n];
  }
  for(int a = 0; a < ngrid; ++a){
    logresult = logresult  + lgamma(Z[a] + qi) - (Z[a] + qi) * log(U * wi[w_indices[a]] + 1);
  }
  // adjust for it being uniform over 1/natoms chunk
  logresult += (naccidentals * log(ngrid));
  return(logresult);
}


// [[Rcpp::export]]
NumericVector heldoutprob_importance(NumericVector contact, int ncategories_contacts, // data
                                                 NumericMatrix w,  NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, // variables to be sampled
                                                 NumericVector q,
                                                 NumericVector w_indices, // the w index to which each pixel is assigned
                                                 NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                                 NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                                                 NumericVector nzoptions, // the number of possible grid for each accidental assignment
                                                 NumericMatrix optiondiffsx, // the differences in x for each of the options
                                                 NumericMatrix optiondiffsy,// the differences in y for each of the options
                                                 int nsim_per = 10, int status=1) {
  int chain_length = q.size();
  NumericMatrix liks(chain_length, nsim_per);
  int naccidentals = zhome.size();
  int ngrid = contact.size();


  NumericVector phi(ncategories_contacts);
  NumericVector wi;
  double qi;
  NumericVector kernelxi;
  NumericVector kernelyi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    wi = w(_, i);
    qi = q[i];
    NumericVector kernelxi = kernelx(_, i);
    NumericVector kernelyi = kernely(_, i);

    NumericVector logwphi(ngrid);
    double sumwqphi = 0;
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = log(wi[w_indices[a]]) + log(phi(contact[a]));
      sumwqphi += wi[w_indices[a]] * phi(contact[a]);
    }
    sumwqphi = sumwqphi *qi;

    std::vector<NumericVector> propprobs;
    std::vector<NumericVector> propoptions;
    NumericVector denom(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){

      NumericVector propoption(nzoptions[zhome[bb]]); //options for the index of the proposal
      NumericVector theprobs(nzoptions[zhome[bb]]); // the corresponding importance probability
      for(int a = 0; a < nzoptions[zhome[bb]]; ++a){
        theprobs[a] = log(kernelxi[optiondiffsx(zhome[bb], a)]) + log(kernelyi[optiondiffsy(zhome[bb], a)]) + logwphi[zoptions(zhome[bb], a)];
        propoption[a] = zoptions(zhome[bb], a);
      }
      propprobs.push_back(theprobs);

      propoptions.push_back(propoption);
      denom[bb] = logsumexp(theprobs);
    }

    for(int j = 0; j < nsim_per; ++j){

      // sampler for U
      double U = R::rgamma(naccidentals, 1.0/sumwqphi);

      NumericVector zs(naccidentals);
      NumericVector Zs(ngrid);

      // sample mixture membership in x and y direction
      for(int n = 0; n < naccidentals; ++n){
        zs[n] = propoptions[n](rdisc_log(propprobs[n]));
        Zs[zs[n]] = Zs[zs[n]] + 1;
      }

      liks(i, j) = integrand(ngrid, U, qi, sumwqphi, wi, contact, phi, denom, Zs, naccidentals, w_indices);
    }
  }

  return(liks);
}

// fixed value for w
// [[Rcpp::export]]
NumericVector heldoutprob_importance_now(NumericVector contact, int ncategories_contacts, // data
                                     NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, // variables to be sampled
                                     NumericVector q,
                                     NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                     NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                                     NumericVector nzoptions, // the number of possible grid for each accidental assignment
                                     NumericMatrix optiondiffsx, // the differences in x for each of the options
                                     NumericMatrix optiondiffsy,// the differences in y for each of the options
                                     int nsim_per = 10, int status=1) {
  int chain_length = q.size();
  NumericMatrix liks(chain_length, nsim_per);
  int naccidentals = zhome.size();
  int ngrid = contact.size();


  NumericVector phi(ncategories_contacts);
  double qi;
  NumericVector kernelxi;
  NumericVector kernelyi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    qi = q[i];
    NumericVector kernelxi = kernelx(_, i);
    NumericVector kernelyi = kernely(_, i);

    NumericVector logwphi(ngrid);
    double sumqphi = 0;
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] =  log(phi(contact[a]));
      sumqphi +=  phi(contact[a]);
    }
    sumqphi = sumqphi *qi;

    std::vector<NumericVector> propprobs;
    std::vector<NumericVector> propoptions;
    NumericVector denom(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){

      NumericVector propoption(nzoptions[zhome[bb]]); //options for the index of the proposal
      NumericVector theprobs(nzoptions[zhome[bb]]); // the corresponding importance probability
      for(int a = 0; a < nzoptions[zhome[bb]]; ++a){
        theprobs[a] = log(kernelxi[optiondiffsx(zhome[bb], a)]) + log(kernelyi[optiondiffsy(zhome[bb], a)]) + logwphi[zoptions(zhome[bb], a)];
        propoption[a] = zoptions(zhome[bb], a);
      }
      propprobs.push_back(theprobs);

      propoptions.push_back(propoption);
      denom[bb] = logsumexp(theprobs);
    }

    for(int j = 0; j < nsim_per; ++j){

      // sampler for U
      double U = R::rgamma(naccidentals, 1.0/sumqphi);

      NumericVector zs(naccidentals);
      NumericVector Zs(ngrid);

      // sample mixture membership in x and y direction
      for(int n = 0; n < naccidentals; ++n){
        zs[n] = propoptions[n](rdisc_log(propprobs[n]));
        Zs[zs[n]] = Zs[zs[n]] + 1;
      }

      double logresult = - ngrid * lgamma(qi);
      // include info about U importance weight
      logresult = logresult - naccidentals * log(sumqphi) + U * sumqphi;

      for(int n = 0; n < naccidentals; ++n){
        logresult += denom[n];
      }
      for(int a = 0; a < ngrid; ++a){
        logresult = logresult  + lgamma(Zs[a] + qi) - (Zs[a] + qi) * log(U * phi[contact[a]] + 1);
      }
      // adjust for it being uniform over 1/natoms chunk
      logresult += (naccidentals * log(ngrid));

      liks(i, j) = logresult;
    }
  }

  return(liks);
}

// [[Rcpp::export]]
NumericVector heldoutprob_importance_noZ(NumericVector contact, int ncategories_contacts, // data
                                     NumericMatrix w,  NumericMatrix Phi, // variables to be sampled
                                     NumericVector q,
                                     NumericVector w_indices, // the w index to which each pixel is assigned
                                     NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                     int nsim_per = 10, int status=1) {
  int chain_length = q.size();
  NumericMatrix liks(chain_length, nsim_per);
  int naccidentals = zhome.size();
  int ngrid = contact.size();


  NumericVector phi(ncategories_contacts);
  NumericVector wi;
  double qi;
  NumericVector kernelxi;
  NumericVector kernelyi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    wi = w(_, i);
    qi = q[i];

    NumericVector logwphi(ngrid);
    double sumwqphi = 0;
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = log(wi[w_indices[a]]) + log(phi(contact[a]));
      sumwqphi += wi[w_indices[a]] * phi(contact[a]);
    }
    sumwqphi = sumwqphi *qi;

    NumericVector denom(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){
      denom[bb] = logwphi[zhome[bb]];
    }

    NumericVector Zs(ngrid);

    for(int n = 0; n < naccidentals; ++n){
      Zs[zhome[n]] = Zs[zhome[n]] + 1;
    }

    for(int j = 0; j < nsim_per; ++j){

      // sampler for U
      double U = R::rgamma(naccidentals, 1.0/sumwqphi);
      liks(i, j) = integrand(ngrid, U, qi, sumwqphi, wi, contact, phi, denom, Zs, naccidentals, w_indices);
    }
    //Rcout << "percentdone" << (i * 1.0)/(1.0 * chain_length) << std::endl;
  }

  return(liks);
}

// [[Rcpp::export]]
NumericVector heldoutprob_importance_noepsilon(NumericVector contact, int ncategories_contacts, // data
                                     NumericMatrix w,  NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, // variables to be sampled
                                     NumericVector w_indices, // the w index to which each pixel is assigned
                                     NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                     NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                                     NumericVector nzoptions, // the number of possible grid for each accidental assignment
                                     NumericMatrix optiondiffsx, // the differences in x for each of the options
                                     NumericMatrix optiondiffsy,// the differences in y for each of the options
                                     int status=1) {
  int chain_length = w.ncol();
  NumericVector liks(chain_length);
  int naccidentals = zhome.size();
  int ngrid = contact.size();


  NumericVector phi(ncategories_contacts);
  NumericVector wi;
  double qi;
  NumericVector kernelxi;
  NumericVector kernelyi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    wi = w(_, i);
    NumericVector kernelxi = kernelx(_, i);
    NumericVector kernelyi = kernely(_, i);

    NumericVector logwphi(ngrid);
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = log(wi[w_indices[a]]) + log(phi(contact[a]));
    }
    double sssum = logsumexp(logwphi);
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = logwphi[a] - sssum;
    }

    NumericVector denom(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){

      NumericVector propoption(nzoptions[zhome[bb]]); //options for the index of the proposal
      NumericVector theprobs(nzoptions[zhome[bb]]); // the corresponding importance probability
      for(int a = 0; a < nzoptions[zhome[bb]]; ++a){
        theprobs[a] = log(kernelxi[optiondiffsx(zhome[bb], a)]) + log(kernelyi[optiondiffsy(zhome[bb], a)]) + logwphi[zoptions(zhome[bb], a)];
      }
      denom[bb] = logsumexp(theprobs);
    }
    liks(i) = sum(denom)+ (naccidentals * log(ngrid));
    }

  return(liks);
}

// [[Rcpp::export]]
NumericVector heldoutprob_importance_noepsilonnoZ(NumericVector contact, int ncategories_contacts, // data
                                         NumericMatrix w,  NumericMatrix Phi, // variables to be sampled
                                         NumericVector w_indices, // the w index to which each pixel is assigned
                                         NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                         int status=1) {
  int chain_length = w.ncol();
  NumericVector liks(chain_length);
  int naccidentals = zhome.size();
  int ngrid = contact.size();


  NumericVector phi(ncategories_contacts);
  NumericVector wi;
  double qi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    phi = Phi(_, i);
    wi = w(_, i);

    NumericVector logwphi(ngrid);
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = log(wi[w_indices[a]]) + log(phi(contact[a]));
    }
    double sssum = logsumexp(logwphi);
    for(int a = 0; a < ngrid; ++a){
      logwphi[a] = logwphi[a] - sssum;
    }

    NumericVector prob(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){
      prob[bb] = logwphi[zhome[bb]];
    }
    liks(i) = sum(prob)+ (naccidentals * log(ngrid));
    }

  return(liks);
}

// [[Rcpp::export]]
NumericVector heldoutprob_importance_noPhi(NumericVector contact, // data
                                     NumericMatrix w, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, // variables to be sampled
                                     NumericVector q,
                                     NumericVector w_indices, // the w index to which each pixel is assigned
                                     NumericVector zhome, // the pixel index of the locations of each of the observed accidentals
                                     NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                                     NumericVector nzoptions, // the number of possible grid for each accidental assignment
                                     NumericMatrix optiondiffsx, // the differences in x for each of the options
                                     NumericMatrix optiondiffsy,// the differences in y for each of the options
                                     int nsim_per = 10, int status=1) {
  int chain_length = q.size();
  NumericMatrix liks(chain_length, nsim_per);
  int naccidentals = zhome.size();
  int ngrid = contact.size();

  NumericVector wi;
  double qi;
  NumericVector kernelxi;
  NumericVector kernelyi;

  for(int i = 0; i < chain_length; ++i){
    // access the parameters in the ith iteration of the chain
    wi = w(_, i);
    qi = q[i];
    NumericVector kernelxi = kernelx(_, i);
    NumericVector kernelyi = kernely(_, i);

    NumericVector logw(ngrid);
    double sumwq = 0;
    for(int a = 0; a < ngrid; ++a){
      logw[a] = log(wi[w_indices[a]]);
      sumwq += wi[w_indices[a]];
    }
    sumwq = sumwq *qi;

    std::vector<NumericVector> propprobs;
    std::vector<NumericVector> propoptions;
    NumericVector denom(naccidentals);


    for(int bb = 0; bb < naccidentals; ++bb){

      NumericVector propoption(nzoptions[zhome[bb]]); //options for the index of the proposal
      NumericVector theprobs(nzoptions[zhome[bb]]); // the corresponding importance probability
      for(int a = 0; a < nzoptions[zhome[bb]]; ++a){
        theprobs[a] = log(kernelxi[optiondiffsx(zhome[bb], a)]) + log(kernelyi[optiondiffsy(zhome[bb], a)]) + logw[zoptions(zhome[bb], a)];
        propoption[a] = zoptions(zhome[bb], a);
      }
      propprobs.push_back(theprobs);

      propoptions.push_back(propoption);
      denom[bb] = logsumexp(theprobs);
    }

    for(int j = 0; j < nsim_per; ++j){

      // sampler for U
      double U = R::rgamma(naccidentals, 1.0/sumwq);

      NumericVector zs(naccidentals);
      NumericVector Zs(ngrid);

      // sample mixture membership in x and y direction
      for(int n = 0; n < naccidentals; ++n){
        zs[n] = propoptions[n](rdisc_log(propprobs[n]));
        Zs[zs[n]] = Zs[zs[n]] + 1;
      }

      liks(i, j) = integrand_noPhi(ngrid, U, qi, sumwq, wi, contact, denom, Zs, naccidentals, w_indices);
    }
  }

  return(liks);
}

