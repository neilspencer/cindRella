#include "rng.h"
#include "slice.h"// contains all of our slice samplers
#include "common.h"
#include "probabilitytracking.h"

// [[Rcpp::export]]
List mcmc_inference(NumericVector shoenums, // the shoe number of each accidental
                        NumericMatrix contacts, // the value of the contact surface at each gridpoint (numbers of shoes by number of atoms)
                        int ncategories_contacts, // the number of categories of contact surface value
                        int nshoes, // the number of shoes in the training set
                        NumericVector z_init, // initialization of the assignments (from 0 to Kx * ky + 1)
                        NumericVector w_init, // the initializations of values of the w parameters
                        NumericVector w_indices, // the w index to which each gridpoint is assigned
                        NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                        NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                        NumericVector nzoptions, // the number of possible grid for each accidental assignment
                        NumericMatrix optiondiffsx, // the differences in x for each of the options
                        NumericMatrix optiondiffsy,// the differences in y for each of the options
                        NumericVector phi_init, // the initialization of the phi variables
                        arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                        int Kx, int Ky, // the kernel width in the x and y direction
                        NumericVector px_init, NumericVector py_init, // the initial unnormalized log kernel values in the x and y direction
                        double sigma, // the hyperparameter for the lognormal for the kernel values
                        NumericVector u_init, // the initialization of the u auxiliary variable
                        double q_init, // the initialization for q hyperparameter
                        double qa, double qb, // the alpha and beta parameters for the prior on q
                        double rhox_init, // the smoothing values for the w's in the x direction
                        double rhoy_init, // the smoothing values for the w's in the y direction
                        int nsim = 10, // the number of simulations to be run
                        int print_output = 1, // if we want additional print statements while running
                        int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                        int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w = w_init;
  int nws = w.size(); // the number of unique w variables
  int ngrid = contacts.ncol();
  arma::vec logw(nws);
  for(int a = 0; a < nws; ++a){
    logw(a) = log(w[a]);
  }
  NumericVector Phi = phi_init;
  NumericVector px = px_init;
  NumericVector py = py_init;
  NumericVector U = u_init;
  double q = q_init;
  double rhox = rhox_init;
  double rhoy = rhoy_init;
  int n = z_init.size();   // number of accidentals

  //create philocs matrix. This separates the grid into ncategories different lists,
  // each containing references to which grid are in each category. It also says which shoe
  std::vector<NumericMatrix> philocs;
  for(int j = 0; j < ncategories_contacts; ++j){
    std::vector<int> indices;
    std::vector<int> shoenums;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        if(contacts(d, a) == j){
          indices.push_back(a);
          shoenums.push_back(d);
        }
      }
    }

    NumericMatrix phij(shoenums.size(), 2);
    for(int i = 0; i < shoenums.size(); ++i){
      phij(i,0) = shoenums[i];
      phij(i,1) = indices[i];
    }
    philocs.push_back(phij);
  }


  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each gridpoint
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = z_init[i];
    Z(shoenums[i], zoptions(zhome[i], z[i]))++;
  }

  //find grid of each point,
  // and keep track of how much they deviate from their "home" gridpoint in each direction
  NumericVector deviationsx(Kx + 1); // the counts of how many of each deviation
  NumericVector deviationsy(Ky + 1);

  for(int i = 0; i < n; ++i){
    deviationsx[optiondiffsx(zhome[i], z[i])] += 1;
    deviationsy[optiondiffsy(zhome[i], z[i])] += 1;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }

  // it is convenient to also maintain some parameters on the log scale as well
  NumericVector logU(nshoes);
  for(int i=0; i < nshoes; ++i){
    logU[i] = log(U[i]);
  }

  NumericVector logPhi(ncategories_contacts);
  for(int i= 0; i < ncategories_contacts; ++i){
    logPhi[i] = log(Phi[i]); // initialize at reasonable intervals
  }

  // the normalization for the kernels
  double sumexppx = 0;
  double sumexppy = 0;
  for(int i = 0; i < (Kx + 1); ++i){
    sumexppx += exp(px[i]);
  }
  for(int i = 0; i < (Ky + 1); ++i){
    sumexppy += exp(py[i]);
  }

  // The kernels on the non-normalized scale
  // unnormalized
  NumericVector kernelx_u(Kx + 1);
  NumericVector kernely_u(Ky + 1);
  // normalized
  NumericVector kernelx(Kx + 1);
  NumericVector kernely(Ky + 1);

  // fill them
  kernelx_u[Kx] = exp(px[Kx])/(2 * Kx + 1);
  kernely_u[Ky] = exp(py[Ky])/(2 * Ky + 1);
  for(int k = 1; k < (Kx + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernelx_u(Kx - k) = kernelx_u(Kx + 1 - k) + exp(px[Kx - k])/(2 * Kx + 1 - 2 * k);
  }

  for(int k = 1; k < (Ky + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernely_u(Ky - k) = kernely_u(Ky + 1 - k) + exp(py[Ky - k])/(2 * Ky + 1 - 2 * k);
  }

  // normalize
  for(int k = 0; k < (Kx + 1); ++k){
    kernelx(k) = kernelx_u(k)/sumexppx;
  }

  for(int k = 0; k < (Ky + 1); ++k){
    kernely(k) = kernely_u(k)/sumexppy;
  }

  arma::mat precw = arma::eye(nws, nws) + rhox * Ax + rhoy * Ay;

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_w;
  std::vector<NumericVector> post_U;
  std::vector<NumericVector> post_z;
  std::vector<NumericVector> post_phi;
  std::vector<double> post_q;
  std::vector<double> post_rhox;
  std::vector<double> post_rhoy;
  std::vector<NumericVector> post_kernelx;
  std::vector<NumericVector> post_kernely;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update u with slice # updated 2018 - note that step out width might need to change
    for(int d = 0; d < nshoes; ++d){

      ld_u dens_u(Phi, q, naccidentals, d, w, Z, contacts, ngrid, w_indices);
      // gave an upper bound because I was running into numerical problems sometimes early in the chain
      U[d] = slice(U[d], &dens_u, 20 * (sqrt(naccidentals[d]))/(ngrid * q), 100, 0, 100000);
      logU[d] = log(U[d]);
    }

    //update ws with elliptical slice
    ld_logw dens_logw(Phi, q, U, Z, contacts, nshoes, ngrid, w_indices);

    logw = ellipticalslice(logw, precw, nws, &dens_logw);

    for(int a = 0; a < nws; ++a){
      w[a] = exp(logw[a]);
    }

    if(print_output == 1){
      for(int a = 0; a < nws; ++a){
        Rcout << "w" << a << " "<< w[a] << std::endl;
      }
    }

    //update Phis with slice (with a uniform prior)
    for(int i = 0; i < ncategories_contacts; ++i){
      ld_phi dens_phi(i, q, philocs, U, Z, w, w_indices);
      Phi[i] = slice(Phi[i], &dens_phi, 0.01, 100, 0, 1);
      logPhi[i] = log(Phi[i]);
      if(print_output == 1){
        Rcout << "phi" << i << " "<< Phi[i] << std::endl;
      }
    }

    //update q with slice (with gamma(qa, qb) prior)
    int highestcount = max(Z);
    NumericVector counts(highestcount + 1);
    for(int ww = 0; ww < (highestcount + 1); ++ww){
      counts[ww] = 0;
    }

    double logsumdenom = 0.0;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        counts[Z(d, a)]++;
        logsumdenom = logsumdenom + log(U[d] * w[w_indices[a]] * Phi[contacts(d,a)] + 1);
      }
    }

    ld_q dens_q(logsumdenom, qa, qb, highestcount, counts);
    q = slice(q, &dens_q, 0.2, 100, 0, 500);
    if(print_output == 1){
      Rcout << "q " << q << std::endl;
    }

    // update kernelx with slice
    for(int i = 0; i < (Kx + 1); ++i){
      double old = px[i];
      ld_p dens_p(deviationsx, i, sumexppx, Kx, sigma, px, n, kernelx_u);
      px[i] = slice(px[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(px[i] - old) - 1.0);
      sumexppx = sumexppx + diff;
      for(int j = 0; j <= i; ++j){
        kernelx_u[j] =  kernelx_u[j] + diff/(1 + 2 * i);
      }
    }

    // update kernely with slice
    for(int i = 0; i < (Ky + 1); ++i){
      double old = py[i];
      ld_p dens_p(deviationsy, i, sumexppy, Ky, sigma, py, n, kernely_u);
      py[i] = slice(py[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(py[i] - old) - 1.0);
      sumexppy = sumexppy + diff;
      for(int j = 0; j <= i; ++j){
        kernely_u[j] =  kernely_u[j] + diff/(1 + 2 * i);
      }
    }

    // normalize the kernels
    for(int k = 0; k < (Kx + 1); ++k){
      kernelx(k) = kernelx_u(k)/sumexppx;
    }

    for(int k = 0; k < (Ky + 1); ++k){
      kernely(k) = kernely_u(k)/sumexppy;
    }
    if(print_output == 1){
      for(int i = 0; i < (Kx + 1); ++i){
        Rcout << "kernelx" << i << " "<< kernelx[i] << std::endl;
      }
      for(int i = 0; i < (Ky + 1); ++i){
        Rcout << "kernely" << i << " "<< kernely[i] << std::endl;
      }
    }


    //Update the gridpoint assignments via Gibbs
    int accind = 0;
    // update allocations shoe by shoe
    for(int d = 0; d < nshoes; ++d){
      // update by accidental
      for(int i = 0; i < naccidentals[d]; ++i){
        // remove accidental from current location
        Z(d, zoptions(zhome[accind], z[accind])) -= 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] -= 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] -= 1;
        NumericVector logwt(nzoptions[zhome[accind]]);
        for(int ind = 0; ind < nzoptions[zhome[accind]]; ++ind){
          logwt[ind] = logw[w_indices[zoptions(zhome[accind], ind)]] + logPhi(contacts(d, zoptions(zhome[accind], ind))) + log(q + Z(d, zoptions(zhome[accind], ind))) - log(U[d] * w[w_indices[zoptions(zhome[accind], ind)]] * Phi[contacts(d, zoptions(zhome[accind], ind))] + 1.0);
          logwt[ind] += log(kernelx[optiondiffsx(zhome[accind], ind)]) + log(kernely[optiondiffsy(zhome[accind], ind)]);
        }

        // sample new index according to normalized conditional probabilities
        z[accind] = rdisc_log(logwt);
        Z(d, zoptions(zhome[accind], z[accind])) += 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] += 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] += 1;
        accind++;
      }
    }


    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_w.push_back(clone(w));
      post_U.push_back(clone(U));
      post_z.push_back(clone(z));
      post_phi.push_back(clone(Phi));
      post_q.push_back(q);
      post_kernelx.push_back(clone(kernelx));
      post_kernely.push_back(clone(kernely));
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob(deviationsx, deviationsy, px, py, kernelx, kernely, sigma, Kx, Ky, ngrid, w, w_indices, logw, U, q, Z, contacts, naccidentals, n, nshoes, logPhi, Phi, qa, qb, precw);
        logliks.push_back(loglike);
      }
      Rcout << "Our model mcmc iteration " << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["kernelx"] = post_kernelx, _["kernely"] = post_kernely,
    _["w"]=post_w,  _["U"]=post_U, _["q"]=post_q,  _["phi"] = post_phi);
}

// [[Rcpp::export]]
List mcmc_inference_noZ(NumericVector shoenums, // the shoe number of each accidental
                            NumericMatrix contacts, // the value of the contact surface at each gridpoint (numbers of shoes by number of atoms)
                            int ncategories_contacts, // the number of categories of contact surface value
                            int nshoes, // the number of shoes in the training set
                            NumericVector w_init, // the initializations of values of the w parameters
                            NumericVector w_indices, // the w index to which each gridpoint is assigned
                            NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                            NumericVector phi_init, // the initialization of the phi variables
                            arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                            NumericVector u_init, // the initialization of the u auxiliary variable
                            double q_init, // the initialization for q hyperparameter
                            double qa, double qb, // the alpha and beta parameters for the prior on q
                            double rhox_init, // the smoothing values for the w's in the x direction
                            double rhoy_init, // the smoothing values for the w's in the y direction
                            int nsim = 10, // the number of simulations to be run
                            int print_output = 1, // if we want additional print statements while running
                            int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                            int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w = w_init;
  int nws = w.size(); // the number of unique w variables
  int ngrid = contacts.ncol();
  arma::vec logw(nws);
  for(int a = 0; a < nws; ++a){
    logw(a) = log(w[a]);
  }
  NumericVector Phi = phi_init;
  NumericVector U = u_init;
  double q = q_init;
  double rhox = rhox_init;
  double rhoy = rhoy_init;
  int n = zhome.size();   // number of accidentals

  //create philocs matrix. This separates the grid into ncategories different lists,
  // each containing references to which grid are in each category. It also says which shoe
  std::vector<NumericMatrix> philocs;
  for(int j = 0; j < ncategories_contacts; ++j){
    std::vector<int> indices;
    std::vector<int> shoenums;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        if(contacts(d, a) == j){
          indices.push_back(a);
          shoenums.push_back(d);
        }
      }
    }

    NumericMatrix phij(shoenums.size(), 2);
    for(int i = 0; i < shoenums.size(); ++i){
      phij(i,0) = shoenums[i];
      phij(i,1) = indices[i];
    }
    philocs.push_back(phij);
  }


  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each atom
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = zhome[i];
    Z(shoenums[i], zhome[i])++;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }

  // it is convenient to also maintain some parameters on the log scale as well
  NumericVector logU(nshoes);
  for(int i=0; i < nshoes; ++i){
    logU[i] = log(U[i]);
  }

  NumericVector logPhi(ncategories_contacts);
  for(int i= 0; i < ncategories_contacts; ++i){
    logPhi[i] = log(Phi[i]); // initialize at reasonable intervals
  }

  arma::mat precw = arma::eye(nws, nws) + rhox * Ax + rhoy * Ay;

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_w;
  std::vector<NumericVector> post_U;
  std::vector<NumericVector> post_phi;
  std::vector<double> post_q;
  std::vector<double> post_rhox;
  std::vector<double> post_rhoy;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update u with slice # updated 2018 - note that step out width might need to change
    for(int d = 0; d < nshoes; ++d){

      ld_u dens_u(Phi, q, naccidentals, d, w, Z, contacts, ngrid, w_indices);
      // gave an upper bound because I was running into numerical problems sometimes early in the chain
      U[d] = slice(U[d], &dens_u, 20 * (sqrt(naccidentals[d]))/(ngrid * q), 100, 0, 100000);
      logU[d] = log(U[d]);
    }

    //update ws with elliptical slice
    ld_logw dens_logw(Phi, q, U, Z, contacts, nshoes, ngrid, w_indices);

    logw = ellipticalslice(logw, precw, nws, &dens_logw);

    for(int a = 0; a < nws; ++a){
      w[a] = exp(logw[a]);
    }

    if(print_output == 1){
      for(int a = 0; a < nws; ++a){
        Rcout << "w" << a << " "<< w[a] << std::endl;
      }
    }

    //update Phis with slice (with a uniform prior)
    for(int i = 0; i < ncategories_contacts; ++i){
      ld_phi dens_phi(i, q, philocs, U, Z, w, w_indices);
      Phi[i] = slice(Phi[i], &dens_phi, 0.01, 100, 0, 1);
      logPhi[i] = log(Phi[i]);
      if(print_output == 1){
        Rcout << "phi" << i << " "<< Phi[i] << std::endl;
      }
    }

    //update q with slice (with gamma(qa, qb) prior)
    int highestcount = max(Z);
    NumericVector counts(highestcount + 1);
    for(int ww = 0; ww < (highestcount + 1); ++ww){
      counts[ww] = 0;
    }

    double logsumdenom = 0.0;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        counts[Z(d, a)]++;
        logsumdenom = logsumdenom + log(U[d] * w[w_indices[a]] * Phi[contacts(d,a)] + 1);
      }
    }

    ld_q dens_q(logsumdenom, qa, qb, highestcount, counts);
    q = slice(q, &dens_q, 0.2, 100, 0, 500);
    if(print_output == 1){
      Rcout << "q " << q << std::endl;
    }


    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_w.push_back(clone(w));
      post_U.push_back(clone(U));
      post_phi.push_back(clone(Phi));
      post_q.push_back(q);
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob_noZ(ngrid, w, w_indices, logw, U, q, Z, contacts, naccidentals, n, nshoes, logPhi, Phi, qa, qb, precw);
        logliks.push_back(loglike);
      }
      Rcout << "Our model no kernel mcmc iteration " << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["w"]=post_w,  _["U"]=post_U, _["q"]=post_q,  _["phi"] = post_phi);
}

// [[Rcpp::export]]
List mcmc_inference_now(NumericVector shoenums, // the shoe number of each accidental
                            NumericMatrix contacts, // the value of the contact surface at each gridpoint (numbers of shoes by number of atoms)
                            int ncategories_contacts, // the number of categories of contact surface value
                            int nshoes, // the number of shoes in the training set
                            NumericVector z_init, // initialization of the assignments (from 0 to Kx * ky + 1)
                            NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                            NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                            NumericVector nzoptions, // the number of possible grid for each accidental assignment
                            NumericMatrix optiondiffsx, // the differences in x for each of the options
                            NumericMatrix optiondiffsy,// the differences in y for each of the options
                            NumericVector phi_init, // the initialization of the phi variables
                            arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                            int Kx, int Ky, // the kernel width in the x and y direction
                            NumericVector px_init, NumericVector py_init, // the initial unnormalized log kernel values in the x and y direction
                            double sigma, // the hyperparameter for the lognormal for the kernel values
                            NumericVector u_init, // the initialization of the u auxiliary variable
                            double q_init, // the initialization for q hyperparameter
                            double qa, double qb, // the alpha and beta parameters for the prior on q
                            int nsim = 10, // the number of simulations to be run
                            int print_output = 1, // if we want additional print statements while running
                            int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                            int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w(1, 1.0);
  int ngrid = contacts.ncol();
  NumericVector w_indices(ngrid);
  NumericVector Phi = phi_init;
  NumericVector px = px_init;
  NumericVector py = py_init;
  NumericVector U = u_init;
  double q = q_init;
  int n = z_init.size();   // number of accidentals

  //create philocs matrix. This separates the grid into ncategories different lists,
  // each containing references to which grid are in each category. It also says which shoe
  std::vector<NumericMatrix> philocs;
  for(int j = 0; j < ncategories_contacts; ++j){
    std::vector<int> indices;
    std::vector<int> shoenums;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        if(contacts(d, a) == j){
          indices.push_back(a);
          shoenums.push_back(d);
        }
      }
    }

    NumericMatrix phij(shoenums.size(), 2);
    for(int i = 0; i < shoenums.size(); ++i){
      phij(i,0) = shoenums[i];
      phij(i,1) = indices[i];
    }
    philocs.push_back(phij);
  }


  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each atom
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = z_init[i];
    Z(shoenums[i], zoptions(zhome[i], z[i]))++;
  }

  //find grid of each point,
  // and keep track of how much they deviate from their "home" gridpoint in each direction
  NumericVector deviationsx(Kx + 1); // the counts of how many of each deviation
  NumericVector deviationsy(Ky + 1);

  for(int i = 0; i < n; ++i){
    deviationsx[optiondiffsx(zhome[i], z[i])] += 1;
    deviationsy[optiondiffsy(zhome[i], z[i])] += 1;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }

  // it is convenient to also maintain some parameters on the log scale as well
  NumericVector logU(nshoes);
  for(int i=0; i < nshoes; ++i){
    logU[i] = log(U[i]);
  }

  NumericVector logPhi(ncategories_contacts);
  for(int i= 0; i < ncategories_contacts; ++i){
    logPhi[i] = log(Phi[i]); // initialize at reasonable intervals
  }

  // the normalization for the kernels
  double sumexppx = 0;
  double sumexppy = 0;
  for(int i = 0; i < (Kx + 1); ++i){
    sumexppx += exp(px[i]);
  }
  for(int i = 0; i < (Ky + 1); ++i){
    sumexppy += exp(py[i]);
  }

  // The kernels on the non-normalized scale
  // unnormalized
  NumericVector kernelx_u(Kx + 1);
  NumericVector kernely_u(Ky + 1);
  // normalized
  NumericVector kernelx(Kx + 1);
  NumericVector kernely(Ky + 1);

  // fill them
  kernelx_u[Kx] = exp(px[Kx])/(2 * Kx + 1);
  kernely_u[Ky] = exp(py[Ky])/(2 * Ky + 1);
  for(int k = 1; k < (Kx + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernelx_u(Kx - k) = kernelx_u(Kx + 1 - k) + exp(px[Kx - k])/(2 * Kx + 1 - 2 * k);
  }

  for(int k = 1; k < (Ky + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernely_u(Ky - k) = kernely_u(Ky + 1 - k) + exp(py[Ky - k])/(2 * Ky + 1 - 2 * k);
  }

  // normalize
  for(int k = 0; k < (Kx + 1); ++k){
    kernelx(k) = kernelx_u(k)/sumexppx;
  }

  for(int k = 0; k < (Ky + 1); ++k){
    kernely(k) = kernely_u(k)/sumexppy;
  }

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_U;
  std::vector<NumericVector> post_z;
  std::vector<NumericVector> post_phi;
  std::vector<double> post_q;
  std::vector<double> post_rhox;
  std::vector<double> post_rhoy;
  std::vector<NumericVector> post_kernelx;
  std::vector<NumericVector> post_kernely;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update u with slice # updated 2018 - note that step out width might need to change
    for(int d = 0; d < nshoes; ++d){

      ld_u dens_u(Phi, q, naccidentals, d, w, Z, contacts, ngrid, w_indices);
      // gave an upper bound because I was running into numerical problems sometimes early in the chain
      U[d] = slice(U[d], &dens_u, 20 * (sqrt(naccidentals[d]))/(ngrid * q), 100, 0, 100000);
      logU[d] = log(U[d]);
    }

    //update Phis with slice (with a uniform prior)
    for(int i = 0; i < ncategories_contacts; ++i){
      ld_phi dens_phi(i, q, philocs, U, Z, w, w_indices);
      Phi[i] = slice(Phi[i], &dens_phi, 0.01, 100, 0, 1);
      logPhi[i] = log(Phi[i]);
      if(print_output == 1){
        Rcout << "phi" << i << " "<< Phi[i] << std::endl;
      }
    }

    //update q with slice (with gamma(qa, qb) prior)
    int highestcount = max(Z);
    NumericVector counts(highestcount + 1);
    for(int ww = 0; ww < (highestcount + 1); ++ww){
      counts[ww] = 0;
    }

    double logsumdenom = 0.0;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        counts[Z(d, a)]++;
        logsumdenom = logsumdenom + log(U[d] * Phi[contacts(d,a)] + 1);
      }
    }

    ld_q dens_q(logsumdenom, qa, qb, highestcount, counts);
    q = slice(q, &dens_q, 0.2, 100, 0, 500);
    if(print_output == 1){
      Rcout << "q " << q << std::endl;
    }

    // update kernelx with slice # updated 2018
    for(int i = 0; i < (Kx + 1); ++i){
      double old = px[i];
      ld_p dens_p(deviationsx, i, sumexppx, Kx, sigma, px, n, kernelx_u);
      px[i] = slice(px[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(px[i] - old) - 1.0);
      sumexppx = sumexppx + diff;
      for(int j = 0; j <= i; ++j){
        kernelx_u[j] =  kernelx_u[j] + diff/(1 + 2 * i);
      }
    }

    // update kernely with slice # updated 2018
    for(int i = 0; i < (Ky + 1); ++i){
      double old = py[i];
      ld_p dens_p(deviationsy, i, sumexppy, Ky, sigma, py, n, kernely_u);
      py[i] = slice(py[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(py[i] - old) - 1.0);
      sumexppy = sumexppy + diff;
      for(int j = 0; j <= i; ++j){
        kernely_u[j] =  kernely_u[j] + diff/(1 + 2 * i);
      }
    }

    // normalize the kernels
    for(int k = 0; k < (Kx + 1); ++k){
      kernelx(k) = kernelx_u(k)/sumexppx;
    }

    for(int k = 0; k < (Ky + 1); ++k){
      kernely(k) = kernely_u(k)/sumexppy;
    }
    if(print_output == 1){
      for(int i = 0; i < (Kx + 1); ++i){
        Rcout << "kernelx" << i << " "<< kernelx[i] << std::endl;
      }
      for(int i = 0; i < (Ky + 1); ++i){
        Rcout << "kernely" << i << " "<< kernely[i] << std::endl;
      }
    }


    //Update the gridpoint assignments via Gibbs # updated 2018
    int accind = 0;
    // update allocations shoe by shoe
    for(int d = 0; d < nshoes; ++d){
      // update by accidental
      for(int i = 0; i < naccidentals[d]; ++i){
        // remove accidental from current location
        Z(d, zoptions(zhome[accind], z[accind])) -= 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] -= 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] -= 1;
        NumericVector logwt(nzoptions[zhome[accind]]);
        for(int ind = 0; ind < nzoptions[zhome[accind]]; ++ind){
          logwt[ind] = logPhi(contacts(d, zoptions(zhome[accind], ind))) + log(q + Z(d, zoptions(zhome[accind], ind))) - log(U[d] * Phi[contacts(d, zoptions(zhome[accind], ind))] + 1.0);
          logwt[ind] += log(kernelx[optiondiffsx(zhome[accind], ind)]) + log(kernely[optiondiffsy(zhome[accind], ind)]);
        }

        // sample new index according to normalized conditional probabilities
        z[accind] = rdisc_log(logwt);
        Z(d, zoptions(zhome[accind], z[accind])) += 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] += 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] += 1;
        accind++;
      }
    }


    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_U.push_back(clone(U));
      post_z.push_back(clone(z));
      post_phi.push_back(clone(Phi));
      post_q.push_back(q);
      post_kernelx.push_back(clone(kernelx));
      post_kernely.push_back(clone(kernely));
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob_now(deviationsx, deviationsy, px, py, kernelx, kernely, sigma, Kx, Ky, ngrid, U, q, Z, contacts, naccidentals, n, nshoes, logPhi, Phi, qa, qb);
        logliks.push_back(loglike);
      }
      Rcout << "Our model uniform w mcmc iteration "  << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["kernelx"] = post_kernelx, _["kernely"] = post_kernely,  _["U"]=post_U, _["q"]=post_q,  _["phi"] = post_phi);

}


// [[Rcpp::export]]
List mcmc_inference_noepsilon(NumericVector shoenums, // the shoe number of each accidental
                            NumericMatrix contacts, // the value of the contact surface at each gridpoint (numbers of shoes by number of atoms)
                            int ncategories_contacts, // the number of categories of contact surface value
                            int nshoes, // the number of shoes in the training set
                            NumericVector z_init, // initialization of the assignments (from 0 to Kx * ky + 1)
                            NumericVector w_init, // the initializations of values of the w parameters
                            NumericVector w_indices, // the w index to which each gridpoint is assigned
                            NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                            NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                            NumericVector nzoptions, // the number of possible grid for each accidental assignment
                            NumericMatrix optiondiffsx, // the differences in x for each of the options
                            NumericMatrix optiondiffsy,// the differences in y for each of the options
                            NumericVector phi_init, // the initialization of the phi variables
                            arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                            int Kx, int Ky, // the kernel width in the x and y direction
                            NumericVector px_init, NumericVector py_init, // the initial unnormalized log kernel values in the x and y direction
                            double sigma, // the hyperparameter for the lognormal for the kernel values
                            double rhox_init, // the smoothing values for the w's in the x direction
                            double rhoy_init, // the smoothing values for the w's in the y direction
                            int nsim = 10, // the number of simulations to be run
                            int print_output = 1, // if we want additional print statements while running
                            int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                            int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w = w_init;
  int nws = w.size(); // the number of unique w variables
  int ngrid = contacts.ncol();
  arma::vec logw(nws);
  for(int a = 0; a < nws; ++a){
    logw(a) = log(w[a]);
  }
  NumericVector Phi = phi_init;
  NumericVector px = px_init;
  NumericVector py = py_init;
  double rhox = rhox_init;
  double rhoy = rhoy_init;
  int n = z_init.size();   // number of accidentals

  //create philocs matrix. This separates the grid into ncategories different lists,
  // each containing references to which grid are in each category. It also says which shoe
  std::vector<NumericMatrix> philocs;
  for(int j = 0; j < ncategories_contacts; ++j){
    std::vector<int> indices;
    std::vector<int> shoenums;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        if(contacts(d, a) == j){
          indices.push_back(a);
          shoenums.push_back(d);
        }
      }
    }

    NumericMatrix phij(shoenums.size(), 2);
    for(int i = 0; i < shoenums.size(); ++i){
      phij(i,0) = shoenums[i];
      phij(i,1) = indices[i];
    }
    philocs.push_back(phij);
  }


  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each atom
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = z_init[i];
    Z(shoenums[i], zoptions(zhome[i], z[i]))++;
  }

  //find grid of each point,
  // and keep track of how much they deviate from their "home" gridpoint in each direction
  NumericVector deviationsx(Kx + 1); // the counts of how many of each deviation
  NumericVector deviationsy(Ky + 1);

  for(int i = 0; i < n; ++i){
    deviationsx[optiondiffsx(zhome[i], z[i])] += 1;
    deviationsy[optiondiffsy(zhome[i], z[i])] += 1;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }


  NumericVector logPhi(ncategories_contacts);
  for(int i= 0; i < ncategories_contacts; ++i){
    logPhi[i] = log(Phi[i]); // initialize at reasonable intervals
  }

  // the normalization for the kernels
  double sumexppx = 0;
  double sumexppy = 0;
  for(int i = 0; i < (Kx + 1); ++i){
    sumexppx += exp(px[i]);
  }
  for(int i = 0; i < (Ky + 1); ++i){
    sumexppy += exp(py[i]);
  }

  // The kernels on the non-normalized scale
  // unnormalized
  NumericVector kernelx_u(Kx + 1);
  NumericVector kernely_u(Ky + 1);
  // normalized
  NumericVector kernelx(Kx + 1);
  NumericVector kernely(Ky + 1);

  // fill them
  kernelx_u[Kx] = exp(px[Kx])/(2 * Kx + 1);
  kernely_u[Ky] = exp(py[Ky])/(2 * Ky + 1);
  for(int k = 1; k < (Kx + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernelx_u(Kx - k) = kernelx_u(Kx + 1 - k) + exp(px[Kx - k])/(2 * Kx + 1 - 2 * k);
  }

  for(int k = 1; k < (Ky + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernely_u(Ky - k) = kernely_u(Ky + 1 - k) + exp(py[Ky - k])/(2 * Ky + 1 - 2 * k);
  }

  // normalize
  for(int k = 0; k < (Kx + 1); ++k){
    kernelx(k) = kernelx_u(k)/sumexppx;
  }

  for(int k = 0; k < (Ky + 1); ++k){
    kernely(k) = kernely_u(k)/sumexppy;
  }

  arma::mat precw = arma::eye(nws, nws) + rhox * Ax + rhoy * Ay;

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_w;
  std::vector<NumericVector> post_z;
  std::vector<NumericVector> post_phi;
  std::vector<double> post_rhox;
  std::vector<double> post_rhoy;
  std::vector<NumericVector> post_kernelx;
  std::vector<NumericVector> post_kernely;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update ws with elliptical slice
    ld_logw_noepsilon dens_logw_noepsilon(Phi, Z, contacts, nshoes, ngrid, naccidentals, w_indices);

    logw = ellipticalslice(logw, precw, nws, &dens_logw_noepsilon);

    for(int a = 0; a < nws; ++a){
      w[a] = exp(logw[a]);
    }

    if(print_output == 1){
      for(int a = 0; a < nws; ++a){
        Rcout << "w" << a << " "<< w[a] << std::endl;
      }
    }

    NumericVector denoms(nshoes);
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        denoms[d] = denoms[d] + w[w_indices[a]] * Phi[contacts(d, a)];
      }
    }

    //update Phis with slice (with a uniform prior) # updated 2018
    for(int i = 0; i < ncategories_contacts; ++i){
      double oldphi = Phi[i];
      ld_phi_noepsilon dens_phi_noepsilon(i, philocs, Z, w, w_indices, naccidentals, denoms, oldphi);
      Phi[i] = slice(Phi[i], &dens_phi_noepsilon, 0.01, 100, 0, 1);
      logPhi[i] = log(Phi[i]);
      if(print_output == 1){
        Rcout << "phi" << i << " "<< Phi[i] << std::endl;
      }
      NumericMatrix phinumb = philocs[i];
      for(int j=0; j < phinumb.nrow(); ++j){
        denoms(phinumb(j, 0)) = denoms(phinumb(j, 0)) + (Phi[i] - oldphi) * w[w_indices[phinumb(j, 1)]];
      }
    }

    // update kernelx with slice # updated 2018
    for(int i = 0; i < (Kx + 1); ++i){
      double old = px[i];
      ld_p dens_p(deviationsx, i, sumexppx, Kx, sigma, px, n, kernelx_u);
      px[i] = slice(px[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(px[i] - old) - 1.0);
      sumexppx = sumexppx + diff;
      for(int j = 0; j <= i; ++j){
        kernelx_u[j] =  kernelx_u[j] + diff/(1 + 2 * i);
      }
    }

    // update kernely with slice # updated 2018
    for(int i = 0; i < (Ky + 1); ++i){
      double old = py[i];
      ld_p dens_p(deviationsy, i, sumexppy, Ky, sigma, py, n, kernely_u);
      py[i] = slice(py[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(py[i] - old) - 1.0);
      sumexppy = sumexppy + diff;
      for(int j = 0; j <= i; ++j){
        kernely_u[j] =  kernely_u[j] + diff/(1 + 2 * i);
      }
    }

    // normalize the kernels
    for(int k = 0; k < (Kx + 1); ++k){
      kernelx(k) = kernelx_u(k)/sumexppx;
    }

    for(int k = 0; k < (Ky + 1); ++k){
      kernely(k) = kernely_u(k)/sumexppy;
    }
    if(print_output == 1){
      for(int i = 0; i < (Kx + 1); ++i){
        Rcout << "kernelx" << i << " "<< kernelx[i] << std::endl;
      }
      for(int i = 0; i < (Ky + 1); ++i){
        Rcout << "kernely" << i << " "<< kernely[i] << std::endl;
      }
    }


    //Update the gridpoint assignments via Gibbs # updated 2018
    int accind = 0;
    // update allocations shoe by shoe
    for(int d = 0; d < nshoes; ++d){
      // update by accidental
      for(int i = 0; i < naccidentals[d]; ++i){
        // remove accidental from current location
        Z(d, zoptions(zhome[accind], z[accind])) -= 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] -= 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] -= 1;
        NumericVector logwt(nzoptions[zhome[accind]]);
        for(int ind = 0; ind < nzoptions[zhome[accind]]; ++ind){
          logwt[ind] = logw[w_indices[zoptions(zhome[accind], ind)]] + logPhi(contacts(d, zoptions(zhome[accind], ind)));
          logwt[ind] += log(kernelx[optiondiffsx(zhome[accind], ind)]) + log(kernely[optiondiffsy(zhome[accind], ind)]);
        }

        // sample new index according to normalized conditional probabilities
        z[accind] = rdisc_log(logwt);
        Z(d, zoptions(zhome[accind], z[accind])) += 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] += 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] += 1;
        accind++;
      }
    }


    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_w.push_back(clone(w));
      post_z.push_back(clone(z));
      post_phi.push_back(clone(Phi));
      post_kernelx.push_back(clone(kernelx));
      post_kernely.push_back(clone(kernely));
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob_noepsilon(deviationsx, deviationsy, px, py, kernelx, kernely, sigma, Kx, Ky, ngrid, w, w_indices, logw, zhome, z, zoptions, shoenums, denoms, contacts, naccidentals, n, nshoes, logPhi, Phi, precw);
        logliks.push_back(loglike);
      }
      Rcout << "Our model no scores mcmc iteration " << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["kernelx"] = post_kernelx, _["kernely"] = post_kernely,
    _["w"]=post_w,  _["phi"] = post_phi);
}

// [[Rcpp::export]]
List mcmc_inference_noepsilonnoZ(NumericVector shoenums, // the shoe number of each accidental
                            NumericMatrix contacts, // the value of the contact surface at each gridpoint (numbers of shoes by number of atoms)
                            int ncategories_contacts, // the number of categories of contact surface value
                            int nshoes, // the number of shoes in the training set
                            NumericVector w_init, // the initializations of values of the w parameters
                            NumericVector w_indices, // the w index to which each gridpoint is assigned
                            NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                            NumericVector phi_init, // the initialization of the phi variables
                            arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                            double rhox_init, // the smoothing values for the w's in the x direction
                            double rhoy_init, // the smoothing values for the w's in the y direction
                            int nsim = 10, // the number of simulations to be run
                            int print_output = 1, // if we want additional print statements while running
                            int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                            int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w = w_init;
  int nws = w.size(); // the number of unique w variables
  int ngrid = contacts.ncol();
  arma::vec logw(nws);
  for(int a = 0; a < nws; ++a){
    logw(a) = log(w[a]);
  }
  NumericVector Phi = phi_init;
  double rhox = rhox_init;
  double rhoy = rhoy_init;
  int n = zhome.size();   // number of accidentals

  //create philocs matrix. This separates the grid into ncategories different lists,
  // each containing references to which grid are in each category. It also says which shoe
  std::vector<NumericMatrix> philocs;
  for(int j = 0; j < ncategories_contacts; ++j){
    std::vector<int> indices;
    std::vector<int> shoenums;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        if(contacts(d, a) == j){
          indices.push_back(a);
          shoenums.push_back(d);
        }
      }
    }

    NumericMatrix phij(shoenums.size(), 2);
    for(int i = 0; i < shoenums.size(); ++i){
      phij(i,0) = shoenums[i];
      phij(i,1) = indices[i];
    }
    philocs.push_back(phij);
  }

  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each atom
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = zhome[i];
    Z(shoenums[i], zhome[i])++;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }

  NumericVector logPhi(ncategories_contacts);
  for(int i= 0; i < ncategories_contacts; ++i){
    logPhi[i] = log(Phi[i]); // initialize at reasonable intervals
  }

  arma::mat wprecision = arma::eye(nws, nws) + rhox * Ax + rhoy * Ay;

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_w;
  std::vector<NumericVector> post_phi;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update ws with elliptical slice
    ld_logw_noepsilon dens_logw_noepsilon(Phi, Z, contacts, nshoes, ngrid, naccidentals, w_indices);

    logw = ellipticalslice(logw, wprecision, nws, &dens_logw_noepsilon);

    for(int a = 0; a < nws; ++a){
      w[a] = exp(logw[a]);
    }

    if(print_output == 1){
      for(int a = 0; a < nws; ++a){
        Rcout << "w" << a << " "<< w[a] << std::endl;
      }
    }

    NumericVector denoms(nshoes);
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        denoms[d] = denoms[d] + w[w_indices[a]] * Phi[contacts(d, a)];
      }
    }

    //update Phis with slice (with a uniform prior) # updated 2018
    for(int i = 0; i < ncategories_contacts; ++i){
      double oldphi = Phi[i];
      ld_phi_noepsilon dens_phi_noepsilon(i, philocs, Z, w, w_indices, naccidentals, denoms, oldphi);
      Phi[i] = slice(Phi[i], &dens_phi_noepsilon, 0.01, 100, 0, 1);
      logPhi[i] = log(Phi[i]);
      if(print_output == 1){
        Rcout << "phi" << i << " "<< Phi[i] << std::endl;
      }
      NumericMatrix phinumb = philocs[i];
      for(int j=0; j < phinumb.nrow(); ++j){
        denoms(phinumb(j, 0)) = denoms(phinumb(j, 0)) + (Phi[i] - oldphi) * w[w_indices[phinumb(j, 1)]];
      }
    }

    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_w.push_back(clone(w));
      post_phi.push_back(clone(Phi));
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob_noepsilonnoZ(w_indices, logw, zhome, shoenums, denoms, contacts, naccidentals, n, logPhi, Phi, wprecision);
        logliks.push_back(loglike);
      }
      Rcout << "Our model no kernel no scores mcmc iteration " << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["w"]=post_w,  _["phi"] = post_phi);
}

// [[Rcpp::export]]
List mcmc_inference_noPhi(NumericVector shoenums, // the shoe number of each accidental
                        int nshoes, // the number of shoes in the training set
                        NumericVector z_init, // initialization of the assignments (from 0 to Kx * ky + 1)
                        NumericVector w_init, // the initializations of values of the w parameters
                        NumericVector w_indices, // the w index to which each gridpoint is assigned
                        NumericVector zhome, // the gridpoint index of the locations of each of the observed accidentals
                        NumericMatrix zoptions, // the possible grid for each accidental to be assigned to
                        NumericVector nzoptions, // the number of possible grid for each accidental assignment
                        NumericMatrix optiondiffsx, // the differences in x for each of the options
                        NumericMatrix optiondiffsy,// the differences in y for each of the options
                        arma::mat Ax, arma::mat Ay, // the graph defining connections between ws in x and y space, respectively
                        int Kx, int Ky, // the kernel width in the x and y direction
                        NumericVector px_init, NumericVector py_init, // the initial unnormalized log kernel values in the x and y direction
                        double sigma, // the hyperparameter for the lognormal for the kernel values
                        NumericVector u_init, // the initialization of the u auxiliary variable
                        double q_init, // the initialization for q hyperparameter
                        double qa, double qb, // the alpha and beta parameters for the prior on q
                        double rhox_init, // the smoothing values for the w's in the x direction
                        double rhoy_init, // the smoothing values for the w's in the y direction
                        int nsim = 10, // the number of simulations to be run
                        int print_output = 1, // if we want additional print statements while running
                        int nburn = 0, int nthin = 1, // the burn-in and thinning levels
                        int status = 1 // the thinning of the joint probability calculation)
){
  // initialize each variable
  NumericVector w = w_init;
  int nws = w.size(); // the number of unique w variables
  int ngrid = w_indices.size();
  arma::vec logw(nws);
  for(int a = 0; a < nws; ++a){
    logw(a) = log(w[a]);
  }
  NumericVector px = px_init;
  NumericVector py = py_init;
  NumericVector U = u_init;
  double q = q_init;
  double rhox = rhox_init;
  double rhoy = rhoy_init;
  int n = z_init.size();   // number of accidentals

  NumericVector z(n); // the gridpoint assignments for each shoe
  NumericMatrix Z(nshoes, ngrid); //the number of accidentals from each shoe assigned to each atom
  // initialize the above to objects
  for(int i = 0; i < n; ++i){
    z[i] = z_init[i];
    Z(shoenums[i], zoptions(zhome[i], z[i]))++;
  }

  //find grid of each point,
  // and keep track of how much they deviate from their "home" gridpoint in each direction
  NumericVector deviationsx(Kx + 1); // the counts of how many of each deviation
  NumericVector deviationsy(Ky + 1);

  for(int i = 0; i < n; ++i){
    deviationsx[optiondiffsx(zhome[i], z[i])] += 1;
    deviationsy[optiondiffsy(zhome[i], z[i])] += 1;
  }

  NumericVector naccidentals(nshoes); // the number of accidentals on each shoe

  for(int i=0; i < shoenums.size(); ++i){
    naccidentals[shoenums[i]]++;
  }

  // it is convenient to also maintain some parameters on the log scale as well
  NumericVector logU(nshoes);
  for(int i=0; i < nshoes; ++i){
    logU[i] = log(U[i]);
  }

  // the normalization for the kernels
  double sumexppx = 0;
  double sumexppy = 0;
  for(int i = 0; i < (Kx + 1); ++i){
    sumexppx += exp(px[i]);
  }
  for(int i = 0; i < (Ky + 1); ++i){
    sumexppy += exp(py[i]);
  }

  // The kernels on the non-normalized scale
  // unnormalized
  NumericVector kernelx_u(Kx + 1);
  NumericVector kernely_u(Ky + 1);
  // normalized
  NumericVector kernelx(Kx + 1);
  NumericVector kernely(Ky + 1);

  // fill them
  kernelx_u[Kx] = exp(px[Kx])/(2 * Kx + 1);
  kernely_u[Ky] = exp(py[Ky])/(2 * Ky + 1);
  for(int k = 1; k < (Kx + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernelx_u(Kx - k) = kernelx_u(Kx + 1 - k) + exp(px[Kx - k])/(2 * Kx + 1 - 2 * k);
  }

  for(int k = 1; k < (Ky + 1); ++k){
    //iteratively build up the kernel probabilities starting with the smallest
    kernely_u(Ky - k) = kernely_u(Ky + 1 - k) + exp(py[Ky - k])/(2 * Ky + 1 - 2 * k);
  }

  // normalize
  for(int k = 0; k < (Kx + 1); ++k){
    kernelx(k) = kernelx_u(k)/sumexppx;
  }

  for(int k = 0; k < (Ky + 1); ++k){
    kernely(k) = kernely_u(k)/sumexppy;
  }

  arma::mat precw = arma::eye(nws, nws) + rhox * Ax + rhoy * Ay;

  // traces
  // these store posterior samples
  std::vector<NumericVector> post_w;
  std::vector<NumericVector> post_U;
  std::vector<NumericVector> post_z;
  std::vector<double> post_q;
  std::vector<double> post_rhox;
  std::vector<double> post_rhoy;
  std::vector<NumericVector> post_kernelx;
  std::vector<NumericVector> post_kernely;
  std::vector<NumericVector> logliks;

  std::time_t  timev = time(0);

  for(size_t t=0; t<nsim; ++t) {
    Rcpp::checkUserInterrupt();

    //update u with slice # updated 2018 - note that step out width might need to change
    for(int d = 0; d < nshoes; ++d){

      ld_u_noPhi dens_u(q, naccidentals, d, w, Z, ngrid, w_indices);
      // gave an upper bound because I was running into numerical problems sometimes early in the chain
      U[d] = slice(U[d], &dens_u, 20 * (sqrt(naccidentals[d]))/(ngrid * q), 100, 0, 100000);
      logU[d] = log(U[d]);
    }

    //update ws with elliptical slice
    ld_logw_noPhi dens_logw(q, U, Z, nshoes, ngrid, w_indices);

    logw = ellipticalslice(logw, precw, nws, &dens_logw);

    for(int a = 0; a < nws; ++a){
      w[a] = exp(logw[a]);
    }

    if(print_output == 1){
      for(int a = 0; a < nws; ++a){
        Rcout << "w" << a << " "<< w[a] << std::endl;
      }
    }

    //update q with slice (with gamma(qa, qb) prior)
    int highestcount = max(Z);
    NumericVector counts(highestcount + 1);
    for(int ww = 0; ww < (highestcount + 1); ++ww){
      counts[ww] = 0;
    }

    double logsumdenom = 0.0;
    for(int d = 0; d < nshoes; ++d){
      for(int a = 0; a < ngrid; ++a){
        counts[Z(d, a)]++;
        logsumdenom = logsumdenom + log(U[d] * w[w_indices[a]] + 1);
      }
    }

    ld_q dens_q(logsumdenom, qa, qb, highestcount, counts);
    q = slice(q, &dens_q, 0.2, 100, 0, 500);
    if(print_output == 1){
      Rcout << "q " << q << std::endl;
    }

    // update kernelx with slice # updated 2018
    for(int i = 0; i < (Kx + 1); ++i){
      double old = px[i];
      ld_p dens_p(deviationsx, i, sumexppx, Kx, sigma, px, n, kernelx_u);
      px[i] = slice(px[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(px[i] - old) - 1.0);
      sumexppx = sumexppx + diff;
      for(int j = 0; j <= i; ++j){
        kernelx_u[j] =  kernelx_u[j] + diff/(1 + 2 * i);
      }
    }

    // update kernely with slice # updated 2018
    for(int i = 0; i < (Ky + 1); ++i){
      double old = py[i];
      ld_p dens_p(deviationsy, i, sumexppy, Ky, sigma, py, n, kernely_u);
      py[i] = slice(py[i], &dens_p, 1, 100, -30, 30);
      double diff = exp(old) * (exp(py[i] - old) - 1.0);
      sumexppy = sumexppy + diff;
      for(int j = 0; j <= i; ++j){
        kernely_u[j] =  kernely_u[j] + diff/(1 + 2 * i);
      }
    }

    // normalize the kernels
    for(int k = 0; k < (Kx + 1); ++k){
      kernelx(k) = kernelx_u(k)/sumexppx;
    }

    for(int k = 0; k < (Ky + 1); ++k){
      kernely(k) = kernely_u(k)/sumexppy;
    }
    if(print_output == 1){
      for(int i = 0; i < (Kx + 1); ++i){
        Rcout << "kernelx" << i << " "<< kernelx[i] << std::endl;
      }
      for(int i = 0; i < (Ky + 1); ++i){
        Rcout << "kernely" << i << " "<< kernely[i] << std::endl;
      }
    }


    //Update the gridpoint assignments via Gibbs
    int accind = 0;
    // update allocations shoe by shoe
    for(int d = 0; d < nshoes; ++d){
      // update by accidental
      for(int i = 0; i < naccidentals[d]; ++i){
        // remove accidental from current location
        Z(d, zoptions(zhome[accind], z[accind])) -= 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] -= 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] -= 1;
        NumericVector logwt(nzoptions[zhome[accind]]);
        for(int ind = 0; ind < nzoptions[zhome[accind]]; ++ind){
          logwt[ind] = logw[w_indices[zoptions(zhome[accind], ind)]] + log(q + Z(d, zoptions(zhome[accind], ind))) - log(U[d] * w[w_indices[zoptions(zhome[accind], ind)]] + 1.0);
          logwt[ind] += log(kernelx[optiondiffsx(zhome[accind], ind)]) + log(kernely[optiondiffsy(zhome[accind], ind)]);
        }

        // sample new index according to normalized conditional probabilities
        z[accind] = rdisc_log(logwt);
        Z(d, zoptions(zhome[accind], z[accind])) += 1;
        deviationsx[optiondiffsx(zhome[accind], z[accind])] += 1;
        deviationsy[optiondiffsy(zhome[accind], z[accind])] += 1;
        accind++;
      }
    }


    // store trace
    if((t >= nburn) & (t % nthin==0)) {
      post_w.push_back(clone(w));
      post_U.push_back(clone(U));
      post_z.push_back(clone(z));
      post_q.push_back(q);
      post_kernelx.push_back(clone(kernelx));
      post_kernely.push_back(clone(kernely));
    }
    //monitor the trace
    if(t % status == 0){
      //calculate full joint probability
      if(t % 10 == 0){
        NumericVector loglike = logprob_noPhi(deviationsx, deviationsy, px, py, kernelx, kernely, sigma, Kx, Ky, ngrid, w, w_indices, logw, U, q, Z, naccidentals, n, nshoes, qa, qb, precw);
        logliks.push_back(loglike);
      }
      Rcout << "Our model constant phi mcmc iteration " << t+1 << std::endl;
    }

  }
  return List::create(
    _["logliks"]= logliks, _["kernelx"] = post_kernelx, _["kernely"] = post_kernely,
    _["w"]=post_w,  _["U"]=post_U, _["q"]=post_q);
}
