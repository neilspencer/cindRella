// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// heldoutprob_importance
NumericVector heldoutprob_importance(NumericVector contact, int ncategories_contacts, NumericMatrix w, NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, NumericVector q, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, int nsim_per, int status);
RcppExport SEXP _cindRella_heldoutprob_importance(SEXP contactSEXP, SEXP ncategories_contactsSEXP, SEXP wSEXP, SEXP PhiSEXP, SEXP KxSEXP, SEXP KySEXP, SEXP kernelxSEXP, SEXP kernelySEXP, SEXP qSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP nsim_perSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernelx(kernelxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernely(kernelySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< int >::type nsim_per(nsim_perSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance(contact, ncategories_contacts, w, Phi, Kx, Ky, kernelx, kernely, q, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, nsim_per, status));
    return rcpp_result_gen;
END_RCPP
}
// heldoutprob_importance_now
NumericVector heldoutprob_importance_now(NumericVector contact, int ncategories_contacts, NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, NumericVector q, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, int nsim_per, int status);
RcppExport SEXP _cindRella_heldoutprob_importance_now(SEXP contactSEXP, SEXP ncategories_contactsSEXP, SEXP PhiSEXP, SEXP KxSEXP, SEXP KySEXP, SEXP kernelxSEXP, SEXP kernelySEXP, SEXP qSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP nsim_perSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernelx(kernelxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernely(kernelySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< int >::type nsim_per(nsim_perSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance_now(contact, ncategories_contacts, Phi, Kx, Ky, kernelx, kernely, q, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, nsim_per, status));
    return rcpp_result_gen;
END_RCPP
}
// heldoutprob_importance_noZ
NumericVector heldoutprob_importance_noZ(NumericVector contact, int ncategories_contacts, NumericMatrix w, NumericMatrix Phi, NumericVector q, NumericVector w_indices, NumericVector zhome, int nsim_per, int status);
RcppExport SEXP _cindRella_heldoutprob_importance_noZ(SEXP contactSEXP, SEXP ncategories_contactsSEXP, SEXP wSEXP, SEXP PhiSEXP, SEXP qSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP nsim_perSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< int >::type nsim_per(nsim_perSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance_noZ(contact, ncategories_contacts, w, Phi, q, w_indices, zhome, nsim_per, status));
    return rcpp_result_gen;
END_RCPP
}
// heldoutprob_importance_noepsilon
NumericVector heldoutprob_importance_noepsilon(NumericVector contact, int ncategories_contacts, NumericMatrix w, NumericMatrix Phi, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, int status);
RcppExport SEXP _cindRella_heldoutprob_importance_noepsilon(SEXP contactSEXP, SEXP ncategories_contactsSEXP, SEXP wSEXP, SEXP PhiSEXP, SEXP KxSEXP, SEXP KySEXP, SEXP kernelxSEXP, SEXP kernelySEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernelx(kernelxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernely(kernelySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance_noepsilon(contact, ncategories_contacts, w, Phi, Kx, Ky, kernelx, kernely, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, status));
    return rcpp_result_gen;
END_RCPP
}
// heldoutprob_importance_noepsilonnoZ
NumericVector heldoutprob_importance_noepsilonnoZ(NumericVector contact, int ncategories_contacts, NumericMatrix w, NumericMatrix Phi, NumericVector w_indices, NumericVector zhome, int status);
RcppExport SEXP _cindRella_heldoutprob_importance_noepsilonnoZ(SEXP contactSEXP, SEXP ncategories_contactsSEXP, SEXP wSEXP, SEXP PhiSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance_noepsilonnoZ(contact, ncategories_contacts, w, Phi, w_indices, zhome, status));
    return rcpp_result_gen;
END_RCPP
}
// heldoutprob_importance_noPhi
NumericVector heldoutprob_importance_noPhi(NumericVector contact, NumericMatrix w, int Kx, int Ky, NumericMatrix kernelx, NumericMatrix kernely, NumericVector q, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, int nsim_per, int status);
RcppExport SEXP _cindRella_heldoutprob_importance_noPhi(SEXP contactSEXP, SEXP wSEXP, SEXP KxSEXP, SEXP KySEXP, SEXP kernelxSEXP, SEXP kernelySEXP, SEXP qSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP nsim_perSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernelx(kernelxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernely(kernelySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< int >::type nsim_per(nsim_perSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(heldoutprob_importance_noPhi(contact, w, Kx, Ky, kernelx, kernely, q, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, nsim_per, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference
List mcmc_inference(NumericVector shoenums, NumericMatrix contacts, int ncategories_contacts, int nshoes, NumericVector z_init, NumericVector w_init, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, NumericVector phi_init, arma::mat Ax, arma::mat Ay, int Kx, int Ky, NumericVector px_init, NumericVector py_init, double sigma, NumericVector u_init, double q_init, double qa, double qb, double rhox_init, double rhoy_init, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference(SEXP shoenumsSEXP, SEXP contactsSEXP, SEXP ncategories_contactsSEXP, SEXP nshoesSEXP, SEXP z_initSEXP, SEXP w_initSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP phi_initSEXP, SEXP AxSEXP, SEXP AySEXP, SEXP KxSEXP, SEXP KySEXP, SEXP px_initSEXP, SEXP py_initSEXP, SEXP sigmaSEXP, SEXP u_initSEXP, SEXP q_initSEXP, SEXP qaSEXP, SEXP qbSEXP, SEXP rhox_initSEXP, SEXP rhoy_initSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contacts(contactsSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_init(w_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type px_init(px_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py_init(py_initSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_init(u_initSEXP);
    Rcpp::traits::input_parameter< double >::type q_init(q_initSEXP);
    Rcpp::traits::input_parameter< double >::type qa(qaSEXP);
    Rcpp::traits::input_parameter< double >::type qb(qbSEXP);
    Rcpp::traits::input_parameter< double >::type rhox_init(rhox_initSEXP);
    Rcpp::traits::input_parameter< double >::type rhoy_init(rhoy_initSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference(shoenums, contacts, ncategories_contacts, nshoes, z_init, w_init, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, phi_init, Ax, Ay, Kx, Ky, px_init, py_init, sigma, u_init, q_init, qa, qb, rhox_init, rhoy_init, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference_noZ
List mcmc_inference_noZ(NumericVector shoenums, NumericMatrix contacts, int ncategories_contacts, int nshoes, NumericVector w_init, NumericVector w_indices, NumericVector zhome, NumericVector phi_init, arma::mat Ax, arma::mat Ay, NumericVector u_init, double q_init, double qa, double qb, double rhox_init, double rhoy_init, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference_noZ(SEXP shoenumsSEXP, SEXP contactsSEXP, SEXP ncategories_contactsSEXP, SEXP nshoesSEXP, SEXP w_initSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP phi_initSEXP, SEXP AxSEXP, SEXP AySEXP, SEXP u_initSEXP, SEXP q_initSEXP, SEXP qaSEXP, SEXP qbSEXP, SEXP rhox_initSEXP, SEXP rhoy_initSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contacts(contactsSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_init(w_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_init(u_initSEXP);
    Rcpp::traits::input_parameter< double >::type q_init(q_initSEXP);
    Rcpp::traits::input_parameter< double >::type qa(qaSEXP);
    Rcpp::traits::input_parameter< double >::type qb(qbSEXP);
    Rcpp::traits::input_parameter< double >::type rhox_init(rhox_initSEXP);
    Rcpp::traits::input_parameter< double >::type rhoy_init(rhoy_initSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference_noZ(shoenums, contacts, ncategories_contacts, nshoes, w_init, w_indices, zhome, phi_init, Ax, Ay, u_init, q_init, qa, qb, rhox_init, rhoy_init, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference_now
List mcmc_inference_now(NumericVector shoenums, NumericMatrix contacts, int ncategories_contacts, int nshoes, NumericVector z_init, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, NumericVector phi_init, arma::mat Ax, arma::mat Ay, int Kx, int Ky, NumericVector px_init, NumericVector py_init, double sigma, NumericVector u_init, double q_init, double qa, double qb, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference_now(SEXP shoenumsSEXP, SEXP contactsSEXP, SEXP ncategories_contactsSEXP, SEXP nshoesSEXP, SEXP z_initSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP phi_initSEXP, SEXP AxSEXP, SEXP AySEXP, SEXP KxSEXP, SEXP KySEXP, SEXP px_initSEXP, SEXP py_initSEXP, SEXP sigmaSEXP, SEXP u_initSEXP, SEXP q_initSEXP, SEXP qaSEXP, SEXP qbSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contacts(contactsSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type px_init(px_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py_init(py_initSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_init(u_initSEXP);
    Rcpp::traits::input_parameter< double >::type q_init(q_initSEXP);
    Rcpp::traits::input_parameter< double >::type qa(qaSEXP);
    Rcpp::traits::input_parameter< double >::type qb(qbSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference_now(shoenums, contacts, ncategories_contacts, nshoes, z_init, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, phi_init, Ax, Ay, Kx, Ky, px_init, py_init, sigma, u_init, q_init, qa, qb, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference_noepsilon
List mcmc_inference_noepsilon(NumericVector shoenums, NumericMatrix contacts, int ncategories_contacts, int nshoes, NumericVector z_init, NumericVector w_init, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, NumericVector phi_init, arma::mat Ax, arma::mat Ay, int Kx, int Ky, NumericVector px_init, NumericVector py_init, double sigma, double rhox_init, double rhoy_init, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference_noepsilon(SEXP shoenumsSEXP, SEXP contactsSEXP, SEXP ncategories_contactsSEXP, SEXP nshoesSEXP, SEXP z_initSEXP, SEXP w_initSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP phi_initSEXP, SEXP AxSEXP, SEXP AySEXP, SEXP KxSEXP, SEXP KySEXP, SEXP px_initSEXP, SEXP py_initSEXP, SEXP sigmaSEXP, SEXP rhox_initSEXP, SEXP rhoy_initSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contacts(contactsSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_init(w_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type px_init(px_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py_init(py_initSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rhox_init(rhox_initSEXP);
    Rcpp::traits::input_parameter< double >::type rhoy_init(rhoy_initSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference_noepsilon(shoenums, contacts, ncategories_contacts, nshoes, z_init, w_init, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, phi_init, Ax, Ay, Kx, Ky, px_init, py_init, sigma, rhox_init, rhoy_init, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference_noepsilonnoZ
List mcmc_inference_noepsilonnoZ(NumericVector shoenums, NumericMatrix contacts, int ncategories_contacts, int nshoes, NumericVector w_init, NumericVector w_indices, NumericVector zhome, NumericVector phi_init, arma::mat Ax, arma::mat Ay, double rhox_init, double rhoy_init, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference_noepsilonnoZ(SEXP shoenumsSEXP, SEXP contactsSEXP, SEXP ncategories_contactsSEXP, SEXP nshoesSEXP, SEXP w_initSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP phi_initSEXP, SEXP AxSEXP, SEXP AySEXP, SEXP rhox_initSEXP, SEXP rhoy_initSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contacts(contactsSEXP);
    Rcpp::traits::input_parameter< int >::type ncategories_contacts(ncategories_contactsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_init(w_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< double >::type rhox_init(rhox_initSEXP);
    Rcpp::traits::input_parameter< double >::type rhoy_init(rhoy_initSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference_noepsilonnoZ(shoenums, contacts, ncategories_contacts, nshoes, w_init, w_indices, zhome, phi_init, Ax, Ay, rhox_init, rhoy_init, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_inference_noPhi
List mcmc_inference_noPhi(NumericVector shoenums, int nshoes, NumericVector z_init, NumericVector w_init, NumericVector w_indices, NumericVector zhome, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy, arma::mat Ax, arma::mat Ay, int Kx, int Ky, NumericVector px_init, NumericVector py_init, double sigma, NumericVector u_init, double q_init, double qa, double qb, double rhox_init, double rhoy_init, int nsim, int print_output, int nburn, int nthin, int status);
RcppExport SEXP _cindRella_mcmc_inference_noPhi(SEXP shoenumsSEXP, SEXP nshoesSEXP, SEXP z_initSEXP, SEXP w_initSEXP, SEXP w_indicesSEXP, SEXP zhomeSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP, SEXP AxSEXP, SEXP AySEXP, SEXP KxSEXP, SEXP KySEXP, SEXP px_initSEXP, SEXP py_initSEXP, SEXP sigmaSEXP, SEXP u_initSEXP, SEXP q_initSEXP, SEXP qaSEXP, SEXP qbSEXP, SEXP rhox_initSEXP, SEXP rhoy_initSEXP, SEXP nsimSEXP, SEXP print_outputSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP statusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shoenums(shoenumsSEXP);
    Rcpp::traits::input_parameter< int >::type nshoes(nshoesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_init(w_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zhome(zhomeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ax(AxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ay(AySEXP);
    Rcpp::traits::input_parameter< int >::type Kx(KxSEXP);
    Rcpp::traits::input_parameter< int >::type Ky(KySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type px_init(px_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py_init(py_initSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_init(u_initSEXP);
    Rcpp::traits::input_parameter< double >::type q_init(q_initSEXP);
    Rcpp::traits::input_parameter< double >::type qa(qaSEXP);
    Rcpp::traits::input_parameter< double >::type qb(qbSEXP);
    Rcpp::traits::input_parameter< double >::type rhox_init(rhox_initSEXP);
    Rcpp::traits::input_parameter< double >::type rhoy_init(rhoy_initSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type status(statusSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_inference_noPhi(shoenums, nshoes, z_init, w_init, w_indices, zhome, zoptions, nzoptions, optiondiffsx, optiondiffsy, Ax, Ay, Kx, Ky, px_init, py_init, sigma, u_init, q_init, qa, qb, rhox_init, rhoy_init, nsim, print_output, nburn, nthin, status));
    return rcpp_result_gen;
END_RCPP
}
// plot_posterior
List plot_posterior(NumericVector contact, NumericMatrix w, NumericMatrix Phi, NumericMatrix kernelx, NumericMatrix kernely, NumericVector q, NumericVector w_indices, NumericMatrix zoptions, NumericVector nzoptions, NumericMatrix optiondiffsx, NumericMatrix optiondiffsy);
RcppExport SEXP _cindRella_plot_posterior(SEXP contactSEXP, SEXP wSEXP, SEXP PhiSEXP, SEXP kernelxSEXP, SEXP kernelySEXP, SEXP qSEXP, SEXP w_indicesSEXP, SEXP zoptionsSEXP, SEXP nzoptionsSEXP, SEXP optiondiffsxSEXP, SEXP optiondiffsySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type contact(contactSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernelx(kernelxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernely(kernelySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_indices(w_indicesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zoptions(zoptionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzoptions(nzoptionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsx(optiondiffsxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type optiondiffsy(optiondiffsySEXP);
    rcpp_result_gen = Rcpp::wrap(plot_posterior(contact, w, Phi, kernelx, kernely, q, w_indices, zoptions, nzoptions, optiondiffsx, optiondiffsy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cindRella_heldoutprob_importance", (DL_FUNC) &_cindRella_heldoutprob_importance, 17},
    {"_cindRella_heldoutprob_importance_now", (DL_FUNC) &_cindRella_heldoutprob_importance_now, 15},
    {"_cindRella_heldoutprob_importance_noZ", (DL_FUNC) &_cindRella_heldoutprob_importance_noZ, 9},
    {"_cindRella_heldoutprob_importance_noepsilon", (DL_FUNC) &_cindRella_heldoutprob_importance_noepsilon, 15},
    {"_cindRella_heldoutprob_importance_noepsilonnoZ", (DL_FUNC) &_cindRella_heldoutprob_importance_noepsilonnoZ, 7},
    {"_cindRella_heldoutprob_importance_noPhi", (DL_FUNC) &_cindRella_heldoutprob_importance_noPhi, 15},
    {"_cindRella_mcmc_inference", (DL_FUNC) &_cindRella_mcmc_inference, 31},
    {"_cindRella_mcmc_inference_noZ", (DL_FUNC) &_cindRella_mcmc_inference_noZ, 21},
    {"_cindRella_mcmc_inference_now", (DL_FUNC) &_cindRella_mcmc_inference_now, 27},
    {"_cindRella_mcmc_inference_noepsilon", (DL_FUNC) &_cindRella_mcmc_inference_noepsilon, 27},
    {"_cindRella_mcmc_inference_noepsilonnoZ", (DL_FUNC) &_cindRella_mcmc_inference_noepsilonnoZ, 17},
    {"_cindRella_mcmc_inference_noPhi", (DL_FUNC) &_cindRella_mcmc_inference_noPhi, 28},
    {"_cindRella_plot_posterior", (DL_FUNC) &_cindRella_plot_posterior, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_cindRella(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
