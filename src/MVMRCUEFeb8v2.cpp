#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;



class Options_MVMRCUE{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_MVMRCUE(uword J){
    this -> agM = 0.001*ones(J, 1);
    this -> bgM = 0.001*ones(J, 1);
    this -> atau1 = 0.001;
    this -> btau1 = 0.001;
    this -> atau2 = 0.001;
    this -> btau2 = 0.001;
    this -> a = 1;
    this -> b = 1;
    this -> maxIter = 3000;
    this -> thin = 10;
    this -> burnin = 2000;
  }
  
  Options_MVMRCUE(vec agM, vec bgM, double atau1, double btau1, double atau2, double btau2, 
                  double a, double b, uword maxIter, uword thin, uword burnin){
    
    this -> agM = agM;
    this -> bgM = bgM;
    this -> atau1 = atau1;
    this -> btau1 = btau1;
    this -> atau2 = atau2;
    this -> btau2 = btau2;
    this -> a = a;
    this -> b = b;
    this -> maxIter = maxIter;
    this -> thin = thin;
    this -> burnin = burnin;
    
  }
  
  vec agM;
  vec bgM;
  double atau1;
  double btau1;
  double atau2;
  double btau2;
  double a;
  double b;
  uword maxIter;
  uword thin;
  uword burnin;
};
struct ObjMVMRCUE{
  vec EtaIterRate;
  mat Beta1res;
  mat Beta2res;
  mat Sgga2res;
  vec Tau12res;
  vec Tau22res;
};


//[[Rcpp::export]]
double normal_pdf(double x, double m, double s)
{
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

//[[Rcpp::export]]
List MVMRCUEfun(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                arma::vec &se2){
  
  int J = gammahM.n_cols;
  int p = gammahM.n_rows;
  
  // hyper-prior parameters.
  vec agM = 0.001*ones(J, 1);
  vec bgM = 0.001*ones(J, 1);
  double atau1 = 0.001;
  double btau1 = 0.001;
  double atau2 = 0.001;
  double btau2 = 0.001;
  double a = 1;
  double b = 1;
  int maxIter = 4000; int burnin = 1000; int thin = 10;
  
  // initial values.
  ivec Eta = zeros<ivec>(p, 1);
  vec Beta1 = 0.01*ones(J, 1);
  vec Beta2 = 0.01*ones(J, 1);
  vec Sgga2 = 0.01*ones(J, 1);
  double tau12 = 0.01;
  double tau22 = 0.01;
  double w = 0.1;
  double xi2 = 0.01;
  double tau12xi2 = tau12*xi2;
  
  mat Mu = 0.01*ones(p, J);
  vec mut = 0.01*ones(p, 1);
  mat Bgamma = zeros(p, J);
  mat Bgammat = zeros(p, J);
  
  
  
  int numsave = maxIter / thin;
  mat Beta1res = ones(numsave, J);
  mat Beta2res = ones(numsave, J);
  mat Sgga2res = ones(numsave, J);
  vec Tau12res = ones(numsave, 1);
  vec Tau22res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  vec EtaIterRate = zeros(p, 1);
  
  vec GinvsG2 = Gammah / (se2 % se2);
  mat ginvsg2M = zeros(p, J);
  
  for(int j = 0; j < J; j++){
    ginvsg2M.col(j) = gammahM.col(j) / (se1M.col(j) % se1M.col(j));
  }
  mat Invsg2 = 1. / (se1M % se1M);
  vec invsG2 = 1. / (se2 % se2);
  
  int l = 0;
  // cout << "MU: " << Mu.row(0) << endl;
  
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    
    vec Invsgga2 = 1. / Sgga2;
    double invtau12xi2 = 1. / (tau12xi2);
    double invtau22 = 1. / tau22;
    
    for(int j = 0; j < J; j++){
      // cout << "j: "<< j << "Beta1: " <<Beta1[j] <<"Mu: "<<  sum(Mu.col(j)) << endl;
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    // 
    // cout << "Bgamma: sum: "<< sum(Bgamma, 0) << endl;
    // cout << "Bgammat: sum: "<<  sum(Bgammat, 0) << endl;
    // --------------------------------------- //
    // Update Gamma
    // --------------------------------------- //
    vec v20t = 1. / (invsG2 + invtau12xi2);
    vec v21t = 1. / (invsG2 + invtau22);
    
    vec mutm0t = (GinvsG2 + sum(Bgamma, 1) * invtau12xi2) % v20t;
    vec mutm1t = (GinvsG2 + sum(Bgammat, 1) * invtau22) % v21t;
    
    // cout << "mutm0t: "<< mutm0t.subvec(0, 2).t() << " sum: " << sum(mutm0t) << endl;
    // cout << "mutm1t: "<< mutm1t.subvec(0, 2).t() << " sum: " << sum(mutm1t) << endl;
    
    for(int k = 0; k < p; ++k) {
      if(Eta[k] == 1) {
        // mut[k] = mutm1t[k];
        mut[k] = mutm1t[k] + randn()*sqrt(v21t[k]);
      } else if(Eta[k] == 0) {
        // mut[k] = mutm0t[k];
        mut[k] = mutm0t[k] + randn()*sqrt(v20t[k]);
      }
    }
    
    // cout << "mut: " << sum(mut) << endl;
    // --------------------------------------- //
    // Update gamma for each j in 1,...,J
    // --------------------------------------- //
    mat V20M(p, J, fill::none);
    mat V21M(p, J, fill::none);
    mat Mutm0(p, J, fill::none);
    mat Mutm1(p, J, fill::none);
    
    for(int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      V20M.col(j) = 1.0 / (Invsgga2[j] + Invsg2.col(j) + Beta1[j] * Beta1[j] * invtau12xi2);
      V21M.col(j) = 1.0 / (Invsgga2[j] + Invsg2.col(j) + Beta2[j] * Beta2[j] * invtau22);
      vec tmp1 = Gminsb1g + Beta1[j] * Mu.col(j);
      vec tmp2 = Gminsb2g + Beta2[j] * Mu.col(j);
      Mutm0.col(j) = (ginvsg2M.col(j) + tmp1 * Beta1[j] * invtau12xi2) % V20M.col(j);
      Mutm1.col(j) = (ginvsg2M.col(j) + tmp2 * Beta2[j] * invtau22) % V21M.col(j);
      
      for(int k = 0; k < p; ++k) {
        if(Eta[k] == 1) {
          // Mu(k, j) = Mutm1(k, j);
          Mu(k, j) = Mutm1(k, j) + randn()*sqrt(V21M(k, j));
        } else {
          // Mu(k, j) = Mutm0(k, j);
          Mu(k, j) = Mutm0(k, j) + randn()*sqrt(V20M(k, j));
        }
      }
      
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    
    // cout << "Mu: "<< sum(Mu, 0) << endl; 
    // --------------------------------------- //
    // Update Beta1 and Beta2
    // --------------------------------------- //
    for(int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec muEta0j = (1 - Eta) % Mu.col(j);
      vec muEta1j = Eta % Mu.col(j);
      vec tmpmuj = Gminsb1g + Beta1[j] * Mu.col(j);
      vec tmpmujt = Gminsb2g + Beta2[j] * Mu.col(j);
      
      if(sum(Eta) == p) {
        Beta1[j] = 0;
      } else {
        double sig2b0j = 1.0 / (sum(muEta0j % Mu.col(j)) * invtau12xi2);
        double mub0j = sum(muEta0j % tmpmuj) * invtau12xi2 * sig2b0j;
        // Beta1[j] = mub0j;
        Beta1[j] = mub0j + randn()*sqrt(sig2b0j);
      }
      
      if(sum(Eta) == 0) {
        Beta2[j] = 0;
      } else {
        double sig2b1j = 1.0 / (sum(muEta1j % Mu.col(j)) * invtau22);
        double mub1j = sum(muEta1j % tmpmujt) * invtau22 * sig2b1j;
        // Beta2[j] = mub1j;
        Beta2[j] = mub1j + randn()*sqrt(sig2b1j);
      }
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    // cout << "Beta1: "<< Beta1.t() << endl;
    // cout << "Beta2: " << Beta2.t() << endl;
    // --------------------------------------- //
    // Update the variance terms
    // --------------------------------------- //
    for(int j = 0; j < J; ++j) {
      double tagmj = agM[j] + p / 2.0;
      double tbgmj = bgM[j] + sum(Mu.col(j) % Mu.col(j)) / 2.0;
      Sgga2[j] = 1.0 / randg<double>(distr_param(tagmj, 1.0 / tbgmj));
      // Sgga2[j] = tbgmj;
    }
    
    vec Gminsb1g = mut - sum(Bgamma, 1);
    vec Gminsb2g = mut - sum(Bgammat, 1);
    vec err02 = (1 - Eta) % Gminsb1g % Gminsb1g;
    vec err12 = Eta % Gminsb2g % Gminsb2g;
    
    double tatau1 = atau1 + sum(1 - Eta) / 2.0;
    double tbtau1 = btau1 + sum(err02) / (2.0 * xi2);
    tau12 = 1.0 / randg<double>(distr_param(tatau1, 1.0 / tbtau1));
    // tau12= tbtau1;
    
    
    double tatau2 = atau2 + sum(Eta) / 2.0;
    double tbtau2 = btau2 + sum(err12) / 2.0;
    tau22 = 1.0 / randg<double>(distr_param(tatau2, 1.0 / tbtau2));
    // tau22 = tbtau2;
    
    double taxi2 = 0.5 * sum(1 - Eta);
    double tbxi2 = 0.5 * sum(err02) / tau12;
    xi2 = 1.0 / randg<double>(distr_param(taxi2, 1.0 / tbxi2));
    // xi2 = tbxi2;
    
    tau12xi2 = tau12*xi2;
    if(tau12xi2 < 1e-7){
      tau12xi2 = 1e-7;
    }
    
    // cout << "Sgga2: "<< Sgga2.t() << endl;
    // cout << "tau12: "<< tau12<< endl;
    // cout << "tau22: "<< tau22 << endl;
    // cout <<"xi2: "<< xi2 << endl;
    // --------------------------------------- //
    // Update omega
    // --------------------------------------- //
    double wpa1 = a + sum(Eta);
    double wpa2 = b + p - sum(Eta);
    w = R::rbeta(wpa1, wpa2);
    
    // --------------------------------------- //
    // Update Eta
    // --------------------------------------- //
    double tau1 = sqrt(tau12 * xi2);
    double tau2 = sqrt(tau22);
    vec bmu0 = sum(Bgamma, 1);
    vec bmu1 = sum(Bgammat, 1);
    // vec tmp = zeros(p, 1);
    for(int m = 0; m < p; ++m) {
      double p1 = w*normal_pdf(mut[m], bmu1[m], tau2);
      double p0 = (1 - w) * normal_pdf(mut[m], bmu0[m], tau1);
      double prob = p1 / (p1 + p0);
      Eta[m] = R::rbinom(1, prob);
      // tmp[m] = prob;
    }
    
    // cout << "tmp: "<< tmp.subvec(0, 4).t() << " sum: "<< sum(tmp) << endl;
    // cout << "Eta: " << sum(Eta) << endl;
    // 
    // cout << "=========================================="<<endl;
    
    // Store results
    if(iter >= burnin) {
      if((iter - burnin) % thin == 0) {
        Sgga2res.row(l) = Sgga2.t();
        Beta1res.row(l) = Beta1.t();
        Beta2res.row(l) = Beta2.t();
        Tau12res[l] = tau12xi2;
        Tau22res[l] = tau22;
        EtaAll.col(l) = Eta;
        l++;
      }
      
    }
    
  }
  
  for(int i = 0; i< (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
  
  List output = List::create(
    // Rcpp::Named("beta.hat") = bhat,
    // Rcpp::Named("beta.se") = se,
    // Rcpp::Named("beta.p.value") = pvalue,
    // Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    // Rcpp::Named("WRes") = Rcpp::wrap(obj.WRes),
    // Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("EtaIterRate") = EtaIterRate,
    Rcpp::Named("Beta1res") = Beta1res, // change notation.
    Rcpp::Named("Beta2res") = Beta2res, // change notation.
    Rcpp::Named("Sgga2res") = Sgga2res,
    Rcpp::Named("Tau12res") = Tau12res, // change notation.
    Rcpp::Named("Tau22res") = Tau22res  // change notation.
  
  );
  return output;
  
  
  
  
}

ObjMVMRCUE MVMRCUEObj(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                      arma::vec &se2,
                      Options_MVMRCUE* opts){
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  vec agM = opts -> agM;
  vec bgM = opts -> bgM;
  double atau1 = opts -> atau1;
  double btau1 = opts -> btau1;
  double atau2 = opts -> atau2;
  double btau2 = opts -> btau2;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  // ----------------------------------------------------------------------
  // initial values
  int J = gammahM.n_cols;
  int p = gammahM.n_rows;
  ivec Eta = zeros<ivec>(p, 1);
  vec Beta1 = 0.01*ones(J, 1);
  vec Beta2 = 0.01*ones(J, 1);
  vec Sgga2 = 0.01*ones(J, 1);
  double tau12 = 0.01;
  double tau22 = 0.01;
  double w = 0.1;
  double xi2 = 0.01;
  double tau12xi2 = tau12*xi2;
  
  mat Mu = 0.01*ones(p, J);
  vec mut = 0.01*ones(p, 1);
  mat Bgamma = zeros(p, J);
  mat Bgammat = zeros(p, J);
  
  
  
  int numsave = maxIter / thin;
  mat Beta1res = ones(numsave, J);
  mat Beta2res = ones(numsave, J);
  mat Sgga2res = ones(numsave, J);
  vec Tau12res = ones(numsave, 1);
  vec Tau22res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  vec EtaIterRate = zeros(p, 1);
  
  vec GinvsG2 = Gammah / (se2 % se2);
  mat ginvsg2M = zeros(p, J);
  
  for(int j = 0; j < J; j++){
    ginvsg2M.col(j) = gammahM.col(j) / (se1M.col(j) % se1M.col(j));
  }
  mat Invsg2 = 1. / (se1M % se1M);
  vec invsG2 = 1. / (se2 % se2);
  
  int l = 0;
  // cout << "MU: " << Mu.row(0) << endl;
  
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    
    vec Invsgga2 = 1. / Sgga2;
    double invtau12xi2 = 1. / (tau12xi2);
    double invtau22 = 1. / tau22;
    
    for(int j = 0; j < J; j++){
      // cout << "j: "<< j << "Beta1: " <<Beta1[j] <<"Mu: "<<  sum(Mu.col(j)) << endl;
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    // 
    // cout << "Bgamma: sum: "<< sum(Bgamma, 0) << endl;
    // cout << "Bgammat: sum: "<<  sum(Bgammat, 0) << endl;
    // --------------------------------------- //
    // Update Gamma
    // --------------------------------------- //
    vec v20t = 1. / (invsG2 + invtau12xi2);
    vec v21t = 1. / (invsG2 + invtau22);
    
    vec mutm0t = (GinvsG2 + sum(Bgamma, 1) * invtau12xi2) % v20t;
    vec mutm1t = (GinvsG2 + sum(Bgammat, 1) * invtau22) % v21t;
    
    // cout << "mutm0t: "<< mutm0t.subvec(0, 2).t() << " sum: " << sum(mutm0t) << endl;
    // cout << "mutm1t: "<< mutm1t.subvec(0, 2).t() << " sum: " << sum(mutm1t) << endl;
    
    for(int k = 0; k < p; ++k) {
      if(Eta[k] == 1) {
        // mut[k] = mutm1t[k];
        mut[k] = mutm1t[k] + randn()*sqrt(v21t[k]);
      } else if(Eta[k] == 0) {
        // mut[k] = mutm0t[k];
        mut[k] = mutm0t[k] + randn()*sqrt(v20t[k]);
      }
    }
    
    // cout << "mut: " << sum(mut) << endl;
    // --------------------------------------- //
    // Update gamma for each j in 1,...,J
    // --------------------------------------- //
    mat V20M(p, J, fill::none);
    mat V21M(p, J, fill::none);
    mat Mutm0(p, J, fill::none);
    mat Mutm1(p, J, fill::none);
    
    for(int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      V20M.col(j) = 1.0 / (Invsgga2[j] + Invsg2.col(j) + Beta1[j] * Beta1[j] * invtau12xi2);
      V21M.col(j) = 1.0 / (Invsgga2[j] + Invsg2.col(j) + Beta2[j] * Beta2[j] * invtau22);
      vec tmp1 = Gminsb1g + Beta1[j] * Mu.col(j);
      vec tmp2 = Gminsb2g + Beta2[j] * Mu.col(j);
      Mutm0.col(j) = (ginvsg2M.col(j) + tmp1 * Beta1[j] * invtau12xi2) % V20M.col(j);
      Mutm1.col(j) = (ginvsg2M.col(j) + tmp2 * Beta2[j] * invtau22) % V21M.col(j);
      
      for(int k = 0; k < p; ++k) {
        if(Eta[k] == 1) {
          // Mu(k, j) = Mutm1(k, j);
          Mu(k, j) = Mutm1(k, j) + randn()*sqrt(V21M(k, j));
        } else {
          // Mu(k, j) = Mutm0(k, j);
          Mu(k, j) = Mutm0(k, j) + randn()*sqrt(V20M(k, j));
        }
      }
      
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    
    // cout << "Mu: "<< sum(Mu, 0) << endl; 
    // --------------------------------------- //
    // Update Beta1 and Beta2
    // --------------------------------------- //
    for(int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec muEta0j = (1 - Eta) % Mu.col(j);
      vec muEta1j = Eta % Mu.col(j);
      vec tmpmuj = Gminsb1g + Beta1[j] * Mu.col(j);
      vec tmpmujt = Gminsb2g + Beta2[j] * Mu.col(j);
      
      if(sum(Eta) == p) {
        Beta1[j] = 0;
      } else {
        double sig2b0j = 1.0 / (sum(muEta0j % Mu.col(j)) * invtau12xi2);
        double mub0j = sum(muEta0j % tmpmuj) * invtau12xi2 * sig2b0j;
        // Beta1[j] = mub0j;
        Beta1[j] = mub0j + randn()*sqrt(sig2b0j);
      }
      
      if(sum(Eta) == 0) {
        Beta2[j] = 0;
      } else {
        double sig2b1j = 1.0 / (sum(muEta1j % Mu.col(j)) * invtau22);
        double mub1j = sum(muEta1j % tmpmujt) * invtau22 * sig2b1j;
        // Beta2[j] = mub1j;
        Beta2[j] = mub1j + randn()*sqrt(sig2b1j);
      }
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    // cout << "Beta1: "<< Beta1.t() << endl;
    // cout << "Beta2: " << Beta2.t() << endl;
    // --------------------------------------- //
    // Update the variance terms
    // --------------------------------------- //
    for(int j = 0; j < J; ++j) {
      double tagmj = agM[j] + p / 2.0;
      double tbgmj = bgM[j] + sum(Mu.col(j) % Mu.col(j)) / 2.0;
      Sgga2[j] = 1.0 / randg<double>(distr_param(tagmj, 1.0 / tbgmj));
      // Sgga2[j] = tbgmj;
    }
    
    vec Gminsb1g = mut - sum(Bgamma, 1);
    vec Gminsb2g = mut - sum(Bgammat, 1);
    vec err02 = (1 - Eta) % Gminsb1g % Gminsb1g;
    vec err12 = Eta % Gminsb2g % Gminsb2g;
    
    double tatau1 = atau1 + sum(1 - Eta) / 2.0;
    double tbtau1 = btau1 + sum(err02) / (2.0 * xi2);
    tau12 = 1.0 / randg<double>(distr_param(tatau1, 1.0 / tbtau1));
    // tau12= tbtau1;
    
    
    double tatau2 = atau2 + sum(Eta) / 2.0;
    double tbtau2 = btau2 + sum(err12) / 2.0;
    tau22 = 1.0 / randg<double>(distr_param(tatau2, 1.0 / tbtau2));
    // tau22 = tbtau2;
    
    double taxi2 = 0.5 * sum(1 - Eta);
    double tbxi2 = 0.5 * sum(err02) / tau12;
    xi2 = 1.0 / randg<double>(distr_param(taxi2, 1.0 / tbxi2));
    // xi2 = tbxi2;
    
    tau12xi2 = tau12*xi2;
    if(tau12xi2 < 1e-7){
      tau12xi2 = 1e-7;
    }
    
    // cout << "Sgga2: "<< Sgga2.t() << endl;
    // cout << "tau12: "<< tau12<< endl;
    // cout << "tau22: "<< tau22 << endl;
    // cout <<"xi2: "<< xi2 << endl;
    // --------------------------------------- //
    // Update omega
    // --------------------------------------- //
    double wpa1 = a + sum(Eta);
    double wpa2 = b + p - sum(Eta);
    w = R::rbeta(wpa1, wpa2);
    
    // --------------------------------------- //
    // Update Eta
    // --------------------------------------- //
    double tau1 = sqrt(tau12 * xi2);
    double tau2 = sqrt(tau22);
    vec bmu0 = sum(Bgamma, 1);
    vec bmu1 = sum(Bgammat, 1);
    // vec tmp = zeros(p, 1);
    for(int m = 0; m < p; ++m) {
      double p1 = w*normal_pdf(mut[m], bmu1[m], tau2);
      double p0 = (1 - w) * normal_pdf(mut[m], bmu0[m], tau1);
      double prob = p1 / (p1 + p0);
      Eta[m] = R::rbinom(1, prob);
      // tmp[m] = prob;
    }
    
    // cout << "tmp: "<< tmp.subvec(0, 4).t() << " sum: "<< sum(tmp) << endl;
    // cout << "Eta: " << sum(Eta) << endl;
    // 
    // cout << "=========================================="<<endl;
    
    // Store results
    if(iter >= burnin) {
      if((iter - burnin) % thin == 0) {
        Sgga2res.row(l) = Sgga2.t();
        Beta1res.row(l) = Beta1.t();
        Beta2res.row(l) = Beta2.t();
        Tau12res[l] = tau12xi2;
        Tau22res[l] = tau22;
        EtaAll.col(l) = Eta;
        l++;
      }
      
    }
    
  }
  
  // cout << "check error 00" << endl;
  for(int i = 0; i< (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
  
  ObjMVMRCUE obj;
  
  obj.EtaIterRate = EtaIterRate;
  obj.Beta1res = Beta1res;
  obj.Beta2res = Beta2res;
  obj.Sgga2res = Sgga2res;
  obj.Tau12res = Tau12res;
  obj.Tau22res = Tau22res;
  
  return obj;
  
}


// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List MVMRCUEoverlapfun(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                arma::vec &se2, arma::mat &Re, arma::ivec &Eta){
  
  int J = gammahM.n_cols;
  int p = gammahM.n_rows;
  
  // ---------------------------------------- //
  // construct the Delta matrix
  field<mat> Delta_list(p);
  // double sum = 0;
  for (int k = 0; k < p; ++k) {
    vec tmpsigma = join_vert(se1M.row(k).t(), se2.subvec(k, k));
    mat Sig_k = kron(tmpsigma, tmpsigma.t());
    mat Re_Sig_k = Re % Sig_k;
    Delta_list(k) = inv(Re_Sig_k);
    // sum += accu(Delta_list(k));
  }
  // cout<<"Delta_list::"<<sum << endl;
  // ---------------------------------------- //
  
  // hyper-prior parameters.
  vec agM = 0.001*ones(J, 1);
  vec bgM = 0.001*ones(J, 1);
  double atau1 = 0.001;
  double btau1 = 0.001;
  double atau2 = 0.001;
  double btau2 = 0.001;
  double a = 1;
  double b = 1;
  int maxIter = 4000; int burnin = 1000; int thin = 10;
  // int maxIter = 4; int burnin = 0; int thin = 1;
  // initial values.
  // ivec Eta = zeros<ivec>(p, 1);
  vec Beta1 = 0.01*ones(J, 1);
  vec Beta2 = 0.01*ones(J, 1);
  vec Sgga2 = 0.01*ones(J, 1);

  double tau12 = 0.01;
  double tau22 = 0.01;
  double w = 0.1;
  double xi2 = 0.01;
  double tau12xi2 = tau12*xi2;

  mat Mu = 0.01*ones(p, J);
  vec mut = 0.01*ones(p, 1);
  mat Bgamma = zeros(p, J);
  mat Bgammat = zeros(p, J);



  int numsave = maxIter / thin;
  mat Beta1res = ones(numsave, J);
  mat Beta2res = ones(numsave, J);
  mat Sgga2res = ones(numsave, J);
  vec Tau12res = ones(numsave, 1);
  vec Tau22res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  vec EtaIterRate = zeros(p, 1);
  
  
  int l = 0;

  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    vec Invsgga2 = 1. / Sgga2;
    double invtau12xi2 = 1. / (tau12 * xi2);
    double invtau22 = 1. / tau22;

    for(int j = 0; j < J; j++){
      // cout << "j: "<< j << "Beta1: " <<Beta1[j] <<"Mu: "<<  sum(Mu.col(j)) << endl;
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }

    // --------------------------------------- //
    // Update Gamma
    // --------------------------------------- //
    vec bgItau1 = sum(Bgamma, 1) * invtau12xi2;
    vec bgItau2 = sum(Bgammat, 1) * invtau22;
    mat Gmmu = gammahM - Mu;

    // mat Delta_k;
    // double deltaTT, GinvsG2_k, degmmu_k, v20t_k, v21t_k, mutm0t_k, mutm1t_k;

    for (int k = 0; k < p; ++k) {
      mat Delta_k = Delta_list(k);
      double deltaTT = Delta_k(J, J);
      double GinvsG2_k = Gammah(k) * deltaTT;
      double degmmu_k = sum(Delta_k.submat(0, J, J-1, J) % Gmmu.row(k).t());
      
      double v20t_k = 1.0 / (deltaTT + invtau12xi2);
      double v21t_k = 1.0 / (deltaTT + invtau22);
      double mutm0t_k = (GinvsG2_k + degmmu_k + bgItau1(k)) * v20t_k;
      double mutm1t_k = (GinvsG2_k + degmmu_k + bgItau2(k)) * v21t_k;
      
      if (Eta(k) == 1) {
        mut[k] = mutm1t_k + randn()*sqrt(v21t_k);
        // mut[k] = mutm1t_k; mutm1t[k] + randn()*sqrt(v21t[k]);
        // if(k < 5){
        //   cout << "k:" << k << "__v20t_k:__" << v21t_k << endl;
        // }
      } else {
        mut[k] = mutm0t_k + randn()*sqrt(v20t_k);
        // mut[k] = mutm0t_k;
        // if(k < 5){
        //   cout << "k:" << k << "v21t_k:" << v20t_k << endl;
        // }
      }
   
    }
    
    
    vec Gmmu1 = Gammah - mut;
    for (int j = 0; j < J; ++j) {
      Gmmu = gammahM - Mu;
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec tmp1 = (Gminsb1g + Beta1(j) * Mu.col(j)) * Beta1(j) * invtau12xi2;
      vec tmp2 = (Gminsb2g + Beta2(j) * Mu.col(j)) * Beta2(j) * invtau22;
      
      for (int k = 0; k < p; ++k) {
        mat Delta_k = Delta_list(k);
        double delta_jj = Delta_k(j, j);
        double delta_jJ1 = Delta_k(j, J);
        
        double term3 = sum(Delta_k.row(j).subvec(0, J-1) % Gmmu.row(k)) - delta_jj * Gmmu(k, j);
        
        double V20M_kj = 1.0 / (Invsgga2(j) + delta_jj + Beta1(j) * Beta1(j) * invtau12xi2);
        double V21M_kj = 1.0 / (Invsgga2(j) + delta_jj + Beta2(j) * Beta2(j) * invtau22);
        
        double Mutm0_kj = (delta_jj * gammahM(k, j) + delta_jJ1 * Gmmu1(k) + term3 + tmp1(k)) * V20M_kj;
        double Mutm1_kj = (delta_jj * gammahM(k, j) + delta_jJ1 * Gmmu1(k) + term3 + tmp2(k)) * V21M_kj;
        
        if (Eta(k) == 1) {
          Mu(k, j) = Mutm1_kj + randn()*sqrt(V21M_kj);
          // Mu(k, j) = Mutm1_kj;
          // if(k <4){
          //   cout <<"sd:j_" <<j << "__k:"<< k <<"__Mutm1_kj:" << Mutm1_kj <<"__sd::"<<V21M_kj<< endl;
          // }
        } else {
          Mu(k, j) = Mutm0_kj + randn()*sqrt(V20M_kj);
          // Mu(k, j) = Mutm0_kj;
          // if(k <4){
          //   cout <<"sd:j_" <<j << "__k:"<< k <<"__Mutm0_kj:" << Mutm0_kj <<"__sd::"<<V20M_kj<< endl;
          // }
        }
      }
      
      Bgamma.col(j) = Beta1(j) * Mu.col(j);
      Bgammat.col(j) = Beta2(j) * Mu.col(j);
    }
    
    
    
    for (int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec muEta0j = (1 - Eta) % Mu.col(j);
      vec muEta1j = Eta % Mu.col(j);
      vec tmpmuj = Gminsb1g + Beta1(j) * Mu.col(j);
      vec tmpmujt = Gminsb2g + Beta2(j) * Mu.col(j);
      
      if (sum(Eta) == p) {
        Beta1(j) = 0;
      } else {
        double sig2b0j = 1.0 / (sum(muEta0j % Mu.col(j)) * invtau12xi2);
        double mub0j = sum(muEta0j % tmpmuj) * invtau12xi2 * sig2b0j;
        Beta1[j] = mub0j + randn() * sqrt(sig2b0j);
        // Beta1[j] = mub0j;
        // cout << "j:"<<j << "__mub0j:"<< mub0j<<"__sig2b0j"<<sig2b0j<<endl;
      }
      
      if (sum(Eta) == 0) {
        Beta2(j) = 0;
      } else {
        double sig2b1j = 1.0 / (sum(muEta1j % Mu.col(j)) * invtau22);
        double mub1j = sum(muEta1j % tmpmujt) * invtau22 * sig2b1j;
        Beta2[j] = mub1j + randn() * sqrt(sig2b1j);
        // Beta2[j] = mub1j;
        // cout << "j:"<<j << "__mub1j:"<< mub1j<<"__sig2b1j"<<sig2b1j<<endl;
      }
      
      Bgamma.col(j) = Beta1(j) * Mu.col(j);
      Bgammat.col(j) = Beta2(j) * Mu.col(j);
    }
    
    
    for (int j = 0; j < J; ++j) {
      double tagmj = agM(j) + p / 2.0;
      double tbgmj = bgM(j) + sum(Mu.col(j) % Mu.col(j)) / 2.0;
      Sgga2[j] = 1.0 / randg<double>(distr_param(tagmj, 1.0 / tbgmj));
      // Sgga2[j] = tbgmj;
      // cout << "j:" << j << "__tagmj::" << tagmj << "__tbgmj::" << tbgmj <<endl;
    }
    
    vec Gminsb1g = mut - sum(Bgamma, 1);
    vec Gminsb2g = mut - sum(Bgammat, 1);
    vec err02 = (1 - Eta) % Gminsb1g % Gminsb1g;
    vec err12 = Eta % Gminsb2g % Gminsb2g;
    
    double tatau1 = atau1 + sum(1 - Eta) / 2.0;
    double tbtau1 = btau1 + sum(err02) / (2.0 * xi2);
    tau12 = 1.0 / randg<double>(distr_param(tatau1, 1.0 / tbtau1));
    // tau12 = tatau1;
    
    double tatau2 = atau2 + sum(Eta) / 2.0;
    double tbtau2 = btau2 + sum(err12) / 2.0;
    tau22 = 1.0 / randg<double>(distr_param(tatau2, 1.0 / tbtau2));
    // tau22 = tatau2;
    
    double taxi2 = 0.5 * sum(1 - Eta);
    double tbxi2 = 0.5 * sum(err02) / tau12;
    xi2 = 1.0 / randg<double>(distr_param(taxi2, 1.0 / tbxi2));
    // xi2 = taxi2;
    
    // cout <<"tatau1::"<<tatau1<<"__tbtau1::"<<tbtau1<< "__tatau2::"<<tatau2<<"__tbtau2::"<<tbtau2<<"__taxi2::"<<taxi2<<"__tbxi2::"<<tbxi2<<endl;
    
    // update omega
    double wpa1 = a + sum(Eta);
    double wpa2 = b + p - sum(Eta);
    w = R::rbeta(wpa1, wpa2);
    // w = 0.1;

    // update Eta
    double tau1 = sqrt(tau12 * xi2);
    double tau2 = sqrt(tau22);
    vec bmu0 = sum(Bgamma, 1);
    vec bmu1 = sum(Bgammat, 1);
    for (int m = 0; m < p; ++m) {
      double p1 = w*normal_pdf(mut[m], bmu1[m], tau2);
      double p0 = (1 - w) * normal_pdf(mut[m], bmu0[m], tau1);
      double prob = p1 / (p1 + p0);
      Eta[m] = R::rbinom(1, prob);
    }
    // cout << "mut:" << mean(mut) << "_MU::" << mean(Mu, 0) << endl;
    // cout <<"-----------------------"<<endl;
    // Store results
    if(iter >= burnin) {
      if((iter - burnin) % thin == 0) {
        Sgga2res.row(l) = Sgga2.t();
        Beta1res.row(l) = Beta1.t();
        Beta2res.row(l) = Beta2.t();
        Tau12res[l] = tau12xi2;
        Tau22res[l] = tau22;
        EtaAll.col(l) = Eta;
        l++;
      }
      
    }



  }
  
  for(int i = 0; i< (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
  List output = List::create(
    // Rcpp::Named("beta.hat") = bhat,
    // Rcpp::Named("beta.se") = se,
    // Rcpp::Named("beta.p.value") = pvalue,
    // Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    // Rcpp::Named("WRes") = Rcpp::wrap(obj.WRes),
    // Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("EtaIterRate") = EtaIterRate,
    Rcpp::Named("Beta1res") = Beta1res, // change notation.
    Rcpp::Named("Beta2res") = Beta2res, // change notation.
    Rcpp::Named("Sgga2res") = Sgga2res,
    Rcpp::Named("Tau12res") = Tau12res, // change notation.
    Rcpp::Named("Tau22res") = Tau22res  // change notation.
  
  );
  return output;
  
}


ObjMVMRCUE MVMRCUE1Obj(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                      arma::vec &se2, arma::mat &Re,
                      Options_MVMRCUE* opts){
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  vec agM = opts -> agM;
  vec bgM = opts -> bgM;
  double atau1 = opts -> atau1;
  double btau1 = opts -> btau1;
  double atau2 = opts -> atau2;
  double btau2 = opts -> btau2;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  // ----------------------------------------------------------------------
  int J = gammahM.n_cols;
  int p = gammahM.n_rows;
  
  // ---------------------------------------- //
  // construct the Delta matrix
  field<mat> Delta_list(p);
  // double sum = 0;
  for (int k = 0; k < p; ++k) {
    vec tmpsigma = join_vert(se1M.row(k).t(), se2.subvec(k, k));
    mat Sig_k = kron(tmpsigma, tmpsigma.t());
    mat Re_Sig_k = Re % Sig_k;
    Delta_list(k) = inv(Re_Sig_k);
    // sum += accu(Delta_list(k));
  }
  // cout<<"Delta_list::"<<sum << endl;
  // ---------------------------------------- //
  
  // hyper-prior parameters.
  // vec agM = 0.001*ones(J, 1);
  // vec bgM = 0.001*ones(J, 1);
  // double atau1 = 0.001;
  // double btau1 = 0.001;
  // double atau2 = 0.001;
  // double btau2 = 0.001;
  // double a = 1;
  // double b = 1;
  // int maxIter = 4000; int burnin = 1000; int thin = 10;
  // int maxIter = 4; int burnin = 0; int thin = 1;
  // initial values.
  ivec Eta = zeros<ivec>(p, 1);
  vec Beta1 = 0.01*ones(J, 1);
  vec Beta2 = 0.01*ones(J, 1);
  vec Sgga2 = 0.01*ones(J, 1);
  
  double tau12 = 0.01;
  double tau22 = 0.01;
  double w = 0.1;
  double xi2 = 0.01;
  double tau12xi2 = tau12*xi2;
  
  mat Mu = 0.01*ones(p, J);
  vec mut = 0.01*ones(p, 1);
  mat Bgamma = zeros(p, J);
  mat Bgammat = zeros(p, J);
  
  int numsave = maxIter / thin;
  mat Beta1res = ones(numsave, J);
  mat Beta2res = ones(numsave, J);
  mat Sgga2res = ones(numsave, J);
  vec Tau12res = ones(numsave, 1);
  vec Tau22res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  vec EtaIterRate = zeros(p, 1);
  
  
  int l = 0;
  
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    vec Invsgga2 = 1. / Sgga2;
    double invtau12xi2 = 1. / (tau12 * xi2);
    double invtau22 = 1. / tau22;
    
    for(int j = 0; j < J; j++){
      // cout << "j: "<< j << "Beta1: " <<Beta1[j] <<"Mu: "<<  sum(Mu.col(j)) << endl;
      Bgamma.col(j) = Beta1[j] * Mu.col(j);
      Bgammat.col(j) = Beta2[j] * Mu.col(j);
    }
    
    // --------------------------------------- //
    // Update Gamma
    // --------------------------------------- //
    vec bgItau1 = sum(Bgamma, 1) * invtau12xi2;
    vec bgItau2 = sum(Bgammat, 1) * invtau22;
    mat Gmmu = gammahM - Mu;
    
    // mat Delta_k;
    // double deltaTT, GinvsG2_k, degmmu_k, v20t_k, v21t_k, mutm0t_k, mutm1t_k;
    
    for (int k = 0; k < p; ++k) {
      mat Delta_k = Delta_list(k);
      double deltaTT = Delta_k(J, J);
      double GinvsG2_k = Gammah(k) * deltaTT;
      double degmmu_k = sum(Delta_k.submat(0, J, J-1, J) % Gmmu.row(k).t());
      
      double v20t_k = 1.0 / (deltaTT + invtau12xi2);
      double v21t_k = 1.0 / (deltaTT + invtau22);
      double mutm0t_k = (GinvsG2_k + degmmu_k + bgItau1(k)) * v20t_k;
      double mutm1t_k = (GinvsG2_k + degmmu_k + bgItau2(k)) * v21t_k;
      
      if (Eta(k) == 1) {
        mut[k] = mutm1t_k + randn()*sqrt(v21t_k);
        // mut[k] = mutm1t_k; mutm1t[k] + randn()*sqrt(v21t[k]);
        // if(k < 5){
        //   cout << "k:" << k << "__v20t_k:__" << v21t_k << endl;
        // }
      } else {
        mut[k] = mutm0t_k + randn()*sqrt(v20t_k);
        // mut[k] = mutm0t_k;
        // if(k < 5){
        //   cout << "k:" << k << "v21t_k:" << v20t_k << endl;
        // }
      }
      
    }
    
    
    vec Gmmu1 = Gammah - mut;
    for (int j = 0; j < J; ++j) {
      Gmmu = gammahM - Mu;
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec tmp1 = (Gminsb1g + Beta1(j) * Mu.col(j)) * Beta1(j) * invtau12xi2;
      vec tmp2 = (Gminsb2g + Beta2(j) * Mu.col(j)) * Beta2(j) * invtau22;
      
      for (int k = 0; k < p; ++k) {
        mat Delta_k = Delta_list(k);
        double delta_jj = Delta_k(j, j);
        double delta_jJ1 = Delta_k(j, J);
        
        double term3 = sum(Delta_k.row(j).subvec(0, J-1) % Gmmu.row(k)) - delta_jj * Gmmu(k, j);
        
        double V20M_kj = 1.0 / (Invsgga2(j) + delta_jj + Beta1(j) * Beta1(j) * invtau12xi2);
        double V21M_kj = 1.0 / (Invsgga2(j) + delta_jj + Beta2(j) * Beta2(j) * invtau22);
        
        double Mutm0_kj = (delta_jj * gammahM(k, j) + delta_jJ1 * Gmmu1(k) + term3 + tmp1(k)) * V20M_kj;
        double Mutm1_kj = (delta_jj * gammahM(k, j) + delta_jJ1 * Gmmu1(k) + term3 + tmp2(k)) * V21M_kj;
        
        if (Eta(k) == 1) {
          Mu(k, j) = Mutm1_kj + randn()*sqrt(V21M_kj);
          // Mu(k, j) = Mutm1_kj;
          // if(k <4){
          //   cout <<"sd:j_" <<j << "__k:"<< k <<"__Mutm1_kj:" << Mutm1_kj <<"__sd::"<<V21M_kj<< endl;
          // }
        } else {
          Mu(k, j) = Mutm0_kj + randn()*sqrt(V20M_kj);
          // Mu(k, j) = Mutm0_kj;
          // if(k <4){
          //   cout <<"sd:j_" <<j << "__k:"<< k <<"__Mutm0_kj:" << Mutm0_kj <<"__sd::"<<V20M_kj<< endl;
          // }
        }
      }
      
      Bgamma.col(j) = Beta1(j) * Mu.col(j);
      Bgammat.col(j) = Beta2(j) * Mu.col(j);
    }
    
    
    
    for (int j = 0; j < J; ++j) {
      vec Gminsb1g = mut - sum(Bgamma, 1);
      vec Gminsb2g = mut - sum(Bgammat, 1);
      vec muEta0j = (1 - Eta) % Mu.col(j);
      vec muEta1j = Eta % Mu.col(j);
      vec tmpmuj = Gminsb1g + Beta1(j) * Mu.col(j);
      vec tmpmujt = Gminsb2g + Beta2(j) * Mu.col(j);
      
      if (sum(Eta) == p) {
        Beta1(j) = 0;
      } else {
        double sig2b0j = 1.0 / (sum(muEta0j % Mu.col(j)) * invtau12xi2);
        double mub0j = sum(muEta0j % tmpmuj) * invtau12xi2 * sig2b0j;
        Beta1[j] = mub0j + randn() * sqrt(sig2b0j);
        // Beta1[j] = mub0j;
        // cout << "j:"<<j << "__mub0j:"<< mub0j<<"__sig2b0j"<<sig2b0j<<endl;
      }
      
      if (sum(Eta) == 0) {
        Beta2(j) = 0;
      } else {
        double sig2b1j = 1.0 / (sum(muEta1j % Mu.col(j)) * invtau22);
        double mub1j = sum(muEta1j % tmpmujt) * invtau22 * sig2b1j;
        Beta2[j] = mub1j + randn() * sqrt(sig2b1j);
        // Beta2[j] = mub1j;
        // cout << "j:"<<j << "__mub1j:"<< mub1j<<"__sig2b1j"<<sig2b1j<<endl;
      }
      
      Bgamma.col(j) = Beta1(j) * Mu.col(j);
      Bgammat.col(j) = Beta2(j) * Mu.col(j);
    }
    
    
    for (int j = 0; j < J; ++j) {
      double tagmj = agM(j) + p / 2.0;
      double tbgmj = bgM(j) + sum(Mu.col(j) % Mu.col(j)) / 2.0;
      Sgga2[j] = 1.0 / randg<double>(distr_param(tagmj, 1.0 / tbgmj));
      // Sgga2[j] = tbgmj;
      // cout << "j:" << j << "__tagmj::" << tagmj << "__tbgmj::" << tbgmj <<endl;
    }
    
    vec Gminsb1g = mut - sum(Bgamma, 1);
    vec Gminsb2g = mut - sum(Bgammat, 1);
    vec err02 = (1 - Eta) % Gminsb1g % Gminsb1g;
    vec err12 = Eta % Gminsb2g % Gminsb2g;
    
    double tatau1 = atau1 + sum(1 - Eta) / 2.0;
    double tbtau1 = btau1 + sum(err02) / (2.0 * xi2);
    tau12 = 1.0 / randg<double>(distr_param(tatau1, 1.0 / tbtau1));
    // tau12 = tatau1;
    
    double tatau2 = atau2 + sum(Eta) / 2.0;
    double tbtau2 = btau2 + sum(err12) / 2.0;
    tau22 = 1.0 / randg<double>(distr_param(tatau2, 1.0 / tbtau2));
    // tau22 = tatau2;
    
    double taxi2 = 0.5 * sum(1 - Eta);
    double tbxi2 = 0.5 * sum(err02) / tau12;
    xi2 = 1.0 / randg<double>(distr_param(taxi2, 1.0 / tbxi2));
    // xi2 = taxi2;
    
    // cout <<"tatau1::"<<tatau1<<"__tbtau1::"<<tbtau1<< "__tatau2::"<<tatau2<<"__tbtau2::"<<tbtau2<<"__taxi2::"<<taxi2<<"__tbxi2::"<<tbxi2<<endl;
    
    // update omega
    double wpa1 = a + sum(Eta);
    double wpa2 = b + p - sum(Eta);
    w = R::rbeta(wpa1, wpa2);
    // w = 0.1;
    
    // update Eta
    double tau1 = sqrt(tau12 * xi2);
    double tau2 = sqrt(tau22);
    vec bmu0 = sum(Bgamma, 1);
    vec bmu1 = sum(Bgammat, 1);
    for (int m = 0; m < p; ++m) {
      double p1 = w*normal_pdf(mut[m], bmu1[m], tau2);
      double p0 = (1 - w) * normal_pdf(mut[m], bmu0[m], tau1);
      double prob = p1 / (p1 + p0);
      Eta[m] = R::rbinom(1, prob);
    }
    // cout << "mut:" << mean(mut) << "_MU::" << mean(Mu, 0) << endl;
    // cout <<"-----------------------"<<endl;
    // Store results
    if(iter >= burnin) {
      if((iter - burnin) % thin == 0) {
        Sgga2res.row(l) = Sgga2.t();
        Beta1res.row(l) = Beta1.t();
        Beta2res.row(l) = Beta2.t();
        Tau12res[l] = tau12xi2;
        Tau22res[l] = tau22;
        EtaAll.col(l) = Eta;
        l++;
      }
      
    }
    
    
    
  }
  
  for(int i = 0; i< (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
  
  ObjMVMRCUE obj;
  
  obj.EtaIterRate = EtaIterRate;
  obj.Beta1res = Beta1res;
  obj.Beta2res = Beta2res;
  obj.Sgga2res = Sgga2res;
  obj.Tau12res = Tau12res;
  obj.Tau22res = Tau22res;
  
  return obj;
  
}


// [[Rcpp::export]]
Rcpp::List MVMRCUEIndepSample(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                   arma::vec &se2, SEXP opts = R_NilValue){
  
  
  
  
  
  uword J = gammahM.n_cols;
  Options_MVMRCUE* lp_opt = NULL;
  
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_MVMRCUE(opt["agM"], opt["bgM"], opt["atau1"], opt["btau1"],
                                 opt["atau2"], opt["btau2"], opt["a"], opt["b"],
                                     opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_MVMRCUE(J);
  }
  
  ObjMVMRCUE obj = MVMRCUEObj(gammahM, Gammah, se1M, se2,lp_opt);
  
  
  vec bhats(J);
  vec ses(J);
  vec pvalues(J);
  
  for(int j = 0; j < J; ++j) {
    vec column = obj.Beta1res.col(j);
    double bhat = mean(column);
    double se = stddev(column);
    
    // Calculate p-value using the standard normal distribution
    double pvalue = 2 * (R::pnorm(abs(bhat / se), 0, 1, 0, 0));
    
    bhats[j] = bhat;
    ses[j] = se;
    pvalues[j] = pvalue;
  }
  
  List output = List::create(
    Rcpp::Named("beta.hats") = bhats,
    Rcpp::Named("beta.ses") = ses,
    Rcpp::Named("beta.p.values") = pvalues,
    Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta2res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2res),
    Rcpp::Named("Tau12Res") = Rcpp::wrap(obj.Tau12res),
    Rcpp::Named("Tau22Res") = Rcpp::wrap(obj.Tau22res)
    
  );
  return output;
  
}

// [[Rcpp::export]]
Rcpp::List MVMRCUE(arma::mat &gammahM, arma::vec &Gammah, arma::mat &se1M,
                   arma::vec &se2, arma::mat &Re, SEXP opts = R_NilValue){
  
  uword J = gammahM.n_cols;
  Options_MVMRCUE* lp_opt = NULL;
  
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_MVMRCUE(opt["agM"], opt["bgM"], opt["atau1"], opt["btau1"],
                                 opt["atau2"], opt["btau2"], opt["a"], opt["b"],
                                     opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_MVMRCUE(J);
  }
  
  ObjMVMRCUE obj = MVMRCUE1Obj(gammahM, Gammah, se1M, se2, Re, lp_opt);
  
  
  vec bhats(J);
  vec ses(J);
  vec pvalues(J);
  
  for(int j = 0; j < J; ++j) {
    vec column = obj.Beta1res.col(j);
    double bhat = mean(column);
    double se = stddev(column);
    
    // Calculate p-value using the standard normal distribution
    double pvalue = 2 * (R::pnorm(abs(bhat / se), 0, 1, 0, 0));
    
    bhats[j] = bhat;
    ses[j] = se;
    pvalues[j] = pvalue;
  }
  
  List output = List::create(
    Rcpp::Named("beta.hats") = bhats,
    Rcpp::Named("beta.ses") = ses,
    Rcpp::Named("beta.p.values") = pvalues,
    Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta2res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2res),
    Rcpp::Named("Tau12Res") = Rcpp::wrap(obj.Tau12res),
    Rcpp::Named("Tau22Res") = Rcpp::wrap(obj.Tau22res)
    
  );
  return output;
  
}
