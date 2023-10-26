#include <RcppArmadillo.h>
#include <cmath>        // std::abs
#include <tuple>
#include <iostream>


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double error(arma::colvec y, arma::colvec mu) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += mu[i];
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) - y[i] + mu[i];
    }
  }
  return sum;
}

arma::mat UpdateExp(arma::mat data, arma::mat exposures, arma::mat signatures){

}

// [[Rcpp::export]]
List nmfgen(arma::mat data, int noSignatures, int iter = 5000) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;

  arma::mat exposures(genomes, noSignatures,arma::fill::randu);
  
  arma::mat signatures(noSignatures, mutTypes, arma::fill::randu);
  
  arma::mat estimate = exposures * signatures;
  arma::mat fraq = data/estimate;
  
  for(int t = 0; t < iter; t++){

    signatures = signatures % (arma::trans(exposures) * fraq);
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );

    estimate = exposures * signatures;
    fraq = data/estimate;

    exposures = exposures % (fraq * arma::trans(signatures));

    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    fraq = data/estimate;

  }
  
  double gkl = error(arma::vectorise(data),arma::vectorise(estimate));
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gkl);
  return output;
}


std::tuple<arma::mat, arma::mat, double> nmf1(arma::mat data, int noSignatures, int iter = 5000) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;

  arma::mat exposures(genomes, noSignatures,arma::fill::randu);
  
  arma::mat signatures(noSignatures, mutTypes, arma::fill::randu);
  
  arma::mat estimate = exposures * signatures;
  arma::mat fraq = data/estimate;
  
  
  for(int t = 0; t < iter; t++){
    
    signatures = signatures % (arma::trans(exposures) * fraq);
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );

    estimate = exposures * signatures;
    fraq = data/estimate;

    exposures = exposures % (fraq * arma::trans(signatures));

    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    fraq = data/estimate;
    
    if(t - floor(t/5)*5 == 0){
      for(int r = 0; r<noSignatures; r++) {
        if(range(signatures.row(r)) <= 1e-10){
          Rcout << "collaps";
          exposures.col(r) = arma::randg(genomes);
          signatures.row(r) = arma::randg<arma::rowvec>(mutTypes);
        }
      }
    }
    
  }
  double gkl = error(arma::vectorise(data),arma::vectorise(estimate));
  
  return {exposures, signatures, gkl};
}

// [[Rcpp::export]]
List nmftrain(arma::mat data, arma::mat exposures, arma::mat signatures, arma::mat weight, int iter = 5000) {

  arma::mat estimate = exposures * signatures;
  arma::mat fraq = data/estimate;

  for(int t = 0; t < iter; t++){

    signatures = signatures % (arma::trans(exposures) * fraq);
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );

    estimate = exposures * signatures;
    fraq = data/estimate;

    exposures = exposures % (fraq * arma::trans(signatures));
    arma::colvec exp_sum = sum(exposures,1);
    exposures = exposures.each_col() / exp_sum;
    exposures = weight * exposures;
    exposures = exposures.each_col() % exp_sum;

    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    fraq = data/estimate;

  }
  double gkl = error(arma::vectorise(data),arma::vectorise(estimate));

  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gkl);
  return output;
}

// [[Rcpp::export]]
List nmfspatial(arma::mat data, int noSignatures, arma::mat weight, int maxiter = 10000, double tolerance = 1e-8, int initial = 100, int smallIter = 500) {
  

  arma::colvec obs_sum = sum(data,1);
  
  auto res = nmf1(data, noSignatures, smallIter);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto gklValue = std::get<2>(res);
  
  for(int i = 1; i < initial; i++){
    auto res = nmf1(data, noSignatures, smallIter);
    auto gklNew = std::get<2>(res);
    
    if(gklNew < gklValue){
      gklValue = gklNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
      
    }
  }
  
  arma::mat estimate = exposures * signatures;
  arma::mat fraq = data/estimate;
  
  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);

  for(int t = 0; t < maxiter; t++){

    signatures = signatures % (arma::trans(exposures) * fraq);
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );

    estimate = exposures * signatures;
    fraq = data/estimate;

    exposures = exposures % (fraq * arma::trans(signatures));
    arma::colvec exp_sum = sum(exposures,1);
    exposures = exposures.each_col() / exp_sum;
    exposures = weight * exposures;
    exposures = exposures.each_col() % exp_sum;

    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    fraq = data/estimate;

    gklvalues.at(t) = gklOld;
    
    gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
    
    if (2*std::abs(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance){
      Rcout << "Total iterations:";
      Rcout << t;
      Rcout << "\n";
      break;
    }
    gklOld = gklNew;

    //if(t - floor(t/10)*10 == 0){
    //  Rcout << "Iteration ";
    //  Rcout << t << "\n";
    //}
    
    // if(t - floor(t/10)*10 == 0 || t < 200){
    //   for(int r = 0; r<noSignatures; r++) {
    //     if(range(exposures.col(r)) <= 1e-10){
    //       Rcout << "collaps";
    //       exposures.col(r) = arma::randu(genomes);
    //       signatures.row(r) = arma::randu(mutTypes);
    //       
    //     }
    //     
    //   }
    //}
    
  }
  
  //arma::colvec rsum = sum(signatures,1);
  //exposures = exposures.each_row() % arma::trans(rsum);
  //signatures = signatures.each_col() / rsum;
  
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;
}
