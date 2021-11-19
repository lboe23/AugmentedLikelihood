#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cmat(NumericVector id, NumericVector time, IntegerVector result, double phi1, double phi0, double negpred) {
  NumericVector utime = unique(time);
  utime.sort();
  int J = utime.size(), nsub = unique(id).size(), nobs = id.size(), i, j;
  NumericVector rid(nobs);
  NumericMatrix Cm(nsub, J + 1), Dm(nsub, J + 1);
  utime.push_back(utime[J-1] + 1);
  //initiate with 1 for Cm
  std::fill(Cm.begin(), Cm.end(), 1);
  //Calculate Cm
  j = 0;
  for (i = 1; i < nobs; i++) {
    if (id[i] != id[i-1]) j++;
    rid[i] = j;
  }
  for (j = 0; j <= J; j++) {
    for (i = 0; i < nobs; i++) {
      if (result[i]==0) {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= 1 - phi1;
        } else {
          Cm(rid[i], j) *= phi0;
        }
      } else {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= phi1;
        } else {
          Cm(rid[i], j) *= 1 - phi0;
        }
      }
    }
  }
  return Cm;
}

