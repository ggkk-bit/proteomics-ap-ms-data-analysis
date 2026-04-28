#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export(.GetDICE)]]
NumericMatrix GetDICE(NumericMatrix mat) {
  int nc = mat.ncol();
  int rstart = 0;
  int rend = mat.nrow();
  NumericMatrix rmat(nc, nc);
  for (int c1 = 0; c1 < nc; c1 ++){
    for(int c2 = 0; c2 < c1; c2++){
      double q = 0;
      double x = 0;
      double y = 0;
      for(int r = rstart; r < rend; r++){
        if((mat(r, c1) == 1) & (mat(r, c2) == 1)){
          q += 1;
        } else if((mat(r, c1) == 1) & (mat(r, c2) == 0)){
          x += 1;
        } else if((mat(r, c1) == 0) & (mat(r, c2) == 1)){
          y += 1;
        }
      }
      rmat(c1, c2) = (2*q)/(2*q + x + y);
    }
  }
  return(rmat);
}
