#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export(.GetPPN)]]
NumericMatrix GetPPN(NumericMatrix mat) {
  int nc = mat.ncol();
  int rstart = 0;
  int rend = mat.nrow();
  NumericMatrix rmat(nc, nc);
  for (int c1 = 0; c1 < nc; c1 ++){
    for(int c2 = 0; c2 < c1; c2++){
      int sMinXY = 0;
      for(int r = rstart; r < rend; r++){
        int minXY = 0;
        if(mat(r, c1) <= mat(r, c2)){
          minXY = mat(r, c1);
        }else{
          minXY = mat(r, c2);
        }
        sMinXY += minXY;
      }
      rmat(c1, c2) = sMinXY;
    }
  }
  return(rmat);
}
