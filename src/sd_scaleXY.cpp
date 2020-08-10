#include <Rcpp.h>
using namespace Rcpp;

//' @title Column-wise standard deviation and mean for a matrix
//' @param X num matrix
//' @keyword internal
//' @return list: 1. sd (num vec), 2. mean (sd vec)
//' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
// [[Rcpp::export(.sdRcpp)]]
List sd_rcpp(NumericMatrix X) {

  NumericVector sd_out(X.cols());
  NumericVector mean_out(X.cols());

  for (int i = 0; i < X.cols(); i++) {
    sd_out[i] = sd(X( _ , i ));
    mean_out[i] = mean(X( _ , i ));
  }

  return List::create(Named("sd") = sd_out, Named("mean") = mean_out);

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


// test=scale_rcpp(X, idc=0:(nrow(X)-1), center=T, scale_type=0)


//' @title Column-wise matrix scaling
//' @export
//' @param X num matrix
//' @param idc int row indices of X
//' @param center bool mean centering
//' @param scale_type int 0: no scaling, 1: SD scaling, 2: Pareto scaling
//' @return list: 1. scale X matrix, 2. mean (sd vec), 3: sd (num vec)
//' @keyword internal
//' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
// [[Rcpp::export(.sdMatRcpp)]]
List scale_rcpp(NumericMatrix X, IntegerVector idc, bool center,  int scale_type) {

  // scale_type: 0 - none, 1 - UV, 2 - Pareto

  NumericVector sd_out(X.cols());
  NumericVector mean_out(X.cols());

  NumericMatrix X1(idc.length(), X.ncol());

  for (int i = 0; i < idc.length(); i++) {
    X1.row(i) =  X.row(idc(i));
  }

  //X1 = X( idc , _ );

  NumericMatrix X_out(X.nrow(), X.ncol());

  if (center==true) {

    if (scale_type == 1) {

      for (int i = 0; i < X.cols(); i++) {

        sd_out[i] = sd(X1( _ , i ));
        mean_out[i] = mean(X1( _ , i ));

        if(sd_out[i] ==0) {sd_out[i] = 1;} else{
          X_out(_,i) = ( X(_,i) - mean_out[i ] ) / sd_out[i];
        }

      }

    }

    if (scale_type == 2) {

      for (int i = 0; i < X.cols(); i++) {

        sd_out[i] = sd(X1( _ , i ));
        mean_out[i] = mean(X1( _ , i ));

        if(sd_out[i] ==0) {sd_out[i];} else{
          X_out(_,i) = ( X(_,i) - mean_out[i] ) / pow(sd_out[i], 0.5);
        }

      }

    }

    if (scale_type == 0) {

      for (int i = 0; i < X.cols(); i++) {

        mean_out[i] = mean(X1( _ , i ));

        X_out(_,i) = ( X(_,i) - mean_out[i] );

      }

    }

  }else{
    if (scale_type == 1) {

      for (int i = 0; i < X.cols(); i++) {

        sd_out[i] = sd(X1( _ , i ));

        if(sd_out[i] ==0) {continue;}
        X_out(_,i) = ( X(_,i) / sd_out[i] );

      }

    }

    if (scale_type == 2) {

      for (int i = 0; i < X.cols(); i++) {

        sd_out[i] = sd(X1( _ , i ));

        if(sd_out[i] ==0) {continue;}
        X_out(_,i) = ( X(_,i) / pow(sd_out[i], 0.5) );

      }

    }

  }

  return List::create(Named("X_prep") = X_out, Named("mean") = mean_out, Named("sd") = sd_out);

}
