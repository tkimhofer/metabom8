#include <Rcpp.h>
using namespace Rcpp;



//' @title Column-wise standard deviation and mean for a matrix
//'  Welford algorithm (single loop for mean/sd)
//' @param X num matrix
//' @keywords internal
//' @return list: 1. sd (num vec), 2. mean (sd vec)
 // [[Rcpp::export(.sdRcpp)]]
 Rcpp::List sd_rcpp(const Rcpp::NumericMatrix& X) {

   int n = X.nrow();
   int p = X.ncol();

   Rcpp::NumericVector mean(p);
   Rcpp::NumericVector var(p);

   for (int j = 0; j < p; j++) {

     double m = 0.0;
     double M2 = 0.0;

     for (int i = 0; i < n; i++) {
       double x = X(i,j);
       double delta = x - m;
       m += delta / (i + 1);
       M2 += delta * (x - m);
     }

     mean[j] = m;

     if (n > 1)
       var[j] = std::sqrt(M2 / (n - 1));
     else
       var[j] = NA_REAL;
   }

   return Rcpp::List::create(
     Rcpp::Named("sd")   = var,
     Rcpp::Named("mean") = mean
   );
 }

//
// List sd_rcpp(NumericMatrix X) {
//
//   int n = X.nrow();
//   int p = X.ncol();
//
//   NumericVector sd_out(p);
//   NumericVector mean_out(p);
//
//   for (int j = 0; j < p; j++) {
//
//     double m = 0.0;
//     for (int i = 0; i < n; i++)
//       m += X(i,j);
//     m /= n;
//
//     mean_out[j] = m;
//
//     double s = 0.0;
//     for (int i = 0; i < n; i++) {
//       double d = X(i,j) - m;
//       s += d*d;
//     }
//
//     sd_out[j] = std::sqrt(s/(n-1));
//   }
//
//   return List::create( _["sd"]   = sd_out, _["mean"] = mean_out
//   );
// }



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


// test=scale_rcpp(X, idc=0:(nrow(X)-1), center=T, scale_type=0)


//' @title Column-wise matrix scaling  (C++ backend)
//' @param X num matrix
//' @param idc int row indices of X
//' @param center logical mean centering
//' @param scale_type int 0: no scaling, 1: SD scaling, 2: Pareto scaling
//' @return list: 1. scale X matrix, 2. mean (sd vec), 3: sd (num vec)
//' @keywords internal
// [[Rcpp::export(.scaleMatRcpp)]]
Rcpp::List scale_rcpp(const Rcpp::NumericMatrix& X,
                      const Rcpp::IntegerVector& idc,
                      bool center,
                      int scale_type) {

  int n  = X.nrow();
  int p  = X.ncol();
  int nt = idc.length();

  Rcpp::NumericVector mean_out(p);
  Rcpp::NumericVector sd_out(p);
  Rcpp::NumericMatrix X_out(n, p);

  const double eps = 1e-12;

  for(int j = 0; j < p; j++){

    double m  = 0.0;
    double M2 = 0.0;

    // --- Welford on TRAINING rows only ---
    if(center || scale_type > 0){

      for(int k = 0; k < nt; k++){

        double x = X(idc[k], j);
        double delta = x - m;

        m += delta / (k + 1);
        M2 += delta * (x - m);
      }

      mean_out[j] = m;
    }

    double denom = 1.0;

    // --- scaling ---
    if(scale_type > 0){

      double var = (nt > 1) ? M2/(nt - 1) : 0.0;

      if(var < eps){

        denom     = 1.0;   // do not scale flat vars
        sd_out[j] = 0.0;

      } else {

        double sd = std::sqrt(var);
        sd_out[j] = sd;

        if(scale_type == 1)       // UV
          denom = sd;
        else if(scale_type == 2)  // Pareto
          denom = std::sqrt(sd);
      }
    }

    // --- apply to ALL rows ---
    for(int i = 0; i < n; i++){

      double val = X(i,j);

      if(center)
        val -= m;

      X_out(i,j) = val / denom;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("X_prep") = X_out,
    Rcpp::Named("mean")   = mean_out,
    Rcpp::Named("sd")     = sd_out
  );
}

// List scale_rcpp(NumericMatrix X,
//                IntegerVector idc,
//                bool center,
//                int scale_type){
//
//  int n  = X.nrow();
//  int p  = X.ncol();
//  int nt = idc.length();
//
//  NumericVector sd_out(p);
//  NumericVector mean_out(p);
//  NumericMatrix X_out(n,p);
//
//  for(int j=0; j<p; j++){
//
//    double m = 0.0;
//
//    // compute mean on training rows
//    for(int i=0; i<nt; i++)
//      m += X(idc[i], j);
//
//    m /= nt;
//    mean_out[j] = m;
//
//    double s = 0.0;
//
//    if(scale_type > 0){
//      for(int i=0; i<nt; i++){
//        double d = X(idc[i],j) - m;
//        s += d*d;
//      }
//      s = std::sqrt(s/(nt-1));
//      if(s == 0) s = 1;
//      sd_out[j] = s;
//    }
//
//    // apply transform to full matrix
//    for(int i=0; i<n; i++){
//
//      if(center){
//        if(scale_type == 1)
//          X_out(i,j) = (X(i,j) - m)/s;
//        else if(scale_type == 2)
//          X_out(i,j) = (X(i,j) - m)/std::sqrt(s);
//        else
//          X_out(i,j) = (X(i,j) - m);
//      } else {
//        if(scale_type == 1)
//          X_out(i,j) = X(i,j)/s;
//        else if(scale_type == 2)
//          X_out(i,j) = X(i,j)/std::sqrt(s);
//        else
//          X_out(i,j) = X(i,j);
//      }
//
//    }
//  }
//
//  return List::create(
//    _["X_prep"] = X_out,
//    _["mean"]   = mean_out,
//    _["sd"]     = sd_out
//  );
// }
//
// List scale_rcpp(NumericMatrix X, IntegerVector idc, bool center,  int scale_type) {
//
//   // scale_type: 0 - none, 1 - UV, 2 - Pareto
//
//   NumericVector sd_out(X.cols());
//   NumericVector mean_out(X.cols());
//
//   NumericMatrix X1(idc.length(), X.ncol());
//
//   for (int i = 0; i < idc.length(); i++) {
//     X1.row(i) =  X.row(idc(i));
//   }
//
//   //X1 = X( idc , _ );
//
//   NumericMatrix X_out(X.nrow(), X.ncol());
//
//   if (center==true) {
//
//     if (scale_type == 1) {
//
//       for (int i = 0; i < X.cols(); i++) {
//
//         sd_out[i] = sd(X1( _ , i ));
//         mean_out[i] = mean(X1( _ , i ));
//
//         if(sd_out[i] ==0) {sd_out[i] = 1;} else{
//           X_out(_,i) = ( X(_,i) - mean_out[i ] ) / sd_out[i];
//         }
//
//       }
//
//     }
//
//     if (scale_type == 2) {
//
//       for (int i = 0; i < X.cols(); i++) {
//
//         sd_out[i] = sd(X1( _ , i ));
//         mean_out[i] = mean(X1( _ , i ));
//
//         if(sd_out[i] ==0) {sd_out[i];} else{
//           X_out(_,i) = ( X(_,i) - mean_out[i] ) / pow(sd_out[i], 0.5);
//         }
//
//       }
//
//     }
//
//     if (scale_type == 0) {
//
//       for (int i = 0; i < X.cols(); i++) {
//
//         mean_out[i] = mean(X1( _ , i ));
//
//         X_out(_,i) = ( X(_,i) - mean_out[i] );
//
//       }
//
//     }
//
//   }else{
//     if (scale_type == 1) {
//
//       for (int i = 0; i < X.cols(); i++) {
//
//         sd_out[i] = sd(X1( _ , i ));
//
//         if(sd_out[i] ==0) {continue;}
//         X_out(_,i) = ( X(_,i) / sd_out[i] );
//
//       }
//
//     }
//
//     if (scale_type == 2) {
//
//       for (int i = 0; i < X.cols(); i++) {
//
//         sd_out[i] = sd(X1( _ , i ));
//
//         if(sd_out[i] ==0) {continue;}
//         X_out(_,i) = ( X(_,i) / pow(sd_out[i], 0.5) );
//
//       }
//
//     }
//
//   }
//
//   return List::create(Named("X_prep") = X_out, Named("mean") = mean_out, Named("sd") = sd_out);
//
// }
