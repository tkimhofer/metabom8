#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// --- Helper: NIPALS PCA Component ---
// [[Rcpp::export(.nipPcaCompRcpp)]]
Rcpp::List nip_pca_comp_rcpp(Eigen::MatrixXd X) {
  Eigen::VectorXd t_x = X.col(0);
  Eigen::VectorXd t_xold = t_x;
  Eigen::RowVectorXd p_x(X.cols());

  double dd = 1.0;
  int count = 0;

  while(dd > 1e-6 && count < 1000) {
    p_x = (t_x.transpose() * X) / t_x.squaredNorm();
    p_x.normalize();
    t_x = (X * p_x.transpose()) / p_x.squaredNorm();

    if(count > 0) {
      dd = (t_x - t_xold).squaredNorm();
    }
    t_xold = t_x;
    count++;
  }

  Eigen::MatrixXd X_res = X - (t_x * p_x);

  return Rcpp::List::create(
    Rcpp::_["X_res"] = X_res,
    Rcpp::_["t"] = t_x,
    Rcpp::_["p"] = p_x
  );
}

// --- Helper: Multi-Y Weights ---
// [[Rcpp::export]]
Eigen::MatrixXd multiY_Tw_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y, int it_max, double eps) {
  Eigen::MatrixXd W_x(X.cols(), Y.cols());
  for(int i = 0; i < Y.cols(); i++){
    W_x.col(i) = (Y.col(i).transpose() * X).transpose() / Y.col(i).squaredNorm();
  }

  Eigen::MatrixXd T_w(X.cols(), 0);
  Eigen::MatrixXd W_res = W_x;
  double ss_W = W_x.squaredNorm();
  double qss = 1.0;
  int count = 0;

  while(qss > eps && count < it_max) {
    Rcpp::List w_PCA = nip_pca_comp_rcpp(W_res);
    Eigen::VectorXd t_w = Rcpp::as<Eigen::VectorXd>(w_PCA["t"]);

    T_w.conservativeResize(Eigen::NoChange, T_w.cols() + 1);
    T_w.col(T_w.cols() - 1) = t_w;

    W_res = Rcpp::as<Eigen::MatrixXd>(w_PCA["X_res"]);
    qss = t_w.squaredNorm() / ss_W;
    count++;
  }
  return T_w;
}

// --- Core PLS Logic ---
Rcpp::List nip_PLS_comp_core(const Eigen::Ref<const Eigen::MatrixXd>& X,
                             const Eigen::Ref<const Eigen::MatrixXd>& Y,
                             int it_max, double eps) {
  Eigen::VectorXd t_y = Y.col(0);
  Eigen::VectorXd t_x, t_yold, w_x, w_y;
  double dd = 1.0;
  int count = 0;

  while(dd > eps && count < it_max) {
    w_x = (X.transpose() * t_y).normalized();
    t_x = X * w_x;
    w_y = (Y.transpose() * t_x).normalized();
    t_y = Y * w_y;

    if(count > 0) dd = (t_y - t_yold).squaredNorm();
    t_yold = t_y;
    count++;
  }

  Eigen::VectorXd p_x = (X.transpose() * t_x) / t_x.squaredNorm();
  Eigen::VectorXd p_y = (Y.transpose() * t_y) / t_y.squaredNorm();
  Eigen::MatrixXd b = (t_y.transpose() * t_x) / t_x.squaredNorm();

  Eigen::MatrixXd y_pred = (b(0,0) * t_x) * p_y.transpose();

  return Rcpp::List::create(
    Rcpp::_["w_x"] = w_x, Rcpp::_["t_x"] = t_x, Rcpp::_["p_x"] = p_x,
            Rcpp::_["w_y"] = w_y, Rcpp::_["t_y"] = t_y, Rcpp::_["p_y"] = p_y,
                    Rcpp::_["y_pred"] = y_pred,
                    Rcpp::_["x_res"] = X - (t_x * p_x.transpose()), Rcpp::_["b"] = b
  );
}

// [[Rcpp::export(.nipPlsCompRcpp)]]
Rcpp::List nip_PLS_comp_rcpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y, int it_max, double eps) {
  return nip_PLS_comp_core(X, Y, it_max, eps);
}

// --- Orthogonalization ---
// [[Rcpp::export(.orthoGramSchmidtRcpp)]]
Eigen::MatrixXd ortho_gram_schmidt_rcpp(Eigen::VectorXd u, Eigen::MatrixXd v) {
  Eigen::MatrixXd v1 = v;
  double denom = u.squaredNorm();
  if (denom < 1e-12) return v;

  for (int i = 0; i < v.cols(); ++i) {
    double fac = u.dot(v.col(i)) / denom;
    v1.col(i) = v.col(i) - fac * u;
  }
  return v1;
}

// [[Rcpp::export(.nipOplsRcpp)]]
Rcpp::List nip_opls_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y) {
  Rcpp::List rcomp = nip_PLS_comp_core(X, Y, 800, 1e-6);
  Eigen::VectorXd p_x = Rcpp::as<Eigen::VectorXd>(rcomp["p_x"]);
  Eigen::MatrixXd T_w;

  if (Y.cols() > 1) {
    T_w = multiY_Tw_rcpp(X, Y, 800, 1e-6);
  } else {
    T_w = Rcpp::as<Eigen::VectorXd>(rcomp["w_x"]);
  }

  // Orthogonalize
  for (int i = 0; i < T_w.cols(); i++) {
    p_x = ortho_gram_schmidt_rcpp(T_w.col(i), p_x);
  }

  Eigen::VectorXd w_o = p_x.normalized();
  Eigen::VectorXd t_o = (X * w_o) / w_o.squaredNorm();
  Eigen::RowVectorXd p_o = (t_o.transpose() * X) / t_o.squaredNorm();
  Eigen::MatrixXd X_res = X - (t_o * p_o);

  Rcpp::List comp_pred = nip_PLS_comp_core(X_res, Y, 800, 1e-6);

  return Rcpp::List::create(
    Rcpp::_["X_res"] = X_res,
    Rcpp::_["t_p"] = comp_pred["t_x"],
                              Rcpp::_["p_p"] = comp_pred["p_x"],
                                                        Rcpp::_["w_p"] = comp_pred["w_x"],
                                                                                  Rcpp::_["t_o"] = t_o,
                                                                                  Rcpp::_["p_o"] = p_o,
                                                                                  Rcpp::_["w_o"] = w_o,
                                                                                  Rcpp::_["p_y"] = rcomp["p_y"]
  );
}

// [[Rcpp::export(.tssRcpp)]]
double tss_rcpp(const Eigen::Map<Eigen::MatrixXd>& X) {
  return X.squaredNorm();
}



// #include <RcppEigen.h>
//
//
// // [[Rcpp::depends(RcppEigen)]]
//
// using Eigen::Map;                       // 'maps' rather than copies
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;                  // variable size vector, double precision
//
//
// // [[Rcpp::export(.nipPcaCompRcpp)]]
// Rcpp::List nip_pca_comp_rcpp(Eigen::MatrixXd X) {
//
//   Eigen::MatrixXd t_x(X.rows(),1); // scores
//   Eigen::MatrixXd t_xold(X.rows(),1); //  scores (i-1)
//   Eigen::MatrixXd p_x(1, X.rows()); // loadings
//   Eigen::MatrixXd X_res(X.rows(), X.cols()); // X residual
//
//   double dd = 1;
//   int count = 0;
//   Eigen::MatrixXd inter;
//
//   t_x= X.col(0);
//   while( dd > 1e-6 ) {
//     p_x = t_x.transpose() * X / t_x.squaredNorm();
//     p_x = p_x / p_x.norm();
//     t_x = X * p_x.transpose() / p_x.squaredNorm();
//     if( count > 0 ) {
//       inter = t_x - t_xold;
//       dd = inter.squaredNorm();
//     }
//     t_xold = t_x;
//     if ( count > 1000 ) {
//       Rcpp::Rcout << "NIPALS failed to converge after 1000 iterations!\n";
//       break;
//     }
//     count += 1;
//   }
//   X_res = X - (t_x * p_x);
//
//   return Rcpp::List::create(Rcpp::_["X_res"] = X_res,
//                             Rcpp::_["t"] = t_x,
//                             Rcpp::_["p"] = p_x);
// }
//
// //' Summarise Multivariate Y Weights via PCA
// //'
// //' Computes X-weight vectors for each response variable in a multivariate outcome matrix `Y`,
// //' then performs PCA on these weights to extract shared predictive structure.
// //' Useful when `Y` has multiple columns (e.g., for multiclass problems).
// //' Conecptually similar to Canonical Correlation Analysis (CCA) and Multiblock PLS.
// //'
// //' @param X A numeric matrix (n × p): predictor matrix.
// //' @param Y A numeric matrix (n × k): response matrix with k variables.
// //' @param it_max Maximum nb of iterations for NIPALS converge.
// //' @param eps Threshold for sum of squares quotient below which NIPALS is considered converged
// //'
// //' @return A matrix (p × r) of PCA scores, where r is the number of retained components.
// //' @details
// //' The function first computes weights for each column in `Y` using least-squares projection,
// //' forming a weight matrix `W_x`. PCA is then applied to this matrix using NIPALS,
// //' and score components are extracted until the explained variance ratio drops below a threshold.
// //' @examples
// //' # Simulate data: 30 samples, 10 predictors, 3 response variables
// //' set.seed(123)
// //' X <- matrix(rnorm(30 * 10), nrow = 30)
// //' Y <- matrix(rnorm(30 * 3), nrow = 30)
// //' # Compute multivariate Y weights via PCA
// //' T_w <- multiY_Tw_rcpp(X, Y, it_max = 50, eps = 1e-4)
// //' @keywords internal
// //' @export
// // [[Rcpp::export]]
// Eigen::MatrixXd multiY_Tw_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y, int it_max, double eps) {
//
//   Eigen::MatrixXd W_x(X.cols(), Y.cols());
//   Eigen::MatrixXd T_w(X.cols(), 0);
//   Eigen::MatrixXd t_w(X.cols(),1);
//
//   Rcpp::List w_PCA;
//   Eigen::MatrixXd W_res(X.cols(), Y.cols());
//
//   double qss = 1;
//   int count = 0;
//   double ss_T;
//   double ss_W;
//
//   for(int i=0; i<Y.cols(); i++){
//     W_x.col(i) = (Y.col(i).transpose() * X).transpose() /  Y.col(i).squaredNorm();
//   }
//
//   ss_W = ss_T  = W_x.squaredNorm();
//
//   while( qss > eps ) {
//
//     if( count > it_max ) {
//       Rcpp::Rcout << "NIPALS failed to converge, increase itermation max!\n";
//       break;
//     }
//
//     W_res = (count == 0) ? W_x : Rcpp::as<Eigen::MatrixXd>(w_PCA["X_res"]);
//     w_PCA = nip_pca_comp_rcpp(W_res);
//
//     if (!w_PCA.containsElementNamed("t")) {
//       Rcpp::stop("PCA output missing 't'");
//     }
//
//     t_w = w_PCA["t"];
//     T_w.conservativeResize(Eigen::NoChange, (T_w.cols()+1));
//     T_w.col(T_w.cols()-1) = t_w;
//
//     ss_T = t_w.squaredNorm();
//     if (ss_W == 0.0) Rcpp::stop("Cannot divide by zero: total SS of W is zero.");
//     qss = abs(ss_T / ss_W);
//
//     count++;
//   }
//
//   return T_w;
// }
//
//
// Rcpp::List nip_PLS_comp_core(
//     const Eigen::Ref<const Eigen::MatrixXd>& X,
//     const Eigen::Ref<const Eigen::MatrixXd>& Y,
//     int it_max, double eps) {
//
//   Eigen::VectorXd w_x(X.cols());
//   Eigen::VectorXd w_y(Y.cols());
//
//   Eigen::VectorXd t_x(X.rows());
//   Eigen::VectorXd t_y(Y.rows());
//   Eigen::VectorXd t_yold(X.rows());
//
//   Eigen::VectorXd p_x(X.cols());
//   Eigen::VectorXd p_y(Y.cols());
//
//   Eigen::MatrixXd b; // covariance x and y scores (inner relation)
//
//   Eigen::MatrixXd y_pred(X.rows(),1);    // x scores
//   Eigen::MatrixXd inter; // store interim results
//
//   Eigen::MatrixXd x_pred(X.rows(),X.cols());  ; // t*p
//   Eigen::MatrixXd x_res(X.rows(), X.cols());  ; // X where component rm
//
//   double dd = 1;
//   int count = 0;
//
//   t_y = Y.col(0);
//
//   while( dd > eps ) {
//     w_x = X.transpose() * t_y;
//     w_x.normalize();
//
//     t_x = X * w_x;
//
//     w_y = Y.transpose() * t_x;
//     w_y.normalize();
//
//
//     t_y = Y * w_y;
//     t_y.normalize();
//
//     if( count > 0 ) {
//       dd = (t_y - t_yold).squaredNorm() / t_y.squaredNorm();
//     }
//
//     t_yold=t_y;
//
//     if( count > it_max ) {
//       Rcpp::Rcout << "NIPALS failed to converge, increase it_max!\n";
//       break;
//     }
//     count ++;
//   }
//
//   w_x = X.transpose() * t_x;
//   w_x.normalize();
//
//   w_y = Y.transpose() * t_y;
//   w_y.normalize();
//
//   p_x = X.transpose() * t_x / t_x.squaredNorm();
//   p_y = Y.transpose() * t_y / t_y.squaredNorm();
//
//
//   b = t_y.transpose() * t_x / t_x.squaredNorm(); // linear inner relation (regression form): t_y = b * t_x
//
//   y_pred = b(0,0) * t_x * p_y.transpose();
//
//   x_pred = t_x * p_x.transpose();   // (n×1)(1×p)
//   x_res = X - x_pred;
//
//   // rescaling once again for pls consistency (vip)
//   double scale = (t_y.transpose() * t_x)(0,0) / t_x.squaredNorm();
//   t_y = scale * t_x;
//
//
//   return  Rcpp::List::create(
//     Rcpp::_["w_x"] = w_x,
//     Rcpp::_["t_x"] = t_x,
//     Rcpp::_["p_x"] = p_x,
//     Rcpp::_["w_y"] = w_y,
//     Rcpp::_["t_y"] = t_y,
//     Rcpp::_["p_y"] = p_y,
//     Rcpp::_["y_pred"] = y_pred,
//     Rcpp::_["x_res"] = x_res,
//     Rcpp::_["b"] = b);
//
// }
//
//
//
// // [[Rcpp::export(.nipPlsCompRcpp)]]
// Rcpp::List nip_PLS_comp_rcpp(
//     const Eigen::Map<Eigen::MatrixXd> X,
//     const Eigen::Map<Eigen::MatrixXd> Y,
//     int it_max,
//     double eps
// ){
//   return nip_PLS_comp_core(X, Y, it_max, eps);
// }
//
//
//
//
// // gram schmidt orthogonalisation
// // [[Rcpp::export(.orthoGramSchmidtRcpp)]]
// Eigen::MatrixXd ortho_gram_schmidt_rcpp(Eigen::MatrixXd u, Eigen::MatrixXd v) {
//   // Projects each column of v orthogonal to u
//   // for vectors u and v: projects v on u and subtracts from v -> result is u_ortho (u*u_ortho = 0)
//   if (u.cols() != 1)
//     Rcpp::stop("u must be a column vector");
//
//   if (u.rows() != v.rows())
//     Rcpp::stop("Dimension mismatch between u and v");
//
//   Eigen::MatrixXd v1 = v;
//   double denom = u.squaredNorm();
//
//   if (denom == 0)
//     Rcpp::stop("Cannot orthogonalise with zero-norm vector");
//
//   for (int i = 0; i < v.cols(); ++i) {
//     double fac = (u.transpose() * v.col(i))(0, 0) / denom;
//     v1.col(i) = v.col(i) - fac * u;
//   }
//
//   return v1;
// }
//
// // [[Rcpp::export(.nipOplsRcpp)]]
// Rcpp::List nip_opls_rcpp(SEXP X_, SEXP Y_) {
//
//   if(!Rf_isMatrix(X_) || !Rf_isReal(X_))
//     Rcpp::stop("X must be numeric matrix");
//
//   if(!Rf_isMatrix(Y_) || !Rf_isReal(Y_))
//     Rcpp::stop("Y must be numeric matrix");
//
//   Rcpp::NumericMatrix Xmat(X_);
//   Rcpp::NumericMatrix Ymat(Y_);
//
//   Eigen::Map<Eigen::MatrixXd> X(Xmat.begin(),
//                                 Xmat.nrow(),
//                                 Xmat.ncol());
//
//   Eigen::Map<Eigen::MatrixXd> Y(Ymat.begin(),
//                                 Ymat.nrow(),
//                                 Ymat.ncol());
//
//   const int n = X.rows();
//   const int p = X.cols();
//   const int k = Y.cols();
//
//   Eigen::MatrixXd t_o(n, 1);
//   Eigen::MatrixXd p_o(1, p);
//   Eigen::MatrixXd X_res(n, p);
//
//   Rcpp::List rcomp = nip_PLS_comp_core(X, Y, 800, 1e-6);
//
//   if (!rcomp.containsElementNamed("p_x")) Rcpp::stop("rcomp missing 'p_x'");
//   // Eigen::MatrixXd p_x = rcomp["p_x"];
//   // p_x = p_x.transpose();
//   // Eigen::VectorXd p_x = rcomp["p_x"];
//   Eigen::VectorXd p_x = Rcpp::as<Eigen::VectorXd>(rcomp["p_x"]);
//   Eigen::MatrixXd T_w;
//
//   if (Y.cols() > 1) {
//     T_w = multiY_Tw_rcpp(X, Y, 800, 1e-6);
//   } else {
//     // Cast here too!
//     T_w = Rcpp::as<Eigen::MatrixXd>(rcomp["w_x"]);
//   }
//
//   Eigen::MatrixXd w_o(p_x.rows(), 1);
//
//   if (Y.cols() > 1) {
//     T_w = multiY_Tw_rcpp(X, Y, 800, 1e-6);
//   } else {
//     if (!rcomp.containsElementNamed("w_x")) Rcpp::stop("rcomp missing 'w_x'");
//     T_w = rcomp["w_x"];
//   }
//
//   for (int i = 0; i < T_w.cols(); i++) {
//     p_x = ortho_gram_schmidt_rcpp(T_w.col(i), p_x);
//   }
//
//   w_o = p_x;
//
//   if (w_o.norm() == 0) Rcpp::stop("Orthogonal weight vector has zero norm");
//
//   w_o = w_o / w_o.norm();
//   t_o = X * w_o / w_o.squaredNorm();
//   p_o = t_o.transpose() * X / t_o.squaredNorm();
//
//   // X_res = X - (t_o * p_o);
//   X_res = (X - (t_o * p_o));
//
//   Eigen::MatrixXd X_res_eval = X_res;   // already eval'd but safe
//   Eigen::MatrixXd Y_eval     = Y;       // ← THIS is the key fix
//
//
//   Rcpp::List comp_pred = nip_PLS_comp_core(X_res_eval, Y_eval, 800, 1e-6);
//
//   if (!comp_pred.containsElementNamed("t_x")) Rcpp::stop("comp_pred missing 't_x'");
//   if (!comp_pred.containsElementNamed("p_x")) Rcpp::stop("comp_pred missing 'p_x'");
//   if (!comp_pred.containsElementNamed("w_x")) Rcpp::stop("comp_pred missing 'w_x'");
//
//
//   return Rcpp::List::create(
//     Rcpp::_["X_res"] = X_res,
//     Rcpp::_["t_p"] = comp_pred["t_x"],
//     Rcpp::_["p_p"] = comp_pred["p_x"],
//     Rcpp::_["w_p"] = comp_pred["w_x"],
//     Rcpp::_["t_o"] = t_o,
//     Rcpp::_["p_o"] = p_o,
//     Rcpp::_["w_o"] = w_o,
//     Rcpp::_["p_y"] = rcomp["p_y"]
//   );
// }
// // Rcpp::List nip_opls_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y) {
// //   Eigen::MatrixXd t_o(X.rows(), 1); // orth scores
// //   Eigen::MatrixXd p_o(1, X.cols()); // orth loadings
// //   Eigen::MatrixXd X_res(X.rows(), X.cols()); // orth filtered X
// //
// //   Rcpp::List rcomp = nip_PLS_comp_rcpp(X, Y);
// //   Eigen::MatrixXd p_x = rcomp["p_x"];
// //   Eigen::MatrixXd w_o(p_x.rows(), 1);
// //   Eigen::MatrixXd T_w;
// //
// //   if(Y.cols()>1) {
// //     T_w = multiY_Tw_rcpp(X, Y);
// //   }else{
// //     T_w = rcomp["w_x"];
// //   }
// //
// //   // orthogonalise p to each column in T
// //   for ( int i = 0; i < T_w.cols() ; i++ ) {
// //     p_x = ortho_gram_schmidt_rcpp(T_w.col(i), p_x.transpose());
// //   }
// //   w_o = p_x;
// //
// //   // calc loadings
// //   w_o = w_o / w_o.norm();
// //   t_o = X * w_o / w_o.squaredNorm();
// //   p_o = t_o.transpose() * X / t_o.squaredNorm();
// //
// //   X_res = X - ( t_o * p_o );
// //
// //   // calc predictive component with orthogonal-filtered X data
// //   Rcpp::List comp_pred = nip_PLS_comp_rcpp(X_res, Y);
// //
// //
// //
// //   return  Rcpp::List::create( Rcpp::_["X_res"] = X_res,
// //                               Rcpp::_["t_p"] = comp_pred["t_x"],
// //                                                         Rcpp::_["p_p"] = comp_pred["p_x"],
// //                                                                                   Rcpp::_["w_p"] = comp_pred["w_x"],
// //                                                                                                             Rcpp::_["t_o"] = t_o,
// //                                                                                                             Rcpp::_["p_o"] = p_o,
// //                                                                                                             Rcpp::_["w_o"] =  w_o,
// //                                                                                                             Rcpp::_["p_y"] =  rcomp["p_y"]);
// //
// // }
//
//
// // tss_rcpp(X=X)
// // [[Rcpp::export(.tssRcpp)]]
// Rcpp::NumericVector tss_rcpp(Eigen::MatrixXd X) {
//
//   Rcpp::NumericVector out(1);
//   out = X.squaredNorm();
//
//   return out;
// }
//
//
// // opls_mod input for this is a list with w_xo, p_xo, (matrix with columns, rows equating to the number of orthogonal components, resp) w_xp, p_xp, p_yp (a single one for predictive component)
// // // this is for OPLS_pred where data has been filtered and there is only a single component (predictive one in OPLS, not suitable for pls with more than one component)
//
// [[Rcpp::export(.oplsPredRcpp)]]
Rcpp::List opls_pred_rcpp(Rcpp::List opls_mod, Rcpp::List pred_mod, Eigen::MatrixXd Xnew) {


  // Rcpp::print(Rcpp::wrap(Xnew.rows()));
  // Rcpp::print(Rcpp::wrap(Xnew.cols()));

  // 1. Extract Orthogonal Components
  // Use MatrixXd if there could be multiple orthogonal components
  Eigen::MatrixXd w_xo = Rcpp::as<Eigen::MatrixXd>(opls_mod["w_o"]);
  Eigen::MatrixXd p_xo = Rcpp::as<Eigen::MatrixXd>(opls_mod["p_o"]);

  // Rcpp::Rcout << "p_xo-row: "<< p_xo.rows()<< std::endl;
  // Rcpp::Rcout << "p_xo-col: "<< p_xo.cols()<< std::endl;

  int n_samples = Xnew.rows();
  int n_ortho   = w_xo.cols(); // Number of orthogonal components

  Eigen::MatrixXd t_xo_new(n_samples, n_ortho);
  Eigen::MatrixXd x_res = Xnew;

  // if (w_xo.rows() != Xnew.cols())
  //   Rcpp::stop("w_o has wrong number of rows");
  //
  // if (p_xo.rows() != Xnew.cols())
  //   Rcpp::stop("p_o has wrong number of rows");

  // Rcpp::Rcout << "starting loop\n";
  // 2. Filter with previous orthogonal components
  for(int i = 0; i < n_ortho; i++) {
    // Score: t = (X * w) / (w' * w)
    // Note: Use .col(i) for both w and p if they are stored column-wise (standard in R)
    double w_norm_sq = w_xo.col(i).squaredNorm();
    t_xo_new.col(i) = (x_res * w_xo.col(i)) / w_norm_sq;
    // Rcpp::Rcout << "t_xo_new-row: "<< t_xo_new.rows()<< std::endl;
    // Rcpp::Rcout << "t_xo_new-col: "<< t_xo_new.cols()<< std::endl;

    // Update residual: X_res = X_res - (t * p')
    // Check if p_xo stores components as columns (usually true in R)
    x_res -= t_xo_new.col(i) * p_xo.row(i);
    // Rcpp::Rcout << "starting next it of loop\n";
  }

  // Rcpp::Rcout << "done looping\n";

  // 3. Predictive Component extraction
  Eigen::MatrixXd w_p = Rcpp::as<Eigen::MatrixXd>(pred_mod["w_p"]);
  Eigen::MatrixXd p_y = Rcpp::as<Eigen::MatrixXd>(pred_mod["p_y"]);
  double b = Rcpp::as<double>(pred_mod["b"]);

  // 4. Calculate prediction
  // Rcpp::Rcout << "Calculate prediction\n";
  // Use MatrixXd for t_pred to be safe, even if it's 1 column
  Eigen::MatrixXd t_pred = x_res * w_p;

  // Rcpp::Rcout << "Calculate prediction1\n";
  // y_pred = b * t_pred * p_y'
  Eigen::MatrixXd y_pred = (b * t_pred) * p_y.transpose();

  // Rcpp::Rcout << "Returning\n";
  return Rcpp::List::create(
    Rcpp::_["Xres"]     = x_res,
    Rcpp::_["t_pred"]   = t_pred,
    Rcpp::_["y_pred"]   = y_pred,
    Rcpp::_["t_xo_new"] = t_xo_new
  );
}


// [[Rcpp::export(.plsPredRcpp)]]
Rcpp::List pls_pred_rcpp(Rcpp::List pls_mod, Eigen::MatrixXd Xnew) {
  // Safe casting
  Eigen::MatrixXd w_x = Rcpp::as<Eigen::MatrixXd>(pls_mod["w_x"]);
  Eigen::MatrixXd p_x = Rcpp::as<Eigen::MatrixXd>(pls_mod["p_x"]);
  Eigen::MatrixXd p_y = Rcpp::as<Eigen::MatrixXd>(pls_mod["p_y"]);
  Eigen::MatrixXd b_obj = Rcpp::as<Eigen::MatrixXd>(pls_mod["b"]);
  double b = b_obj(0,0);

  Eigen::MatrixXd t_pred = Xnew * w_x;
  Eigen::MatrixXd E_h = Xnew - t_pred * p_x;
  Eigen::MatrixXd res = b * t_pred * p_y.transpose();

  return Rcpp::List::create(Rcpp::_["Xres"] = E_h,
                            Rcpp::_["t_pred"] = t_pred,
                            Rcpp::_["y_pred"] = res);
}
