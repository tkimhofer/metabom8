#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

// [[Rcpp::export(.oplsPredRcpp)]]
Rcpp::List opls_pred_rcpp(Rcpp::List opls_mod, Rcpp::List pred_mod, Eigen::MatrixXd Xnew) {

  Eigen::MatrixXd w_xo = Rcpp::as<Eigen::MatrixXd>(opls_mod["w_o"]);
  Eigen::MatrixXd p_xo = Rcpp::as<Eigen::MatrixXd>(opls_mod["p_o"]);

  int n_samples = Xnew.rows();
  int n_ortho   = w_xo.cols(); // Number of orthogonal components

  Eigen::MatrixXd t_xo_new(n_samples, n_ortho);
  Eigen::MatrixXd x_res = Xnew;

  for(int i = 0; i < n_ortho; i++) {

    double w_norm_sq = w_xo.col(i).squaredNorm();
    t_xo_new.col(i) = (x_res * w_xo.col(i)) / w_norm_sq;
    x_res -= t_xo_new.col(i) * p_xo.row(i);
  }

  Eigen::MatrixXd w_p = Rcpp::as<Eigen::MatrixXd>(pred_mod["w_p"]);
  Eigen::MatrixXd p_y = Rcpp::as<Eigen::MatrixXd>(pred_mod["p_y"]);
  double b = Rcpp::as<double>(pred_mod["b"]);

  Eigen::MatrixXd t_pred = x_res * w_p;
  Eigen::MatrixXd y_pred = (b * t_pred) * p_y.transpose();

  return Rcpp::List::create(
    Rcpp::_["Xres"]     = x_res,
    Rcpp::_["t_pred"]   = t_pred,
    Rcpp::_["y_pred"]   = y_pred,
    Rcpp::_["t_xo_new"] = t_xo_new
  );
}


// [[Rcpp::export(.plsPredRcpp)]]
Rcpp::List pls_pred_rcpp(Rcpp::List pls_mod, Eigen::MatrixXd Xnew) {
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
