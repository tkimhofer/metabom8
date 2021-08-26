#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision



// data(iris)
// X=as.matrix(iris[,c(1:4)])
// Y=cbind(as.numeric(as.factor(iris[,5])) )
// idx=which(Y<3)
// X1=scale(X[idx,])
// Y1=scale(cbind(Y[idx,]))
// out = nip_pca_comp_rcpp(X1)

// [[Rcpp::export(.nipPcaCompRcpp)]]
Rcpp::List nip_pca_comp_rcpp(Eigen::MatrixXd X) {

  Eigen::MatrixXd t_x(X.rows(),1); // scores
  Eigen::MatrixXd t_xold(X.rows(),1); //  scores (i-1)
  Eigen::MatrixXd p_x(1, X.rows()); // loadings
  Eigen::MatrixXd X_res(X.rows(), X.cols()); // X residual

  double dd = 1;
  int count = 0;
  Eigen::MatrixXd inter;

  t_x= X.col(0);
  while( dd > 1e-6 ) {
    //Rcpp::Rcout << "Iter ";
    //Rcpp::Rcout << count;
    //Rcpp::Rcout << '\t';
    //Rcpp::Rcout << dd;
    //Rcpp::Rcout << '\n';

    p_x = t_x.transpose() * X / t_x.squaredNorm();
    p_x = p_x / p_x.norm();
    t_x = X * p_x.transpose() / p_x.squaredNorm();
    if( count > 0 ) {
      inter = t_x - t_xold;
      dd = inter.squaredNorm();
    }
    t_xold = t_x;
    if ( count > 1000 ) {
      Rcpp::Rcout << "NIPALS failed to converge after 1000 iterations!\n";
      break;
    }
    count += 1;
  }
  X_res = X - (t_x * p_x);

  return Rcpp::List::create(Rcpp::_["X_res"] = X_res,
                            Rcpp::_["t"] = t_x,
                            Rcpp::_["p"] = p_x);
}





// data(iris)
// X=as.matrix(iris[,c(1:4)])
// Y=cbind(as.numeric(as.factor(iris[,5])) )
// idx=which(Y<3)
// X1=scale(X[idx,])
// Y1=MetaboMate:::create_dummy_Y(iris[idx,5])[[1]]
// out = multiY_Tw_rcpp(X1, Y1)
// [[Rcpp::export(.multiY_TwRcpp)]]
Eigen::MatrixXd multiY_Tw_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y) {

  // if Y has multiple columns, prouduce weights for each Y column and summarise these by PCA
  Eigen::MatrixXd W_x(X.cols(), Y.cols()); // X weights matrix
  Eigen::MatrixXd T_w(X.cols(), 0); // PCA scores matrix (multiple compoments) of PCA of W
  Eigen::MatrixXd t_w(X.cols(),1); // PCA scores (single component) of PCA of W

  Rcpp::List w_PCA; // results PCA
  Eigen::MatrixXd W_res(X.cols(), Y.cols()); // residuals in Weight matrix from PCA

  double qss = 1;
  int count = 0;
  double ss_T;
  double ss_W;

  // for multilevel, produce weights
  for(int i=0; i<Y.cols(); i++){
    W_x.col(i) = Y.col(i).transpose() * X /  Y.col(i).squaredNorm();
  }

  ss_W = ss_T  = W_x.squaredNorm();

  // estimate PCA the principal components of W as long as the ratio of SS of current score vector t divided by SS of W is larger than given threshold, typicaly 10e-10
  while( qss > 10e-10 ) {
    if( count == 0 ) { W_res = W_x; }else { W_res = w_PCA["X_res"]; }
    w_PCA = nip_pca_comp_rcpp(W_res);
    t_w = w_PCA["t"];
    T_w.conservativeResize(Eigen::NoChange, (T_w.cols()+1));
    T_w.col(T_w.cols()-1) = t_w;
    if( count > 300 ) {
      Rcpp::Rcout << "NIPALS failed to converge after 300 iterations!\n";
      break;
    }
    ss_T = t_w.squaredNorm();
    qss = abs(ss_T / ss_W);
    count += 1 ;
  }

  return T_w;

}


// tt=nip_reg_PLS_comp(X, Y)
// ss=nipcomp_rcpp(X, Y)


// data(iris)
// X=as.matrix(iris[,c(1:4)])
// Y=cbind(as.numeric(as.factor(iris[,5])) )
// idx=which(Y<4)
// X1=scale(X[idx,])
// Y1=MetaboMate:::create_dummy_Y(iris[idx,5])[[1]]
// out = .nipPlsCompRcpp(X1, Y1)
// [[Rcpp::export(.nipPlsCompRcpp)]]
Rcpp::List nip_PLS_comp_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y) {

  Eigen::MatrixXd w_x(1, X.cols()); // x weights
  Eigen::MatrixXd w_y(1, Y.cols()); // y weights

  Eigen::MatrixXd t_x(X.rows(),1);    // x scores
  Eigen::MatrixXd t_y(Y.rows(),1);  ; // y scores

  Eigen::MatrixXd p_y(Y.cols(),1); // x loadings
  Eigen::MatrixXd p_x(1, X.cols()); // y loadings
  Eigen::MatrixXd b; // covariance x and y scores (inner relation)

  Eigen::MatrixXd y_pred(X.rows(),1);    // x scores
  Eigen::MatrixXd t_yold(X.rows(),1); // x scores (i-1)
  Eigen::MatrixXd inter;// store interim results

  Eigen::MatrixXd x_pred(X.rows(),X.cols());  ; // t*p
  Eigen::MatrixXd x_res(X.rows(), X.cols());  ; // X where component rm

  double dd = 1;
  int count = 0;

  t_y = Y.col(0); // keep in mind: indexing starts at zero

  while( dd > 1e-10 ) {
    w_x =  t_y.transpose() * X /  t_y.squaredNorm();
    w_x = w_x.transpose() / w_x.norm();
    t_x = X * w_x / w_x.squaredNorm();
    w_y = t_x.transpose() * Y / t_x.squaredNorm();
    t_y = Y * w_y.transpose() / w_y.squaredNorm();
    if( count > 0 ) {
      inter = t_y - t_yold;
      dd = inter.squaredNorm() / t_y.squaredNorm();
    }
    t_yold=t_y;
    if( count > 300 ) {
      Rcpp::Rcout << "NIPALS failed to converge after 300 iterations!\n";
      break;
    }
    count += 1;
  }

  p_x = t_x.transpose() * X / t_x.squaredNorm();
  p_x = p_x / p_x.norm();

  p_y = Y.transpose() * t_y / t_y.squaredNorm();
  b = t_y.transpose() * t_x / t_x.squaredNorm(); // linear inner relation (regression form): t_y = b * t_x

  y_pred = b(0,0) * t_x * p_y.transpose();

  x_pred = t_x * p_x;
  x_res = X - x_pred;


  return  Rcpp::List::create(
    Rcpp::_["w_x"] = w_x,
    Rcpp::_["t_x"] = t_x,
    Rcpp::_["p_x"] = p_x,
    Rcpp::_["w_y"] = w_y,
    Rcpp::_["t_y"] = t_y,
    Rcpp::_["p_y"] = p_y,
    Rcpp::_["y_pred"] = y_pred,
    Rcpp::_["x_res"] = x_res,
    Rcpp::_["b"] = b);

}



// out=ortho_gram_schmidt_rcpp(matrix(c(3,1)), matrix(c(2,2)))

// gram-schmnidt orthogonalisation (proj) for vectors u and v: project v on u and subtract from v -> result is u_ortho (u*u_ortho = 0)
// [[Rcpp::export(.orthoGramSchmidtRcpp)]]
Eigen::MatrixXd ortho_gram_schmidt_rcpp(Eigen::MatrixXd u, Eigen::MatrixXd v) {

  Eigen::MatrixXd v1(u.rows(), 1);
  Eigen::MatrixXd facs(1,1);

  facs = (( u.transpose() * v ) / u.squaredNorm()) ;
  v1 =  v - (u*facs(0,0));

  return v1;

}

// data(iris)
// X=scale(as.matrix(iris[,1:4]))
// Y=scale(MetaboMate:::create_dummy_Y(iris[,5])[[1]])
//
  // tt=nip_opls_mlevcomp_rcpp(X, Y)
// ts=NIPALS_OPLS_component_mulitlevel(X, Y)

// // amat = matrix(rnorm(1000), ncol=50)
// // bmat = matrix(scale(rep(c(1:20), 2)), ncol=2)
// // nip_opls_mlevcomp_rcpp(amat, bmat)


// data(iris)
// X=as.matrix(iris[,c(1:4)])
// Y=cbind(as.numeric(as.factor(iris[,5])) )
// idx=which(Y<4)
// X1=scale(X[idx,])
// Y1=MetaboMate:::create_dummy_Y(iris[idx,5])[[1]]
// out = nip_opls_rcpp(X1, Y1)

// [[Rcpp::export(.nipOplsRcpp)]]
Rcpp::List nip_opls_rcpp(Eigen::MatrixXd X, Eigen::MatrixXd Y) {
  Eigen::MatrixXd t_o(X.rows(), 1); // orth scores
  Eigen::MatrixXd p_o(1, X.cols()); // orth loadings
  Eigen::MatrixXd X_res(X.rows(), X.cols()); // orth filtered X

  Rcpp::List rcomp = nip_PLS_comp_rcpp(X, Y);
  Eigen::MatrixXd p_x = rcomp["p_x"];
  Eigen::MatrixXd w_o(p_x.rows(), 1);
  Eigen::MatrixXd T_w;

  if(Y.cols()>1) {
    T_w = multiY_Tw_rcpp(X, Y);
  }else{
    T_w = rcomp["w_x"];
  }

  // orthogonalise p to each column in T
  for ( int i = 0; i < T_w.cols() ; i++ ) {
    p_x = ortho_gram_schmidt_rcpp(T_w.col(i), p_x.transpose());
  }
  w_o = p_x;

  // calc loadings
  w_o = w_o / w_o.norm();
  t_o = X * w_o / w_o.squaredNorm();
  p_o = t_o.transpose() * X / t_o.squaredNorm();

  X_res = X - ( t_o * p_o );

  // calc predictive component with orthogonal-filtered X data
  Rcpp::List comp_pred = nip_PLS_comp_rcpp(X_res, Y);



  return  Rcpp::List::create( Rcpp::_["X_res"] = X_res,
                              Rcpp::_["t_p"] = comp_pred["t_x"],
                                                        Rcpp::_["p_p"] = comp_pred["p_x"],
                                                                                  Rcpp::_["w_p"] = comp_pred["w_x"],
                                                                                                            Rcpp::_["t_o"] = t_o,
                                                                                                            Rcpp::_["p_o"] = p_o,
                                                                                                            Rcpp::_["w_o"] =  w_o,
                                                                                                            Rcpp::_["p_y"] =  rcomp["p_y"]);

}






// tss_rcpp(X=X)

// [[Rcpp::export(.tssRcpp)]]
Rcpp::NumericVector tss_rcpp(Eigen::MatrixXd X) {

  Rcpp::NumericVector out(1);
  out = X.squaredNorm();

  return out;
}


  // opls_mod input for this is a list with w_xo, p_xo, (matrix with columns, rows equating to the number of orthogonal components, resp) w_xp, p_xp, p_yp (a single one for predictive component)
  // // this is for OPLS_pred where data has been filtered and there is only a single component (predictive one in OPLS, not suitable for pls with more than one component)

  // [[Rcpp::export(.oplsPredRcpp)]]
  Rcpp::List opls_pred_rcpp(Rcpp::List opls_mod, Rcpp::List pred_mod, Eigen::MatrixXd Xnew) {


    // filter with previous orthogonal components

    Eigen::MatrixXd w_xo = opls_mod["w_o"]; // colvec
    Eigen::MatrixXd p_xo =  opls_mod["p_o"]; // rowvec

    Eigen::MatrixXd t_xo_new(Xnew.rows(), w_xo.cols()); // check dimsnsions of p_xo, correct here and below
    Eigen::MatrixXd x_res = Xnew;

    for( int i = 0; i < w_xo.cols(); i++) {
      t_xo_new.col(i) = x_res * w_xo.col(i) / w_xo.col(i).squaredNorm();
      x_res = x_res - ( t_xo_new.col(i) * p_xo );
    }



    // perform perdiction with predictive component using orthogoanl filtered data

    Eigen::MatrixXd w_x = pred_mod["w_x"]; // colvec
    Eigen::MatrixXd p_x =  pred_mod["p_x"]; // rowvec
    Eigen::MatrixXd p_y =  pred_mod["p_y"]; // colvec
    double b = pred_mod["b"];

    Eigen::MatrixXd t_pred(Xnew.rows(),1);
    // Eigen::MatrixXd E_h(Xnew.rows(), Xnew.cols());
    Eigen::MatrixXd y_pred(Xnew.rows(), 1);

    t_pred = x_res * w_x;
    // E_h =   x_res - t_pred * p_x;

    y_pred = b * t_pred * p_y.transpose();


    return Rcpp::List::create( Rcpp::_["Xres"] = x_res, // filtered data, w/o predictive component
                               Rcpp::_["t_pred"] = t_pred,
                               Rcpp::_["y_pred"] = y_pred,
                               Rcpp::_["t_xo_new"] = t_xo_new
    );

  }


// [[Rcpp::export(.plsPredRcpp)]]
Rcpp::List pls_pred_rcpp(Rcpp::List pls_mod, Eigen::MatrixXd Xnew) {
  // this returns pls predictions from a single component PLS that was calculated before
  Eigen::MatrixXd w_x = pls_mod["w_x"];
  Eigen::MatrixXd p_x =  pls_mod["p_x"];
  Eigen::MatrixXd p_y =  pls_mod["p_y"];
  double b = pls_mod["b"];

  Eigen::MatrixXd t_pred(Xnew.rows(),1);
  Eigen::MatrixXd E_h(Xnew.rows(), Xnew.cols());
  Eigen::MatrixXd res(Xnew.rows(), 1);

  t_pred = Xnew * w_x;
  E_h =   Xnew - t_pred * p_x;

  res = b * t_pred * p_y.transpose();


  return Rcpp::List::create( Rcpp::_["Xres"] = E_h,
                             Rcpp::_["t_pred"] = t_pred,
                             Rcpp::_["y_pred"] = res
  );

}



