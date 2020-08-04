#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec calcPhi(double ph0, double& ph1, int& le) {
  arma::vec v = arma::linspace<arma::vec>(0, le, le+1);
  arma::vec ang =ph0 + ((ph1 * v) / le);
  return ang;
}



// [[Rcpp::export]]
arma::vec phase1d(arma::vec& sp_re, arma::vec& sp_im, double ph0, double ph1){
  int le=sp_re.n_elem-1;
  arma::vec phi=calcPhi(ph0, ph1, le);
  arma::vec out= (sp_re  % cos(phi)) - (sp_im % sin(phi));
  return out;
}


// [[Rcpp::export]]
Rcpp::List phaseTsp(arma::vec& sp_re, arma::vec& sp_im, arma::vec& ppm, arma::vec& ph0, double& ph1, arma::uvec& idx_tsp){

  arma::vec s_ph;

  int tsp_max;

  arma::vec bound;
  arma::vec comp;
  arma::vec iid;
  int start;
  int end;
  arma::uvec iiu;
  arma::uvec iis;
  arma::vec out(ph0.n_elem-1);

  float inter;
  int bmin;

  for(int i=0; i<ph0.n_elem; ++i)
  {
    s_ph=phase1d(sp_re, sp_im, (i), ph1);
    tsp_max=arma::index_max(s_ph.elem(idx_tsp));

  // this needs double checking
    inter = (idx_tsp.n_elem - tsp_max);
    bound << tsp_max << inter << arma::endr;

    bmin=bound.min();
    start=tsp_max-bound.min();
    end=tsp_max+bound.min();

    // // prep spec in tsp reagion with boundary
    iid=arma::linspace<arma::vec>(start, end, end-start);
    iiu=arma::conv_to<arma::uvec>::from(iid);
    arma::cout << (i) << arma::endl;
    //iis=idx_tsp.elem(iiu);

    //arma::sum(abs(s_ph.elem(iiu) - reverse(s_ph.elem(iiu))));
    out[i] = arma::sum(abs(s_ph.elem(iiu) - reverse(s_ph.elem(iiu))));
  }


//  return bound;
  return Rcpp::List::create(Rcpp::Named("s_ph") = s_ph,
                            Rcpp::Named("iiu")       = iiu,
                            Rcpp::Named("tsp_max")       = tsp_max,
                             Rcpp::Named("idx_tsp.max")  = idx_tsp.max(),
                             Rcpp::Named("idx_tsp")  = idx_tsp,
                             Rcpp::Named("bound")=bound,
                             Rcpp::Named("start")=start,
                             Rcpp::Named("end")=end,
                            Rcpp::Named("out")=out);
}


// //
// // // [[Rcpp::export]]
// // Rcpp::List fastLm(const arma::mat& X, const arma::colvec& y) {
// //   int n = X.n_rows, k = X.n_cols;
// //
// //   arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
// //   arma::colvec res  = y - X*coef;           // residuals
// //
// //   // std.errors of coefficients
// //   double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
// //
// //   arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));
// //
// //   return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
// //                             Rcpp::Named("stderr")       = std_err,
// //                             Rcpp::Named("df.residual")  = n - k);
// // }
// //
// //
// //
// //
