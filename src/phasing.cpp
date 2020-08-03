#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec calcPhi(double ph0, double& ph1, int& le) {
  arma::vec v = arma::linspace<arma::vec>(0, le, le+1);
  arma::vec ang =ph0 + ((ph1 * v) / le);
  return ang;
}

//
// [[Rcpp::export]]
arma::vec phase1d(arma::vec& sp_re, arma::vec& sp_im, double ph0, double ph1){
  int le=sp_re.n_elem-1;
  arma::vec phi=calcPhi(ph0, ph1, le);
  arma::vec out= (sp_re  % cos(phi)) - (sp_im % sin(phi));
  return out;
}
//
//
// [[Rcpp::export]]
int phaseTsp(arma::vec& sp_re, arma::vec& sp_im, arma::vec& ppm, arma::vec& ph0, double& ph1, arma::uvec& idx_tsp){

  arma::mat::iterator it     = ph0.begin();
  arma::mat::iterator it_end = ph0.end();

  arma::vec s_ph;

  int tsp_max;
  arma::vec bound;

  arma::vec comp;
  arma::vec iid;
  int start;
  int end;
  arma::uvec iiu;
  arma::uvec iis;
  arma::vec out;

  int bmin;

  for(; it != it_end; ++it)
  {
    s_ph=phase1d(sp_re,sp_im, (*it), ph1);
    tsp_max=arma::index_max(s_ph.elem(idx_tsp));
    //
    bound=(tsp_max, idx_tsp.max());

    bmin=bound.min();
    start=tsp_max-bound.min();
    end=tsp_max+bound.min();
    // // prep spec in tsp reagion with boundary
    //
    // iid=arma::linspace<arma::vec>(start, end, end-start);
    // iiu=arma::conv_to<arma::uvec>::from(iid);
    // // arma::cout << iiu << arma::endl;
    // iis=idx_tsp.elem(iiu);
    //
    // out[(*it)], arma::sum(abs(s_ph.elem(iis) - reverse(s_ph.elem(iis))));
  }

  return bmin;
}


//
//
//
// {
//   range <- sort(range, decreasing = T)
//   which(ppm <= range[1] & ppm >= range[2])
// }

//
//
//
//
//
//
// //   .phase_tsp<-function(sp_re, sp_im, ppm, phi=seq(0, pi, by=0.01), psi=0){
// //   idx=get.idx(c(-0.05, 0.05), ppm)
// //
// //     tsp_max<-which.max(sp_re[idx])
// //
// // # peak symmetry
// //     out=sapply(phi, function(p, ps=psi){
// //       s_ph<-.phase1d(sp_re, sp_im, ang=c(p,ps), demo=F)
// // #browser()
// //       tsp_max<-which.max(s_ph[idx])
// // #print(tsp_max)
// //       cent_cap<-min(tsp_max, length(idx)-tsp_max)-2
// //       iid=(tsp_max-cent_cap):(tsp_max+cent_cap)
// //
// //       sum(abs(s_ph[idx][iid]-rev(s_ph[idx][iid])))^2
// //     })
// //
// //       return(phi[which.min(out)])
// //
// //   }
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
// //
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
