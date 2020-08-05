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
arma::vec phaseTsp(arma::vec& sp_re, arma::vec& sp_im, arma::vec& ppm, arma::vec& ph0, double& ph1, arma::uvec idx_tsp){

  arma::vec s_ph;
  arma::vec bound;
  arma::vec comp;
  arma::vec iid;
  arma::uvec iiu;
  arma::uvec iis;
  arma::vec out(ph0.n_elem-1);
  arma::vec idx_tsp1;

  float inter;
  int bmin;
  int tsp_max;
  int start;
  int end;
  int test;

  for(int i=0; i<ph0.n_elem; ++i)
  {
    s_ph=phase1d(sp_re, sp_im, ph0(i), ph1);
    tsp_max=arma::index_max(abs(s_ph.elem(idx_tsp)));

    test=idx_tsp.n_elem;
    inter = (idx_tsp.n_elem - tsp_max);
    bound << tsp_max << inter << arma::endr;

    bmin=bound.min()-2;
    start=tsp_max-bmin;
    end=tsp_max+bmin;

    iid=arma::linspace<arma::vec>(start, end, end-start+1);
    iiu=arma::conv_to<arma::uvec>::from(iid);
    iis=idx_tsp.elem(iiu);
    out[i] = arma::sum(abs(s_ph.elem(iis) - reverse(s_ph.elem(iis))));
  }

  float ang_best = ph0(arma::index_min(out));
  arma::vec spec_phase=phase1d(sp_re, sp_im, ang_best, 0);

  return spec_phase;

// //  return bound;
//   return Rcpp::List::create(Rcpp::Named("idx_tsp") = idx_tsp,
//                             Rcpp::Named("idx_tsp1") = idx_tsp1,
//                             Rcpp::Named("s_ph") = s_ph,
//                             Rcpp::Named("tsp_max")       = tsp_max,
//                             Rcpp::Named("test")  = test,
//                             Rcpp::Named("inter")  = inter,
//                             Rcpp::Named("bound")=bound,
//                             Rcpp::Named("bmin")=bmin,
//                             Rcpp::Named("start")=start,
//                             Rcpp::Named("end")=end,
//                             Rcpp::Named("iid")=iid,
//                             Rcpp::Named("iiu") = iiu,
//                             Rcpp::Named("iis") = iis,
//                             Rcpp::Named("out")=out,
//                             Rcpp::Named("spec_phase")=spec_phase,
//                             Rcpp::Named("ang_best")=ang_best
//                             );

}

// [[Rcpp::export]]
arma::vec zerofil(arma::vec fid, const int zf, int le_ori){
  int n_zero= (pow(2, log2(le_ori) +zf)) - fid.n_elem;

  arma::vec  ze = arma::zeros<arma::vec>(n_zero);

  arma::vec out = join_cols(fid, ze);
  return out;
}


// [[Rcpp::export]]
arma::cx_vec cplx_fft(arma::vec fid){

  arma::vec idx_re = arma::linspace<arma::vec>(0, fid.n_elem-2,( fid.n_elem/2));
  arma::uvec idx_reu=arma::conv_to<arma::uvec>::from(idx_re);

  arma::vec idx_im = idx_re+1;
  arma::uvec idx_imu=arma::conv_to<arma::uvec>::from(idx_im);

  arma::cx_vec fid_cplx = arma::cx_vec(fid.elem(idx_reu),fid.elem(idx_imu));

  arma::cx_vec spec = fft(fid_cplx);


   // re-arrange axis to get a spectrum from min to max
   int split = spec.n_elem / 2;

   //arma::cout << split << arma::endl;

   arma::vec idx_up = arma::linspace<arma::vec>(split, spec.n_elem-1,(spec.n_elem-split));
   arma::uvec idx_upu=arma::conv_to<arma::uvec>::from(idx_up);

   arma::vec idx_dwn = arma::linspace<arma::vec>(0, split-1, split);
   arma::uvec idx_dwnu=arma::conv_to<arma::uvec>::from(idx_dwn);

   arma::uvec idx_reord = join_cols(idx_upu, idx_dwnu);

   arma::cx_vec out = spec(idx_reord);

   return out;
}




// [[Rcpp::export]]
arma::vec defineChemShiftppm(const float sf_mhz, const float sw_hz, const int n_sp_re, const float dref, bool ref){
  float dist=sw_hz/n_sp_re;
  arma::vec pps = arma::regspace<arma::vec>(0, dist, (sw_hz-(dist/2)));
  arma::vec ppm = (pps - (pps(n_sp_re/2)-(dref*sf_mhz))) / sf_mhz;

  if(ref==true){
    float lower = round(ppm[0])-0.001;
    float upper = trunc(ppm[ppm.n_elem-1])+0.001;
    arma::uvec ids = find((ppm >= lower && ppm <= upper));
    ppm=ppm(ids);
  }

  return ppm;
}

// # define ppm
// ppmDist<-pars$a_SW_h[s]/length(sp_re)
//   pps=seq(0, pars$a_SW_h[s], by=ppmDist)[-1]
// ppm=(pps-(pps[length(sp_re)/2]-(4.79 * pars$a_SFO1[s])))/pars$a_SFO1[s]



// [[Rcpp::export]]
arma::vec calibTsp(arma::vec spec,arma::vec ppm){

// find max close to zero
  arma::uvec ids = find((ppm >= -0.1 && ppm <= 0.1));
  int idx_max= arma::index_max(spec(ids));

  arma::vec ppm_tsp = ppm(ids);
  float shift=ppm_tsp(idx_max);

  arma::vec ppm_shift=ppm-shift;
  return ppm_shift;
}

//
// idx <- get.idx(c(-0.2, 0.2), ppm)
//   zero.ppm <- which.min(abs(ppm[idx]))
//   maxInt <- array()
//   for (i in 1:nrow(X)) {
//     maxInt[i] <- which.max(X[i, idx])
//   }
//   Int.corr <- zero.ppm - maxInt


