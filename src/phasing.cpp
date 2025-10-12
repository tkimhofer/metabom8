#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
// [[Rcpp::export(.calcPhi)]]
arma::vec calcPhi(double ph0, double& ph1, int& le) {
  arma::vec v = arma::linspace<arma::vec>(0, le, le+1);
  arma::vec ang =ph0 + ((ph1 * v) / le);
  return ang;
}

//' @keywords internal
// [[Rcpp::export(.phase1d)]]
arma::vec phase1d(arma::vec& sp_re, arma::vec& sp_im, double ph0, double ph1){
//  Rcpp::Rcout << "phasing 1D..." << std::endl;
  int le=sp_re.n_elem-1;
  arma::vec phi= calcPhi(ph0, ph1, le);
  arma::vec out= (sp_re  % cos(phi)) - (sp_im % sin(phi));
 // Rcpp::Rcout << "done." << std::endl;
  return out;
}

//' @keywords internal
// [[Rcpp::export(.phaseTsp)]]
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

  for(int i=0; i<ph0.n_elem; ++i)
  {
    s_ph=phase1d(sp_re, sp_im, ph0(i), ph1);
    tsp_max=arma::index_max(abs(s_ph.elem(idx_tsp)));

    inter = (idx_tsp.n_elem - tsp_max);
    //bound << tsp_max << inter << arma::endr;
    bound.set_size(2);
    bound(0) = static_cast<double>(tsp_max);
    bound(1) = static_cast<double>(inter);

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
}

//' @keywords internal
// [[Rcpp::export(.zerofill)]]
arma::vec zerofil(arma::vec fid, const int zf, int le_ori){
  int n_zero= (pow(2, log2(le_ori) +zf)) - fid.n_elem;

  arma::vec  ze = arma::zeros<arma::vec>(n_zero);

  arma::vec out = join_cols(fid, ze);
  return out;
}

//' @keywords internal
// [[Rcpp::export(.cplxFft)]]
arma::cx_vec cplxFft(arma::vec fid){

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

//' @keywords internal
// [[Rcpp::export(.defineChemShiftPpm)]]
arma::vec defineChemShiftPpm(const float sf_mhz, const float sw_hz, const int n_sp_re, const float dref, bool ref){
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

//' @title Calibrate
//' @return Shift-adjusted ppm vector
//' @keywords internal
// [[Rcpp::export(.calibTsp)]]
arma::vec calibTsp(arma::vec spec, arma::vec ppm){

// find max close to zero
  arma::uvec ids = find((ppm >= -0.1 && ppm <= 0.1));
  int idx_max= arma::index_max(spec(ids));

  arma::vec ppm_tsp = ppm(ids);
  float shift=ppm_tsp(idx_max);

  arma::vec ppm_shift=ppm-shift;
  return ppm_shift;
}
