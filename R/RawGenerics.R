#  apodisation functions
.em<-function(n, lb){
  idx<-seq(n)
  out<-exp(-( idx * lb * pi) / (length(idx)))
  minmax(out)
}

# .gmod<-function(n, a, b){
#   idx<-seq(n)
#   out<-exp(((-a*idx) - (b*idx^2))  / n)
#   minmax(out)
# }

.cosine<-function(n){
  idx<-seq(n)
  out=cos(pi*idx/n)
  minmax(out)
}

.sine<-function(n){
  idx=seq(n)
  out<-sin(pi*idx/n)
  minmax(out)
}

.sem<-function(n, lb=1.5){
  idx=seq(n)
  out<-sin((pi*idx)/n) * exp(-( idx * lb * pi) / n)
  minmax(out)
}

.fidApodisationFct<-function(n, pars){

  if(is.null(pars$fun) || !pars$fun %in% c('uniform', 'exponential', 'cosine', 'sine', 'sem')){stop('Check apodisation function argument (fun).')}

  switch(pars$fun,
         'uniform'={afun<-rep(1, n)},
         'exponential'={ if( 'lb' %in% names(pars)) {afun<-.em(n, pars$lb)} else{stop('Check aposation fct arguments: Exponention function requires lb parameter')}},
         #'modGauss'={},
         #'expDampJmod'={},
         'cosine'={afun<-.cosine(n)},
         'sine'={afun<-.sine(n)},
         #'sineShift'={if( all(c('offs', 'end', 'exp') %in% names(pars))){afun<-.sineMod(n, pars$offs, pars$end, pars$exp)} else{stop('Check aposation fct arguments: sineShift function requires offs, end, exp parameters')}},
         #'triangle'={},
         'sem'={if( 'lb' %in% names(pars)) {afun<-.sem(n, pars$lb)} else{ stop('Check aposation fct arguments: SEM function requires lb parameter')}}
         )

  if(!is.null(pars$plot) && pars$plot && exists('afun')){ plot(afun, type='l', main=paste('Apodisation function:', pars$fun))}

  return(afun)

  }


#
#
# .phase_tsp<-function(sp_re, sp_im, ppm, phi=seq(0, pi, by=0.01), psi=0){
#   idx=get.idx(c(-0.05, 0.05), ppm)
#
#   tsp_max<-which.max(sp_re[idx])
#
#   # peak symmetry
#   out=sapply(phi, function(p, ps=psi){
#     s_ph<-.phase1d(sp_re, sp_im, ang=c(p,ps), demo=F)
#     #browser()
#     tsp_max<-which.max(s_ph[idx])
#     #print(tsp_max)
#     cent_cap<-min(tsp_max, length(idx)-tsp_max)-2
#     iid=(tsp_max-cent_cap):(tsp_max+cent_cap)
#
#     sum(abs(s_ph[idx][iid]-rev(s_ph[idx][iid])))
#   })
#
#   return(phi[which.min(out)])
#   #return(out)
# }
#
#
#
# .phase1d<-function(sp_re, sp_im, ang=NULL, demo=T){
#
#
#   if(is.null(ang)){
#     ang<-expand.grid(-180:0, 0:360)
#     out=apply(ang, 1, function(x, le=length(sp_re), s_re=sp_re, s_im=sp_im){
#       phi=.calc_phi(x[1], x[2], le)
#       # perform phasing
#       s_phase=(s_re * cos(phi)) - (s_im * sin(phi))
#       return((s_phase))
#
#     })
#
#     if(demo==T){
#       return(list(ang, out))
#     }else{
#       return(cbind(ang, out))
#     }
#
#
#   }else{
#
#     phi<-.calc_phi(ang[1], ang[2], le=length(sp_re))
#     s_phase<-(sp_re * cos(phi)) - (sp_im * sin(phi))
#     return(s_phase)
#   }
#
# }
#
# .calc_phi<-function(ph0, ph1, le){
#   ph0+(ph1*(seq(le)/le))
# }
#
#



#
# # correct for FID offset using a baseline correction when decayed signal is not centered at zero (that is due to quadrature imbalance and is called a direct current (DC) offset)
#
# fid_bc<-function(fid){
#
#   n<-length(fid)
#   idx=(n-n*1/100):n
#   me=mean(fid[idx])
#
#   if(me>1){
#     print(me); stop('implement baseline correction');
#   }
#
#
# }
