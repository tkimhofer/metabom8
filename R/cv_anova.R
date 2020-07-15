#' @title Cross-validated ANOVA
#' @export
#' @description Significance testing for OPLS models
#' @param smod OPLS_metabom8 object of the package \emph{metabom8}.
#' @details The CV-ANOVA diagnostic formally compares the fit of two models to the same data by the size of their residuals. The function tests the residuals of the linear regression between cross-validated scores of the predictive O-PLS component and the response Y, with the variation of Y around its mean. The p value is derived from an F-test with the null hypothesis of equal residuals of the two models. For detailed information  on p value colculation see refrence further below.
#' @references Eriksson, L, et al. (2008) CV-ANOVA for significance testing of PLS and OPLS models. \emph{Journal of Chemometrics}, 22, 594-600.
#' @return \emph{data.frame} describing ANOVA stats incl, p value
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom stats lm
cvanova<-function(smod){

  if(class(smod)[1]!='OPLS_metabom8'){ stop('Function requires OPLS_metabom8 object.')}

  # calculate degrees of freedom
  df_total<-nrow(smod@Y$dummy)-1
  df_reg<-(smod@nPC)*2
  df_res<-df_total-df_reg

  # cacl linear models
  lm2=lm(smod@Y$dummy[,1] ~ 1)
  res2=lm2$residuals

  # calc sum or squares of residuals
  lm1=lm(smod@Y$dummy[,1] ~ 1+ smod@t_pred_cv[,1] )
  res1=lm1$residuals

  ss_total<-sum(res2^2)
  ss_reg<-sum(res1^2)
  ss_res<-ss_total-ss_reg

  df<-c(df_total, df_reg, df_res)
  ss<-c(ss_total, ss_reg, ss_res)

  # ss normalised by df
  ms<-ss/df

  F_val<-ms[2]/ms[3]

  p_val<-1-pf(F_val, df[2], df[3])

  out=data.frame('SS'=c(format(ss, scientific = T), NA), 'DF'=c(df, NA), 'MS'=c(format(ms, scientific = T), NA), 'F_value'=c(rep(NA, 3), F_val),   'p_value'=c(rep(NA, 3),format(p_val, scientific = T)), row.names = c('Total corrected', 'Regression', 'Residual', 'RESULT'))

  return(out)
}
