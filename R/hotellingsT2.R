
#' @title Hotelling's T2 ellipse in 2D
#' @description This function is used to calculate hotellings T2 ellipse
#' @param x num vector descrbing dimension 1
#' @param y, num vector descrbing dimension 2
#' @param alpha num, probability level
#' @references Geladi, P and Kowalski, B.R. (1986), Partial least squares and regression: a tutorial. \emph{Analytica Chimica Acta}, 185, 1-17.
#' @return data.frame containing H.T2 ellipse cooredinates
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom ellipse ellipse
#' @importFrom stats cov
#
.hotellingsT2=function(x, y, alpha=0.95){
  SD <- cov(cbind(x, y), use='complete.obs')
  el <- ellipse(SD, centre = colMeans(cbind(x, y), na.rm = T), level = alpha)
  colnames(el) <- c("V1", "V2")
  xlim <- c(min((c(el[, 1], x))), max((c(el[, 1], x))))
  xlim <- xlim + c(diff(range(xlim)) * -0.05, diff(range(xlim)) *
                     +0.05)
  ylim <- c(min((c(el[, 2], y))), max((c(el[, 2], y))))
  ylim <- ylim + c(diff(range(ylim)) * -0.05, diff(range(ylim)) *
                     +0.05)
  df <- as.data.frame(el)

  return(df)
}

