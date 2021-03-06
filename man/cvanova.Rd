% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvanova.R
\name{cvanova}
\alias{cvanova}
\title{Cross-validated ANOVA for O-PLS models}
\usage{
cvanova(smod)
}
\arguments{
\item{smod}{OPLS_metabom8 object of the package \emph{metabom8}.}
}
\value{
\emph{data.frame} describing ANOVA stats incl, p value
}
\description{
Significance testing for OPLS models
}
\details{
The CV-ANOVA diagnostic formally compares the fit of two models to the same data by the size of their residuals. The function tests the residuals of the linear regression between cross-validated scores of the predictive O-PLS component and the response Y, with the variation of Y around its mean. The p value is derived from an F-test with the null hypothesis of equal residuals of the two models. For detailed information  on p value colculation see refrence further below.
}
\references{
Eriksson, L, et al. (2008) CV-ANOVA for significance testing of PLS and OPLS models. \emph{Journal of Chemometrics}, 22, 594-600.
}
\seealso{
Other NMR ++: 
\code{\link{.hotellingsT2}()},
\code{\link{es_cdelta}()},
\code{\link{minmax}()},
\code{\link{opls_perm}()},
\code{\link{pca}()},
\code{\link{plotscores}()},
\code{\link{pqn}()},
\code{\link{predict_opls}()},
\code{\link{scRange}()}
}
\author{
\email{torben.kimhofer@murdoch.edu.au}
}
\concept{NMR ++}
