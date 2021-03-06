% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opls_perm.R
\name{opls_perm}
\alias{opls_perm}
\title{OPLS model Y-permutations}
\usage{
opls_perm(smod, n = 10, plot = TRUE, mc = FALSE)
}
\arguments{
\item{smod}{OPLS model of the package \emph{metabom8}}

\item{n, }{num Number of permutations}

\item{plot}{logical, indicating if results should be visualised}

\item{mc}{logical, indicating if tasked should be parallelised using multiple cores}
}
\value{
data.frame with perutation indices
}
\description{
Model validation using Y-permutations
}
\section{}{
NA
}

\examples{
data(covid)
model=opls(X, Y=an$type)
perm=opls_perm(model)
}
\seealso{
Other NMR ++: 
\code{\link{.hotellingsT2}()},
\code{\link{cvanova}()},
\code{\link{es_cdelta}()},
\code{\link{minmax}()},
\code{\link{pca}()},
\code{\link{plotscores}()},
\code{\link{pqn}()},
\code{\link{predict_opls}()},
\code{\link{scRange}()}
}
\author{
Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}

Torben Kimhofer \email{torben.kimhofer@murdoch.edu.au}
}
\concept{NMR ++}
