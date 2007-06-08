
\name{plotXRreg}
\alias{plotXRreg}
\title{Plots component linear models in extreme regression model}
\description{ 
\code{plotXRreg} This function returns a plot of the regression function and data points 
as a matrix plot. Rows of matrix are variables and columns correspond to miniumum terms.
This representation of the regression function probably only works well as a graphical display up to 
about 5 unique variables and 5 terms. It also tells how many observations are indicated by each rule
and gived the overall indicator of the observations indicated by the rule.}

\usage{plotXRreg(y,x,xrout,allpoints=F)

}

\arguments{
  \item{y}{response }
  \item{x}{matrix of predictors}
  \item{xrout}{output from XR regression function XReg}
  \item{allpoints}{if allpoints=T, plot all points (but highlight those corresponding to the local univariate function}
}


\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
  # See XReg
# }
\keyword{ regression }
\keyword{ nonlinear }