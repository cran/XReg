
\name{plotXRinvreg}
\alias{plotXRinvreg}
\title{Plots inverse of extreme regression model}
\description{ 
\code{plotXRinvreg} Returns a plot of the inverse regression function, as a matrix plot. Rows of matrix are variables and columns correspond to miniumum terms.
This representation of the regression function probably only works well as a graphical display for up to about 5 unique variables and 5 terms.}

\usage{plotXRinvreg(x,xrout,shading=F,quant=F,both=F,qlev=c(.2,.4,.6,.8),
                      spec.ylim=F,mat.ylim=NULL,level.set=NA)
}

\arguments{
  \item{x}{matrix of predictors}
  \item{xrout}{output from XR regression function XReg}
  \item{shading}{if shading=T shade the directions of the rules in panel associated with sequence of rules {x:f(x)>=q} }  
  \item{quant}{if quant=T, plot rules in terms of quantiles or ranked predictions. Note, models will not appear linear univarite.}
  \item{both}{if both=T, plot in linear form, but put quantiles on the top of plot} 
  \item{qlev}{if both=T, these are the labels for the quantiles put on the top of the plot}
  \item{spec.ylim}{if spec.ylim=F, do not specify the limits of the y axis of the plot. Note, in this case the same variable 
   different plots could have different axes. If spec.ylim=T, then the user must provide the limits 
   y-axis for each of the variables.}
  \item{mat.ylim}{This gives the y-axis limits (min,max) for each of the variables, it is a 2xp matrix, where 
  each column corresponds to a predictor variable. It is only used if spec.ylim=T.}
  \item{level.set}{If this is set to a numeric value, a highlighted line will be added to the 
  plot to indicate associated rule}
}


\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
  # See XReg
}
\keyword{ regression }
\keyword{ nonlinear }