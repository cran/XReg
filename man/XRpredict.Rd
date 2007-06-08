
\name{predictfun}
\alias{predictfun}
\title{Predictions for extreme regression models}
\description{ 
\code{predictfun} This gives predictions but also give assigments to the active sets 
 for data given in the xpred matrix. The default is not to give assignments. 
 }
 
\usage{predictfun(xpred,varlist,beta0list,lower=-Inf,upper=Inf,activeset=F)  

}


\arguments{

 \item{xpred}{xmatrix for predictions}
 \item{varlist}{model form list (described in XR regression function)}
 \item{beta0list}{the beta list from XR regression function}
 \item{lower}{lower bound for predictions}
 \item{upper}{upper bound for predictions}
 \item{activeset}{if T denote where observations are in terms of active set}


}
\value{
  a list containing: 
  \item{predfit}{predicted values}

}

\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
 #See XReg 
}
\keyword{ regression }
\keyword{ nonlinear }