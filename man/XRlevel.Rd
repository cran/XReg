
\name{level.set}
\alias{level.set}
\title{Level sets for extreme regression models}
\description{ 
\code{level.set} This function gives the decision rule for a given threshold. It also tells how many observations are indicated by each rule
and gives the overall indicator of the observations indicated by the rule. Currently this function is only written for >= rules; however, with it
one should be able to specify the desired inequality }

\usage{level.set(thresh, output,pred=F,xpred=NULL,lower=-Inf,upper=Inf)

}

\arguments{
  \item{thresh}{threshold or target output for rule {x:f(x)>=thresh} }
  \item{output}{data structure from XR regression function}
  \item{pred}{determine which observations are consistent with rule}
  \item{xpred}{if pred=T, then sent observations for which the rule should be applied}
  \item{lower}{lower bound for predictions}
  \item{upper}{upper bound for predictions} 


}
\value{
  a list containing: 
  \item{cutpoints}{cutpoints corresponding to level thresh}

}

\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
  # see XReg
}
\keyword{ regression }
\keyword{ nonlinear }