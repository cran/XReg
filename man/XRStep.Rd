
\name{stepXR}
\alias{stepXR}
\title{Forward stepwise function for adding terms to a extreme regression model}
\description{ 
\code{stepXR} 
This functions calls builds extreme regression models in a step-wise fashion }


\usage{stepXR(y, x, ypred = y, xpred = x, kstep = 5, varlist = NULL, outtemp = NULL, niter = 5, step.size = 0.2, 
	niterfull = 20, step.sizefull = 0.2, kbest = 1, limsize = 0.05, step.adapt = F, upper = max(ypred), lower
	 = min(ypred), predict = T, randinit=0,penalty = 2, update.all = T,survival=F,times=y,status=y)
}


\arguments{
\item{y}{response variable }
  \item{x}{matrix of predictors}
  \item{xpred}{xmatrix for predictions}
  \item{ypred}{y values on test data set}
  \item{varlist}{ is a list specifying the form of the extreme regression model. 
   EXAMPLE:  varlist=list(1, c(2,3),c(3,4)) --> max (a1+b1x1,min(a2+b2x2,a3+b3x3),min(a4+b4x3,a5+b5x4)) }
  \item{lower}{ minimum prediction }
  \item{outtemp}{ current output model}
  \item{step.adapt}{If T the step size is adapted so that the sum of squares is reduced at each step.
       This slows down the algorithm.}
  \item{predict}{ make predictions on the xpred values}
  \item{update.all}{ update all coefficients as terms are added}
  \item{upper}{ maximum prediction 
   EXAMPLE lower=min(y) and upper=max(y) predictions will always be in the range of y. }
  \item{limsize}{the targeted smallest number of observations used to estimate a univariate function.}
  \item{niter}{number of estimations steps}
  \  \item{step.adapt}{If T the step size is adapted so that the sum of squares is reduced at each step.
       This slows down the algorithm.}
       \item{survival}{If yes, the expoential survival model is used.}
  \item{times}{If survival=T, time under observation}
  \item{status}{If survival=T, survival status (1=dead, 0=alive)}
  \item{randinit}{If randinit>0 some noise is added to the initial estimates.
    Randinit adds noise of standard deviation form randinit/sqrt(n). This avoids bad initial starts when the 
    same variable is involved in multiple terms. }    
  \item{kstep}{the number of linear components to add}
  \item{penalty}{penalty used in model selection}
    \item{step.size}{step.size for each term considered for addition}
  \item{niterfull}{number of steps after term has been selected}\
  \item{step.sizefull}{ as above}
  \item{kbest}{used in addeval function and is the number of variables considered for selection. 
   It is the kbest variables with highest univariate correlation.}


}
\value{
  a list containing: 
  \item{stepcoef}{list of coefficients for the stepwise process}
  \item{steplist}{list of variables for the stepwise process}
  \item{stepfit}{fit for training data for each model}
  \item{stepfitpred}{fit for test data (if given) for each model}
  \item{bestgcv}{best model by gcv}
}

\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
  # See XReg
}
\keyword{ regression }
\keyword{ nonlinear }