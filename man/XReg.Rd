
\name{XReg}
\alias{XReg}
\title{Extreme Regression}
\description{ 
  \code{XReg}             This function estimates parameters in an extreme regression model }
\usage{XReg(y, x, varlist, lower =  - Inf, upper = Inf, 
limsize = (min(c(50,max(c(0.05*length(y),25)))))/length(y), niter = 5, betainit = NULL, 
var.step = F,step.adapt = F, step.size = 1, adaptcntlim = 6, randinit = 0, 
update.terms = as.matrix(c(0, 0)), survival=F,times=y,status=y)    }

\arguments{
  \item{y}{response variable }
  \item{x}{matrix of predictors}
  \item{varlist}{ is a list specifying the form of the extreme regression model. 
   EXAMPLE:  varlist=list(1, c(2,3),c(3,4)) --> max (a1+b1x1,min(a2+b2x2,a3+b3x3),min(a4+b4x3,a5+b5x4)) }
  \item{lower}{ minimum prediction }
  \item{upper}{ maximum prediction 
   EXAMPLE lower=min(y) and upper=max(y) predictions will always be in the range of y. }
  \item{limsize}{the targeted smallest number of observations used to estimate a univariate function.}
  \item{niter}{number of estimations steps}
  \item{var.step}{variable step-size}
  \item{betainit}{ initial estimates for regression parameters. This is in 
        list form the same as output of the regression parameters.}
  \item{step.adapt}{If T the step size is adapted so that the sum of squares is reduced at each step.
       This slows down the algorithm.}
  \item{adaptcntlim}{This limits the number of reduction steps in the adaptive step selection to reduce sums of squares.
        Larger the numbers can slow algorithm.}
  \item{step.size}{Fixed step size. A size of 1 moves directly to the least squares solution. 
            A value of 1 is almost always too large. Better results are obtained with fixed 
            step sizes about .2-.5. The larger the step size the 
            more rapidly the estimates change, but may also may cause lack of convergence.}
  \item{update.terms}{Only update the specified terms. EXAMPLE: Suppose the 
        varlist=list(1,2,3,c(1,2),c(1,3),c(2,3)). One needs to give it an betainit list 
        corresponding to the varlist above  but then to specify only to update terms c(1,2),c(1,3),c(2,3) 
        above one would specify update.terms=cbind(c(4,1),c(4,2),c(5,1),c(5,2),c(6,2),c(6,3)). 
        Update components of the 4th, 5th and 6th terms in the list. }
  \item{survival}{If yes, the expoential survival model is used.}
  \item{times}{If survival=T, time under observation}
  \item{status}{If survival=T, survival status (1=dead, 0=alive)}
  \item{randinit}{If randinit>0 some noise is added to the initial estimates.
    Randinit adds noise of standard deviation form randinit/sqrt(n). This avoids bad initial starts when the 
    same variable is involved in multiple terms. }
}
\value{
  a list containing: 
  \item{beta0list}{list of parameter estimates} 
  \item{pred0list}{list of predictions from each univariate model}
  \item{prop0list}{proportions of observations in each subregion}
  \item{actset}{observations identified in each subregion}
  \item{actset0}{observation used in the estimation of each univariate model}
  \item{predicted}{overall XR model prediction}
  \item{varlist}{varlist as input to the function}
  \item{outmat}{gives the current activeset and change with each iteration. Useful for 
       diagnosing convergence issues.}




}

\author{ Michael LeBlanc \email{mleblanc@fhcrc.org}}

\examples{
##
  

# Simple MAX-MIN  example

ffsimp=function(x) {
	eta=pmax(pmin(1*x[,1],.5*x[,2]),pmin(x[,3],x[,4]))
	return(eta) }

# training data
set.seed(123)
n=500
x=matrix(rnorm(n*5),ncol=5)
eta=ffsimp(x)
y=eta+rnorm(n)*1

# For survival data would need to load SURVIVAL library

# XR regression
# specify form of model. 
varlist=list(c(1,2),c(3,4))
outxr1=XReg(y,x,varlist,step.size=.5,niter=20) 

# XR regression using adaptive step-size

outxr2=XReg(y,x,varlist,step.size=.5,step.adapt=TRUE,niter=20) 


# ----------------------------------------------------------------------------
# XR Survival analysis

# load packages
library(survival)  
library(splines)

# MAKE SURVIVAL DATA
tt=-log(runif(n))/exp(eta-mean(eta))  # make uncensored survival times
cc=runif(n)*4   # make censoring times

status=1*(tt<cc)
times=status*tt+(1-status)*cc                 

# Exponential extreme regression model

outxr4=XReg(y,x,varlist,step.size=.5,step.adapt=FALSE,niter=20,survival=TRUE,times=times,status=status) 

#------------------------------------------------------------------------------------
# Find the rule correspondint to a cutpoint eta=0

ooo=level.set(0,outxr1,pred=TRUE,xpred=x)

print(ooo$cutpoints)
#----------------------------------------------------------------------------------
# EXAMPLES OF PLOTTING FUNCTIONS
# Plot the regression function

plotXRreg(y,x,outxr1,allpoints=FALSE) # only show active points
plotXRreg(y,x,outxr1,allpoints=TRUE) # show all data points

# Plot the inverse regression function

plotXRinvreg(x,outxr1,shading=TRUE,quant=TRUE) # plot against quantiles
# Plot against regression function (but put quantiles on upper horizontal axis
plotXRinvreg(x,outxr1,shading=TRUE,quant=FALSE,both=TRUE,qlev=(1:9)/10)  

#----------------------------------------------------------------------------------
# STEPWISE BUIDING FOR XR MODEL 
# this example builds models on martingale residuals


outstepxr1=stepXR(y,x,ypred=y,xpred=x,niterfull=10,niter=5,step.size=.3,kstep=8,kbest=4,penalty=3,survival=FALSE)

# pick model with best GCV
varlist=outstepxr1$steplist[[outstepxr1$bestgcv]]
# fit final model with larger number of interations.
outxrs=XReg(y=y,x=x,varlist=varlist,step.adapt=TRUE,niter=20,limsize=.05,step.size=.5,
             survival=FALSE)

# Get predicted values from the model 
etaxr=predictfun(xpred=x,varlist=varlist,outxrs$beta0list)$predicted

# get rules for upper 90% of predicted values 
levelout=level.set(quantile(etaxr,.90),outxrs,pred=TRUE,xpred=x)
print(levelout$cutpoints)

# These are the observations corresponding to that group
group=levelout$selectobs
}
\keyword{ regression }
\keyword{ nonlinear }