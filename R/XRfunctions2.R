
#--------------------------------------------------------------------------
# XReg --- THIS IS THE EXTREME REGRESSION FUNCTION FOR A GIVEN FUNCTIONAL FORM
#-------------------------------------------------------------------------


# XReg_function(y, x, varlist, lower =  - Inf, upper = Inf, 
#      limsize = min(c(50,max(c(0.05*length(y),25))))/length(y), 
#      niter = 5, betainit = NULL, var.step = F,	
#      step.adapt = F, step.size = 1, 
#      adaptcntlim = 6, randinit = 0, update.terms = as.matrix(c(0, 0)), 
#      survival=F,times=y,status=y)


# y: response variable
# x: matrix of predictors
# varlist: is a list specifying the form of the extreme regression model. 
# EXAMPLE  varlist=list(1, c(2,3),c(3,4)) --> max (a1+b1x1,min(a2+b2x2,a3+b3x3),min(a4+b4x3,a5+b5x4))
# lower: mimimum prediction (along with upper limits extreme predictions)
# upper: maximum prediction
# EXAMPLE lower=min(y) and upper=max(y) predictions will always be in the range of y.
# limsize: gives the target smallest number of observations used to estimate 
#          a univariate function
# niter:   number of estimations steps

# betaint: one can specifie initial estimates for regression parameters. This is in 
#          list form the same as output of the regression parameters.

# step.adapt=F: If true the step size is adapted so that the sum of squares is reduced at each step.
#           of course this slows down the algorithm
# adaptcntlim: This limits the number of reduction steps in the adaptive step selection to reduce sums of squares.
#            Larger the number can slow algorithm.
# step.size: Fixed step size. A size of 1 moves directly to least squares solution. This is almost always too
#           large. Better results are obtained with fixed step sizes about .2-.5. The larger the step size the 
#           more rapidly the estimates change -- but also may cause lack of convergence.
# update.terms: Only update these terms
# survival: If yes the expoential survival model is used
# times and status: Time under observation and survival status (1=dead, 0=alive
#          This is only needed if survival=Y.

# randinit: If randinit>0 some noise is added to the initial estimates. 
#          Randinit adds noise of standard deviation form randinit/sqrt(n)

# COMMENT: 
# Note if the same variable is involved in multiple model components then initial
# estimates become important (add some noise to the marginal estimates to get started)
# 
 
# COMMENTS on update terms:
# For instance, suppose the varlist=list(1,2,3,c(1,2),c(1,3),c(2,3))
# One needs to give it an betainit list corresponding to the varlist above
# but then to specify only to update terms c(1,2),c(1,3),c(2,3) above one would specify
# update.terms=cbind(c(4,1),c(4,2),c(5,1),c(5,2),c(6,2),c(6,3)). Update components of the
# 4th, 5th and 6th terms in the list. 


# OUTPUT

# beta0list: list of parameter estimates 
# pred0list: list of predictions from each univariate model
# prop0list: proportions of observations in each subregion
# actset: observations identified in each subregion
# actset0: obersvation used in the estimation of each univariate model
# predicted: over XR model prediction
# varlist: varlist as input to the function -- list of terms
# outmat: gives the current activeset and change with each iteration
#         useful for diagnosing convergence issues
#------------------------------------------------------------------------------------------
# Additional Code needed for using survival data


survcalc<-function(tt,ss,xx,fit=T) {
if (fit) {
aa<-survreg(Surv(tt,ss)~xx,dist="exponential",control=survreg.control(maxiter=10))
err=-2*aa$loglik[2]
coef<-aa$coef   }
else {
aa=survreg(Surv(tt,ss)~xx,init=c(0,1),dist="exponential",control=survreg.control(maxiter=0))
err=-2*aa$loglik[2]
coef=aa$coef   }
survout=list(err=err,coef=coef) 
return(survout) }


# MAIN EXTREME REGRESSION FUNCTION

XReg=function(y, x, varlist, lower =  - Inf, upper = Inf, limsize = (min(c(50,max(c(0.05*length(y),25)))))/length(y), niter = 5, betainit = NULL, var.step = F,
	step.adapt = F, step.size = 1, adaptcntlim = 6, randinit = 0, update.terms = as.matrix(c(0, 0)),
      survival=F,times=y,status=y)
{
	# some set up calculations
	n <- length(y)
	mterm <- length(varlist)
	outmat <- matrix(NA, nrow = niter, ncol = 2)
	# Set up the default active set 
	actset <- actsetn <- list(NULL)
	for(j in 1:mterm) {
		varvec <- as.vector(varlist[[j]])
		nterm <- length(varvec)
		actset[[j]] <- matrix(T, nrow = n, ncol = nterm)
	}
	actsetn <- actset
	#loop
	if(!is.null(betainit)) {
		beta0list <- betainit
	}
	else {
		beta0list <- list(NULL)
	}
	for(kk in 1:niter) {
		actsetchange <- mean(abs(unlist(actset) - unlist(actsetn)))
		actset <- actsetn
		# calculate the coefficients and predicted on active set
		# use initial estimates if supplied
		pred0lista <- pred0list <- list(NULL)
		prop0lista <- prop0list <- list(NULL)
		outerr0lista <- outerr0list <- list(NULL)
		beta0lista <- beta0list
		min0mat <- matrix(NA, nrow = n, ncol = mterm)
		adaptpred <- Inf
		adaptcnt <- 1
		if(kk == 1) {
			priorpred <- Inf
		}
		else {
			priorpred <- outmat[kk - 1, 1]
		}
		# here we loop at least once
		# we loop more if adaptive step size was used to make sure the RSS always goes down.
		while((adaptcnt == 1) | ((kk > 1) & (step.adapt) & (adaptcnt < adaptcntlim) & (adaptpred > 
			priorpred))) {
			for(j in 1:mterm) {
				#Loop over the outterms
				varvec <- as.vector(varlist[[j]])
				nterm <- length(varvec)
				beta0 <- matrix(NA, nrow = 2, ncol = nterm)
				pred0 <- matrix(NA, nrow = n, ncol = nterm)
				outerr0 <- prop0 <- rep(NA, nterm)
				actset0 <- as.matrix(actset[[j]])
				for(i in 1:nterm) {
					# Loop on the inner set
					# on first loop use inital estimates if supplied
					if(kk == 1) {
						if(!is.null(betainit)) {
							beta0[, i] <- betainit[[j]][, i]
						}
						else {
                                   if (!survival) {
							beta0[, i] <- fastreg(x[actset0[, i], varvec[i]], y[actset0[
								, i]])$coef + (rnorm(2) * randinit)/sqrt(n)}
                                    # survival outcome
                                    if (survival) {
                                    survout=survcalc(tt=times[actset0[, i]],ss=status[actset0[, i]],xx=x[actset0[, i], varvec[i]] )
                                    beta0[, i] <- -survout$coef + (rnorm(2) * randinit)/sqrt(n)
                                                                        
                                   }
						}
					}
					# fit on the active 
					# here we have the step.size incorporated into the updating algorithm
					if(kk > 1) {
						#xvec=x[actset0[,i],varvec[i]]
						if(length(actset0[, i]) > 2) {
							m1 <- update.terms[1,  ] == j
							m2 <- update.terms[2,  ] == i
							if((sum(update.terms) == 0) | (sum(m1 * m2) > 0)) {
								if(var.step)
									fact <- 1 - exp( - step.size * kk)
								else fact <- step.size
								# here is where we do the adaptive step size 
								fact <- step.size * 0.75^(adaptcnt - 1)
                               if (!survival) {
					 beta0[, i] <- beta0list[[j]][, i] + fact * (fastreg(
									x[actset0[, i], varvec[i]], y[actset0[
									, i]])$coef - beta0list[[j]][, i]) }
                                   if (survival) {
                                    acts=unlist(actset0[,i])
                                    survout=survcalc(tt=times[acts],ss=status[acts],xx=x[acts, varvec[i]] )
                                    beta0[, i] <- beta0list[[j]][, i] + fact * (-survout$coef-beta0list[[j]][, i])
                                                                          
                                     }
					

							}
							else {
								beta0[, i] <- beta0list[[j]][, i]
							}
						}
					}
					# fit on the active set
					pred0[, i] <- as.vector(cbind(1, x[, varvec[i]]) %*% beta0[, i])
					prop0[i] <- mean(actset0[, i])
					outerr0[i] <- sum((y[actset0[, i]] - pred0[actset0[, i], i])^2)
                              #  survival outcome

                              if (survival) {
                              acts=unlist(actset0[,i])
                              survout=survcalc(tt=times[acts],ss=status[acts], xx=-pred0[acts, i],fit=F)
                              outerr0[i]<-survout$err
                              
                                     }

				}
				# Calculate minimum of all elements in term
				min0mat[, j] <- apply(cbind(pred0, upper), 1, "min")
				# this is the update of all coefficients. So it should be on the outside of the adapt step.
				beta0lista[[j]] <- beta0
				pred0lista[[j]] <- pred0
				prop0lista[[j]] <- prop0
				outerr0lista[[j]] <- outerr0
			}
			# Calculate maximum
			maxvec <- apply(cbind(min0mat, lower), 1, "max")
			adaptpred <- mean((y - maxvec)^2)
                  # outcome
                  if (survival) {
                  survout=survcalc(times,status, maxvec,fit=T)
                  adaptpred<-survout$err                                                                   
                                     }
 
			if(step.adapt)
				cat("STEP= ", adaptcnt, "STEPSIZE= ", round(step.size * 0.75^(adaptcnt - 1), 3),
					"RESID ERROR= ", adaptpred, "\n")
			adaptcnt <- adaptcnt + 1
		}
		beta0list <- beta0lista
		pred0list <- pred0lista
		prop0list <- prop0lista
		outerr0list <- outerr0lista
            curerr=mean((y - maxvec)^2)
            if (survival) {   survout=survcalc(tt=times,ss=status, xx=maxvec,fit=T)
                              curerr <-survout$err     }
           print( c(curerr, actsetchange))
		outmat[kk,  ] <- c(curerr,actsetchange)
		# Calculate the new active set
		actsetn <- actset
		for(j in 1:mterm) {
			varvec <- as.vector(varlist[[j]])
			nterm <- length(varvec)
			for(i in 1:nterm) {
				discrep <- abs(maxvec - pred0list[[j]][, i])
				actsetn[[j]][, i] <- discrep == 0
				if(is.na(actsetn[[j]][1, i]))
					browser()
				if(mean(actsetn[[j]][, i]) < limsize) {
					actsetn[[j]][, i] <- actsetn[[j]][, i] | (discrep <= (quantile(discrep,
						limsize) + 1e-005))
				}
			}
		}
	}
	predicted <- maxvec
	# calculate real active set
	actset0 <- actset
	for(j in 1:mterm) {
		varvec <- as.vector(varlist[[j]])
		nterm <- length(varvec)
		for(i in 1:nterm) {
			discrep <- abs(maxvec - pred0list[[j]][, i])
			actset0[[j]][, i] <- discrep == 0
			prop0list[[j]][i] <- mean(actset0[[j]][, i])
		}
	}
        fit=list(beta0list=beta0list,pred0list=pred0list,prop0list=prop0list,actset=actset,
                 actset0=actset0,predicted=predicted,varlist=varlist,outmat=outmat)
	return(fit)
}


#-------------------------------------------------------------------------
# FASTREG:  This function avoids matrix inversion on all the univariate regresssions
#-------------------------------------------------------------------------
# It doesnt add much speed in the current implementation. However, it would be easy to
# Modify for penalized versions if we try that. 

fastreg=function(x,y,wt=rep(1,length(y)),delta=.00001) {
x=as.vector(x)
y=as.vector(y)
wt=as.vector(wt)
coef=rep(NA,2)
mx=mean(wt*x)
sxy=mean(wt*y*(x-mx))
sxx=mean((x-mx)*(x-mx))
coef[2]=sxy/(sxx+delta)
coef[1]=mean(y*wt)-coef[2]*mx
return(list(coef=coef)) }



# ----------------------------------------------------------------
# LEVEL SET FUNCTION
# ----------------------------------------------------------------
# This function gives the decision rule for a given threshold.
# It also tells how many observations are indicated by each rule
# and give the overall indicator of the observations indicated by
# the rule. Currently this function is only written for >= rules, it
# one should be able to specify the desired inequality

# INPUT
# thresh: threshold or target output for rule {x:f(x)>=thresh}
# output: output from XR regression function
# pred=F: determine which observations are consistent with rule
# xpred: if pred=T, then sent observations for which the rule should be applied
# lower and upper: lower and upper bounds for predictions 

level.set=function(thresh, output,pred=F,xpred=NULL,lower=-Inf,upper=Inf) {

beta0list=output$beta0list
varlist=output$varlist
mterm=length(varlist)
cutpoints=cutpointout=cutpointdir=list(NULL)
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
nterm=length(varvec)
cutpointvec=varcutvec=grltvec=varcutdir=rep(NA,nterm)
for (i in 1:nterm) {
cutpt=(thresh-beta0list[[j]][1,i])/beta0list[[j]][2,i]
signb=(beta0list[[j]][2,i])>0
cutpointvec[i]=cutpt
varcutvec[i]=varvec[i]
varcutdir[i]=signb
grltvec[i]=c(" Less   "," Greater ")[signb+1]
 }
# put in list
cutpoints[[j]]=cbind(varcutvec,grltvec,round(cutpointvec,4))
cutpointout[[j]]=cutpointvec
cutpointdir[[j]]=varcutdir
}

if (pred) {
bb=predictfun(xpred,output$varlist,output$beta0list,lower,upper,activeset=T)
predicted=bb$predicted
pred0list=bb$newpred0list
mterm=length(varlist)
nlist=list(NULL)
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
nterm=length(varvec)
nvec=rep(NA,nterm)
for (i in 1:nterm) {
nvec[i]=sum((pred0list[[j]][,i]==predicted)&(predicted>=thresh)) }
# calculate the number of observations associated with each component
nlist[[j]]=nvec }
selectobs=predicted>=thresh }
out=list(selectobs=selectobs,cutpointdir=cutpointdir,cutpointout=cutpointout,cutpoints=cutpoints,nlist=nlist)
return(out) }

 
# ----------------------------------------------------------------
# PREDICTION FUNCTION 
# ----------------------------------------------------------------
# This gives predictions but also give assigments to the active sets 
# for data given in the xpred matrix. The default is not to give assignments.
# INPUT
# xpred: xmatrix for predictions
# varlist: model form list (described in XR regression function)
# beta0list: the beta list from XR regression function
# lower,upper: lower and upper bounds for predictions (as described in XR regression function)
# active set: if T denote where observations are in terms of active set.

predictfun=function(xpred,varlist,beta0list,lower=-Inf,upper=Inf,activeset=F) {
mterm=length(varlist)
n=dim(xpred)[1]
newpred0list=list(NULL)
minnewmat=matrix(NA,nrow=n,ncol=mterm)
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
nterm=length(varvec)
prednew=matrix(NA,nrow=n,ncol=nterm)
for (i in 1:nterm) {
prednew[,i]=as.vector(cbind(1,xpred[,varvec[i]])%*%beta0list[[j]][,i]) }
# Calculate minimum of all elements in term
newpred0list[[j]]=prednew
minnewmat[,j]=apply(cbind(prednew,upper),1,"min")
}

# Calculate maximum (predicted)
predicted=apply(cbind(minnewmat,lower),1,"max")


# Here we calculate the active set membership and numbers for each of the predicted cases 

actset=list(NULL)
nlist=list(NULL)
# Default is not to calculate active set -- this is to cut down on computation
# if this function is called many times in model building
if (activeset) {
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
nterm=length(varvec)
nlistv=rep(NA,nterm)
actset[[j]]=matrix(NA,nrow=n,ncol=nterm)
for (i in 1:nterm) {
discrep=predicted-newpred0list[[j]][,i]
actset[[j]][,i]=discrep==0
nlistv[i]=sum(discrep==0) }
nlist[[j]]=nlistv } }
predfit=list(predicted=predicted,actset=actset,nlist=nlist,newpred0list=newpred0list)
return(predfit) }



# ------------------------------------------------------------------
# CALCULATE SUM OF SQUARES FOR DROPPING 1-TERM  
#-----------------------------------------------------------------
# This function drops each term in the model and calculates the drop
# in sum of squares associated with each term removal. 
# One can refit or leave parameter estimates as in full model fit.
# This function works on training or test sample data if available 
# The most computationally efficient is to have niter=0



dropeval=function(y,x,ypred=y,xpred=x,varlist,outtemp,niter=5,step.size=.5,
            upper=max(ypred),lower=min(ypred),randinit=.1) {

beta0list=outtemp$beta0list
predfull=predictfun(xpred,varlist,beta0list)$predicted
# OUTCOME
fullfit=mean((ypred-predfull)**2)
dropfit=list(NULL)
mterm=length(varlist)
for(j in 1:mterm) {
varvec=as.vector(varlist[[j]])
nterm=length(varvec)
dropvec=rep(NA,nterm)
for (i in 1:nterm) {
beta1list=beta0list
if (nterm==1) {
vv=varlist[-j]
beta1list=beta0list[-j] } else { vv=varlist
                       vv[j]=vv[[j]][-i] 
                       beta1list[[j]]=as.matrix(beta0list[[j]][,-i]) }
print(vv)
print(beta1list)
# could add a refit component in here
if (niter==0) beta0=beta0list else {

beta0=XReg(y,x,varlist=vv,lower=lower,upper=upper,betainit=beta1list,niter=niter,step.size=step.size,randinit=randinit)$beta0list }
yp=predictfun(xpred,vv,beta0,lower=lower,upper=upper)
# OUTCOME
dropvec[i]=mean((ypred-yp$predicted)**2)-fullfit } 
dropfit[[j]]=dropvec }

evalfit=list(fullfit=fullfit,dropfit=dropfit)

return(evalfit) }


# -----------------------------------------------------------------------
# XR BUILD FUNCTION
#------------------------------------------------------------------------

# This function builds a XR model by adding terms in a computationally 
# way. The function starts with a varlist of interest and adds 1 new term.
# kbest is the number of variables to be considered for adding -- kbest variables
# are chosen by looking at univariate correlations to variables considered for 
# inclusion in the model.

addeval=function(y,x,varlist=NULL,outtemp=NULL,niter=5,step.size=.2,kbest=1,   
            limsize=.05,upper=max(y),lower=min(y),update.all=T,repeat.term=F,randinit=0) {
 if ( !is.null(varlist) ) {
    beta0list=outtemp$beta0list
    predcurr=predictfun(x,varlist,beta0list)$predicted 
    nterm=length(varlist) }
  else {   predcurr=mean(y)
    nterm=0 }
# FIT FOR THE CURRENT MODEL
  resid=y-predcurr
  currfit=mean(resid**2)
  p=dim(x)[2]
  correlat=cor(resid,x)**2
# PICK THE k-best VARIABLES OF INTEREST 
  newvarpos=(1:p)[order(-correlat)<=kbest]
  # If this is the first term -- just a univariate linear model so we are done
  if (nterm==0) {
  newvar=(1:p)[order(-correlat)<=1]
  varlist.sel=list(newvar)
  bestfit=mean( (y-predict(lm(y~x[,newvar])))**2)  }
  else {
  nterm=length(varlist)
# If there are already terms in model it becomes more complicated.

# TWO LOOPS 
#FIRST OVER VARIABLES

bestfit=Inf 

for (kk in 1:kbest) {
newvar=newvarpos[kk]    
addfit=rep(Inf,nterm+1)
# add variable as singleton
      
varlists=list(NULL)
varlist1=varlist
# don't add the same single variable
t10=unlist(lapply(lapply(lapply(varlist,match,newvar),"is.na"),sum))==0   #no missing in the match
t20=unlist(lapply(lapply(varlist,match,newvar),"length"))==1            #full match with new term
if ((sum(t10*t20)<1)|(repeat.term)) {
varlist1[[nterm+1]]=newvar 
# XReg regression with potential new term
# check and see if full updating was wanted
if (update.all) { updterm=as.matrix(c(0,0)) } else { updterm=as.matrix(c(nterm+1,1)) }
beta0=XReg(y,x,varlist=varlist1,niter=niter,step.size=step.size,limsize=limsize,lower=lower,upper=upper,update.terms=updterm,randinit=randinit)$beta0list 
yp=predictfun(x,varlist=varlist1,beta0,lower=lower,upper=upper)
# OUTCOME

addfit[1]=mean((y-yp$predicted)**2)  }
else {
  addfit[1]=Inf }
varlists[[1]]=varlist1

# add as combination (to any of existing terms) 


for (i in 1:nterm) {
vv1=varlist
vv0=c(varlist[[i]],newvar)          # potential new combination variable
# don't add the same term as is currently in the model
t1=unlist(lapply(lapply(lapply(varlist,match,vv0),"is.na"),sum))==0        #no missing in the match
t2=unlist(lapply(lapply(varlist,match,vv0),"length"))==length(vv0)    #full match with new term


# Don't add same variable in the same subterm or same term into the model

if ((!is.na(match(newvar,varlist[[i]])))||(((sum(t1*t2)>0))&(!repeat.term))) {
#if (!is.na(match(newvar,varlist[[i]]))) {

    addfit[1+i]=Inf  } else {
vv1[[i]]=c(varlist[[i]],newvar)      # this just replaces current term with expanded term
varlists[[i+1]]=vv1

# XR fun regression with potential new term
# check and see if full updating was wanted
subterm=length(varlist[[i]])+1
if (update.all) { updterm=as.matrix(c(0,0)) } else { updterm=as.matrix(c(i,subterm)) }

beta0=XReg(y,x,varlist=vv1,niter=niter,step.size=step.size,limsize=limsize,lower=lower,upper=upper,update.terms=updterm,randinit=randinit)$beta0list 
yp=predictfun(x,varlist=vv1,beta0,lower=lower,upper=upper)
addfit[1+i]=mean((y-yp$predicted)**2) 
}
 }
# this is a switch to eliminate add terms
select.terma=((1:(nterm+1))[addfit==min(addfit)])[1] 
varlist.sela=varlists[[select.terma]]
bestfita=addfit[select.terma] 
  
if (bestfita<bestfit) { bestfit=bestfita
                        varlist.sel=varlist.sela
                        select.term=select.terma } 
               } }

if ( bestfit==Inf) { print("no variables can be added")
                     varlist.sel=varlist }
                   
buildfit=list(varlist.sel=varlist.sel,bestfit=bestfit)
return(buildfit) }

#----------------------------------------------------------------------------------------
# stepXR  forward stepwise function for adding terms to a logic model
#---------------------------------------------------------------------
# kstep is number of linear components to add
# ninter, step.size (this refers to each term considered for addition)
# niterfull, step.sizefull (this refers to fitting after term is selected)
# kbest is used in addeval - and is the number of variables considered for selection
# it is the kbest variables with highest univariate correlation.



stepXR=function(y, x, ypred = y, xpred = x, kstep = 5, varlist = NULL, outtemp = NULL, niter = 5, step.size = 0.2, 
	niterfull = 20, step.sizefull = 0.2, kbest = 1, limsize = 0.05, step.adapt = F, upper = max(ypred), lower
	 = min(ypred), predict = T, randinit=0,penalty = 2, update.all = T,survival=F,times=y,status=y)
{
	steplist <- list(NULL)
	stepfitpred <- stepfit <- rep(NA, kstep)
	stepcoef <- list(NULL)
	varlist1 <- varlist
	n <- length(y)
	for(jj in 1:kstep) {
		varlist1 <- addeval(y, x, varlist = varlist1, outtemp = outtemp, kbest = kbest, limsize = limsize,
			lower = lower, upper = upper, update.all = update.all,randinit=randinit)$varlist.sel
		outtemp <- XReg(y, x, varlist = varlist1, niter = niterfull, step.adapt = step.adapt, step.size = 
			step.sizefull, lower = lower, upper = upper,randinit=randinit)
		# switch to parameters
		steplist[[jj]] <- varlist1
		stepfit[jj] <- mean((y - outtemp$predict)^2)
             if (survival) {
                survout=survcalc(tt=times,ss=status, xx=outtemp$predict,fit=F)
                stepfit[jj]=survout$err/n           }
		stepcoef[[jj]] <- outtemp$beta0list
		if(predict) {
			outpred <- predictfun(xpred, varlist = varlist1, outtemp$beta0list, lower = lower, upper = 
				upper)$predicted
            
			stepfitpred[jj] <- mean((ypred - outpred)^2)
               
         }
		
	}
	gcv <- stepfit/(1 - (penalty * (1:kstep))/n)
	bestgcv <- min((1:kstep)[min(gcv, na.rm = T) == gcv], na.rm = T)
        stepfit=list(stepcoef=stepcoef, steplist=steplist, stepfit=stepfit, stepfitpred=stepfitpred, 
        bestgcv=bestgcv)
	return(stepfit)
}

#--------------------------------------------------------------------------
# CVXR - this does k-fold crossvalidation of the stepwise building process.
#--------------------------------------------------------------------------

# fold: number of folds in cross-validation

CVXR=function(fold=5,y,x,kstep=5,varlist=NULL,outtemp=NULL,niter=5,niterfull=15,step.size=.2,
                step.sizefull=.2,kbest=1,
                upper=max(y),lower=min(y))  {
bbfull=stepXR(y,x,ypred=y,xpred=x,kstep=kstep,varlist=varlist,outtemp=outtemp,
              niter=niter,step.size=step.size,
              niterfull=niterfull,step.sizefull=step.sizefull,kbest=kbest,
              upper=upper,lower=lower,predict=T)
bbcv=list(NULL)
n=length(y)
predmat=matrix(NA,ncol=length(bbfull$stepfit),nrow=fold)
rand <- sample(fold,n, replace = T)
for (kk in 1:fold) {
xtrain=x[rand!=kk,]
ytrain=y[rand!=kk]
xtest=x[rand==kk,]
ytest=y[rand==kk]
print(kk)
bbcv0=stepXR(ytrain,xtrain,ypred=ytest,xpred=xtest,kstep=kstep,varlist=varlist,outtemp=outtemp,
                niter=niter,step.size=step.size,
                niterfull=niterfull,step.sizefull=step.sizefull,kbest=kbest,
                upper=upper,lower=lower,predict=T)
predmat[kk,]=bbcv0$stepfitpred
bbcv[[kk]]=bbcv0 }

cvfit=list(bbcv=bbcv,bbfull=bbfull,predmat=predmat)
return(cvfit) }



