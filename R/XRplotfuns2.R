
#-------------------------------------------------------------------
# PLOT THE EXTREME REGRESSION FUNCTION 
# Returns a plot of the regression function and data points as a matrix plot
# Rows of matrix are variables and columns correspond to miniumum terms.
# This representation of the regression function probably only works well
# as a graphical display up to about 5 unique variables and 5 terms.
# INPUT
# y: response
# x: matrix of predictors
# xrout: output from XR regression function XReg
# allpoints:   if allpoints=F only plot points in panel associated with that local 
# univarite function. If allpoints=T plot all points (but highlight those corresponding 
# to the local univariate function
#----------------------------------------------------------------------


plotXRreg=function(y,x,xrout,allpoints=F) {
lab.names=names(x)
p=dim(x)[2]
if (is.null(lab.names)) {
  lab.names=paste("x",1:p,sep="") }

varlist=xrout$varlist
mterm=length(varlist)
varall=NULL
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
varall=unique(c(varvec,varall))
 }
varall=sort(varall)
nvar=length(varall)

par(mfrow=c(nvar,mterm))
par(mar=c(4,4,0,0))

y2=c(min(y),max(y))
for (i in varall) {
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
varind=match(i,varvec)
print(varind)    

if (is.na(varind)) { 
  plot(1,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
     }
else {
     # beta0list here could be looped over some simulation list
     x1a=x[,varvec[varind]]
     x1o=order(x1a)
     x1=x1a[x1o]
     x2=c(min(x1),max(x1))
     pred=xrout$pred0list[[j]][,varind]
     plot(x1,pred[x1o],xlim=x2,ylim=y2,ylab="predicted",xlab=lab.names[varvec[varind]],type="l") 

     indpt=xrout$actset[[j]][,varind]
     if (allpoints) { points(x1a,y,col=5) }
     points(x1a[indpt],y[indpt]) } }}
return() }


#-------------------------------------------------------------------
# PLOT THE INVERSE REGRESSION FUNCTION 
# Returns a plot of the inverse regression function, as a matrix plot
# Rows of matrix are variables and columns correspond to miniumum terms.
# This representation of the regression function probably only works well
# as a graphical display up to about 5 unique variables and 5 terms.
# INPUT
# x: matrix of predictors
# xrout: output from XR regression function XReg
# shading:   if shading=T shade the directions of the rules 
# in panel associated with sequence of rules {x:f(x)>=q}  
# quant=T: plot rules in terms of quantiles or ranked predictions 
# Note, models will not appear linear univarite.
# both=T: plot in linear form, but put quantiles on the top of plot 
# qlev: If both=T, these are the labels for the quantiles put on the top of the plot
# spec.ylim=F: do not specify the limits of the y axis of the plot. Note, in this case the same variable 
# different plots could have different axes. If spec.ylim=T, then the user must provide the limits 
# y-axis for each of the variables.
# mat.ylim: This gives the y-axis limits (min,max) for each of the variables, it is a 2xp matrix, where 
# each column corresponds to a predictor variable. It is only used if spec.ylim=T.
# level.set: If this is set to a numeric value, a highlighted line will be added to the plot to indicate 
#----------------------------------------------------------------------


plotXRinvreg=function(x,xrout,shading=F,quant=F,both=F,qlev=c(.2,.4,.6,.8),
                      spec.ylim=F,mat.ylim=NULL,level.set=NA) {
lab.names=dimnames(x)[[2]]
p=dim(x)[2]
if (is.null(lab.names)) {
  lab.names=paste("x",1:p,sep="") }

varlist=xrout$varlist
mterm=length(varlist)
varall=NULL
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
varall=unique(c(varvec,varall))
 }
varall=sort(varall)
nvar=length(varall)
ff=xrout$predicted 
f2=c(min(ff),max(ff))

par(mfrow=c(nvar,mterm))
par(mar=c(4,4,2,0))

for (i in varall) {
for (j in 1:mterm) {
varvec=as.vector(varlist[[j]])
varind=match(i,varvec)

if (is.na(varind)) { 
  plot(1,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
     }
else {
     # beta0list here could be looped over some simulation list
     ff=xrout$predicted 
     f2=c(min(ff),max(ff))
     x1a=x[,varvec[varind]]
     x1o=order(x1a)
     x1=x1a[x1o]
     x2=c(min(x1),max(x1))
     bv=xrout$beta0list[[j]][,varind]
     off=order(ff)
     gg=(ff-bv[1])/bv[2]
     g2=(f2-bv[1])/bv[2]
     if(!is.na(level.set)) {  
     level.setff1=(quantile(ff,level.set)-bv[1])/bv[2]
     level.setff=(level.set-bv[1])/bv[2] }

     rankff=rank(ff)/length(ff)
     if (quant) { ff=rankff
     f2=c(min(ff),max(ff))
     if(!is.na(level.set)) {  
     level.setff=level.setff1 } }
     if (!spec.ylim) {
     plot(ff[off],gg[off],type="l",xlab="predicted",ylab=lab.names[varvec[varind]]) } else {
     plot(ff[off],gg[off],type="l",xlab="predicted",ylab=lab.names[varvec[varind]],ylim=mat.ylim[,varvec[varind]]) }
     if (both) {
            axis(side=3,at=quantile(ff,qlev),label=qlev) }
    if (shading) {
      if (bv[2]<0) {
            polygon(c(rev(ff[off]),f2[c(1,2)]),c(rev(gg[off]),g2[c(2,2)]),density=20,angle=45,lty=2)
            if(!is.na(level.set)) {  
            lines(c(level.set,level.set), c(level.setff,g2[1]) ) 
                                         }
    } # greater
   else {
            polygon(c(ff[off],f2[c(1,1)]),c(gg[off],g2[c(2,1)]),density=20,angle=45,lty=2) 
          if(!is.na(level.set)) {  
            lines(c(level.set,level.set), c(g2[2],level.setff) )
                                         }




}  # lesser

          } } } }
    return() }

