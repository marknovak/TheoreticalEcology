#######################################################################
#######################################################################
# Maximum likelihood - Poisson example
#######################################################################
#######################################################################
rm(list=ls()) # clears workspace
# Data
Y<-4

# The probability function for the Poisson
dPois<-function(y,lambda){
  (exp(-lambda)*lambda^y)/factorial(y)
  }

# lambda's to search over
lambda<-seq(0,20,0.5)

# Probability
Prob<-dPois(Y,lambda)
# quartz(width=4,height=4)
plot(lambda,Prob,ylab='Probability',xlab='Possible lambdas',cex=0.5)

# Compare to built-in function "dpois"
Prob2<-dpois(Y,lambda)
lines(lambda,Prob2,col='grey')

#########
Y<-c(4,6)
JointProb<-dpois(Y[1],lambda)*dpois(Y[2],lambda)
lines(lambda,JointProb)


#######################################################################
#######################################################################
# Model fitting - Exponential & Ricker model
###########################################
# NOTE:
# We're going to use nls() rather than lm() so that we can accomodate 
# the (nonlinear) Theta logistic model of Problem set #3
#######################################################################
#######################################################################
rm(list=ls()) # clears workspace

# Load Great tit data (from last class)
dat<-read.csv(file.choose())

Years<-dat[,1]
Nobs<-dat[,2]
T<-length(Years)

log_lam<-log(Nobs[-1]/Nobs[-T])

# quartz(height=4,width=8)
par(mfrow=c(1,2),cex=0.8,mgp=c(2,0.2,0),tcl=-0.1,las=1)
plot(Years,Nobs,type='b')
plot(Nobs[-T],log_lam,
     xlab=expression(N[t]),
     ylab=expression(log(N[t+1]/N[t])))

# Shorten variable names
x<-Nobs[-T]
y<-as.vector(log_lam)
n<-length(y)

# Fit density-independent model
# (using nonlinear least squares for kicks)
DIfit<-nls(y~b1+0*x, 
           start=list(b1=2))  # nls requires initial guess
summary(DIfit)

# Fit Ricker model using nonlinear least squares. 
# (Note that data is log(N_{t+1}/N_{t}).)
Rfit<-nls(y~b1*(1-x/b2), 
          start=list(b1=2,b2=1))
summary(Rfit)
p_nls<-coef(Rfit)
names(p_nls)<-c('r','K')
print(p_nls)

# Define Ricker function of expected growth for plotting
RickerCurve<-function(x){
    p_nls[1]*(1-x/p_nls[2])
  }
abline(h=0,lty=3)
curve(RickerCurve,
      c(min(x):max(x)),
      add=TRUE,lty=2)

logLik(DIfit)
logLik(Rfit)
AIC(DIfit,Rfit) # => Ricker model performs better

library(bbmle)
AICtab(DIfit,Rfit,
       weights=TRUE)
################################################
# Repeat by optimizing own logLikelihood function
#################################################
# Define process model of expected growth as function 
# of parameters and data
Ricker<-function(x,p){
    p[1]*(1-x/p[2])
  }

# Create function to calculate negative log-likelihood.
# Assume that residuals expected from log(lambda) are normally distributed.  This should give us the same answer as using nonlinear least squares as above.  Why?
nLL<-function(p){ 
    (n/2)*log(2*pi)+
    (n/2)*log(var(log_lam))+
    (1/(2*var(log_lam)))*sum((y-Ricker(x,p))^2)
  }


# Set initial guesses
p=c(r=1,K=50)

# optimize parameters to minimize nLL function
fit<-optim(p,nLL)

fit.nLL<-fit$value
fit.p<-fit$par

myAIC<-function(fit.nLL,fit.p){
  2*(length(fit.p))+2*fit.nLL
  }

# Compare estimates
print(fit.p)
print(p_nls)

# Compare AIC estimates
myAIC(fit.nLL,length(fit.p))
AIC(Rfit)

#######################################################################
#######################################################################
#######################################################################