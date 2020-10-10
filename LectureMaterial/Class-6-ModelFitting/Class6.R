###################################################################
###################################################################
# Contrasting the effects of assuming observation vs. process error
#             and
# Model fitting and calculation of AIC
###################################################################
rm(list=ls()) # clears workspace
par(mfrow=c(1,4))

###################################################################
# Observation Error only
########################
# Load data for South African fur seals ('safseal1.csv') or Grizzlies ('grizzlies.csv')
dat<-read.csv(file.choose())
# or use: dat<-read.csv('~/Desktop/safseal1.csv')

Year<-dat[,1]; Nobs<-dat[,2];
T = length(Year);

plot(Nobs~Year)
plot(log(Nobs)~Year)

# Estimate lambda and observation error: least-squares fit of log(Nt) vs. time
fit<-lm(log(Nobs)~Year)
abline(fit)

# lambda= e^r, thus corresponds to the exponentiatd 'b' coefficient of 'y=a+b*x'
r_obs<-coef(fit)[2]

# calculate observer error as sum of squared residuals
resids<-resid(fit)
v_est = sum(resids^2)/(T-1) # corrected for finite sample size
# convert to standard deviation
sd_obs = sqrt(v_est)

# Display estimates
r_obs
sd_obs

####################
# Process Error only
####################

# estimate log of year-to-year lambda values from data (assumes constant step size)
log_lambda_obs<-log(Nobs[-1]/Nobs[-T])

plot(Nobs[-T],log_lambda_obs,ylab='log(Nt+1/Nt)',xlab='N')

# fit
fit1<-lm(log_lambda_obs~Nobs[-T])
summary(fit1)
# Note that slope (i.e. density-dependence) is not deemed "significant"

# Hence estimate r by mean
r_proc<-mean(log_lambda_obs)
# ...or with zero-slope model
fit2<-lm(log_lambda_obs~1)
summary(fit2)

# Estimate process error
v_proc<-var(log_lambda_obs)
sd_proc<-sqrt(v_proc)

# Compare estimates
r_obs; r_proc
sd_obs; sd_proc

abline(h=0,lty=1,col='grey',lwd=2)
abline(r_obs,0,lty=2,col='red',lwd=2)
abline(h=r_proc,lty=2,col='blue',lwd=2)


#################################################
# Simulate future population growth for 100 years 
# with parameters assuming process error
# and compare to
# with parameters assuming observation error
#################################################
library(sfsmisc) # for eaxis function - log-scale tick-marks

reps<-100
FutureTime<-100 # Number of time-steps into future
Np<-array(NA,dim=c(reps,FutureTime+1))
Np[,1]<-Nobs[T] # Starting population size is last observed population size
No<-Np

for (i in 1:reps){
    for (t in 1:FutureTime){
        Np[i,t+1]<- Np[i,t]*exp(r_proc+rnorm(1,0,sd_proc))
        No[i,t+1]<- No[i,t]*exp(r_obs+rnorm(1,0,sd_obs))
}}


plot(Nobs~Year,xlim=c(min(Year),max(Year)+FutureTime), ylim=c(min(Nobs),max(Np)),type='l',log='y',axes=FALSE)

xt<-seq(max(Year)+1,max(Year)+(t+1),1)
axis(1); eaxis(2); box(lwd=1);
for (i in 1:reps){
    lines(xt,Np[i,])
}
for (i in 1:reps){
  lines(xt,No[i,],col='red')
}

#######################################################################
#######################################################################
#######################################################################

