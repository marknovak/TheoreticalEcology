##################################################################
######### Theoretical Ecology - Problem set 3a ###################
#########    Observation vs. Process Error     ###################
##################################################################
rm(list=ls()) # clears workspace
#################################################################
#################################################################
# Define parameters
reps<-10 # Number of simulations to run
T<-100 # Number of time points
N0<-10 # Starting population size
lambda<-1.01 # instantaneous growth rate
sd_p <- 0.1 # process error (standard dev. of log-lambda across time-steps
sd_o <- 0.1 # observation error (standard dev. o estimated log popn sizes)

# Define local functions for use in plotting and fitting
# r<-log(lambda) # per capita growth rate
f1<-function(x){log(N0)+log(lambda)*(x)} # exponential popn growth (lambda known)
f2<-function(x){log(lambda)+0*x} # log(Nt+1/Nt) vs. Nt, exponential growth (lambda known)
f3<-function(b,x){log(b[1])+log(b[2])*x} # log(Nt) vs. t, exponential growth, solve for N0 and lambda

print(paste('Simulated population growth with either proceess or observation error N(0)= ',N0,', lambda = ', lambda,', process error = ',sd_p,' and observation error = ',sd_o))


# Initialize variables to store population counts and growth rates
Time<-array(0,dim=c(reps,T))
Time[,1]<-1
N<-array(0,dim=c(reps,T)); N_pe<-N; N_oe<-N
N[,1]<-N0; N_pe[,1]<-N0
N_oe[,1]<-exp(rnorm(reps,log(N[,1]),sd_o))
grate_pe<-array(0,dim=c(reps,(T-1)))

for (i in 1:reps){
	for (t in 1:(T-1)){
		Time[i,t+1]<-t+1
		N[i,t+1]<-N[i,t]*lambda
		N_pe[i,t+1]<-N_pe[i,t]*(lambda*exp(rnorm(1,0,sd_p)))
		N_oe[i,t+1]<-rlnorm(1,log(N[i,t+1]),sd_o) # same as exp(rnorm(1,log(N[i,t+1]),sd_o))
		grate_pe[i,t]<-N_pe[i,t+1]/N_pe[i,t]
	}
}

plot(1:T,log(N_pe[1,]),type='n',ylim=c(min(log(N_pe)),max(log(N_pe))),xlab='Time', ylab='log(N)',main='Population growth with Process error')
for (i in 1:nrow(N_pe)) {lines(1:T,log(N_pe[i,]))}
curve(f1,1:T,add=TRUE,col='red')

plot(1:T,log(N_oe[1,]),type='n',ylim=c(min(log(N_oe)),max(log(N_oe))),xlab='Time', ylab='log(N)',main='Population growth with Observation error')
for (i in 1:nrow(N_oe)) {lines(1:T,log(N_oe[i,]))}
curve(f1,1:T,add=TRUE,col='red')


plot(N_pe[1,1:T-1],log(grate_pe[1,]),main='Log Annual Growth Rate with Process Error', xlab='N',ylab='log(Nt+1/Nt)')
curve(f2,add=TRUE)

# Estimate mu, lambda and process error: variance in log(Nt+1/Nt) vs. Nt
mu_est<-mean(log(grate_pe))
lam_est_pe<-exp(mu_est)
sd_pe_est<-sd(log(as.vector(grate_pe)))
print(paste('With process error, estimated lambda_g = ',round(lam_est_pe,5),' estimated s_p = ', round(sd_pe_est,5)))

#Estimate lambda and observation error: log(Nt) vs. time
x<-Time; y<-log(N_oe)

#Fit exponential growth model using nonlinear least squares:
fit<-nls(as.vector(y)~log(b1)+log(b2)*as.vector(x),start=list(b1=1,b2=1))
N0_est<-coef(fit)[1]
lam_est_oe<-coef(fit)[2]

# calculate observer variance as mean of squared residuals
v_oe_est<-mean(resid(fit)^2)
sd_oe_est<-sqrt(v_oe_est)

print(paste('With observation error, estimated lambda = ',round(lam_est_oe,5),' estimated s_o = ',round(sd_oe_est,5),' and estimated N0 = ',round(N0_est,5)))

###############################################################################################
###############################################################################################
###############################################################################################