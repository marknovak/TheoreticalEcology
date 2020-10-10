#######################
# The discrete logistic
########################

Dlogis<-function(N0=0.01,K=10,T=2000,rd=3){
    N<-numeric();N[1]<-N0
    for (t in 1:T){     N[t+1] <- N[t]+rd*(1-N[t]/K)*N[t]	 }
    N
}

# Examine the effects of increasing rd on the population's dynamics
par(mar=c(5,5,1,1))

out<-Dlogis(rd=0.5,T=50)
plot(out,type='o',xlab='Time',ylab='N',ylim=c(0,1.1*max(out)))
out

##################################################################
##################################################################
##################################################################
# Stop here in class
##################################################################
##################################################################
##################################################################

# Examine the effects of minutely different N0 on the population size at time T=2000
out1<-Dlogis(N0=0.01,rd=3)

out2<-Dlogis(N0=0.01+1E-18, rd=3)
# NOTE:  1E-18 is the minimum distinguishable by computer error

par(mfrow=c(1,2))
plot(out1[1950:2000],type='l',col='red')
points(out2[1950:2000],type='l',col='blue')

par(pty='s')
plot(out1,out2,cex=0.5)

##################################################################
##################################################################
##################################################################
# Stop here in class
##################################################################
##################################################################
##################################################################
###################
# Bifurcation plot
##################
# This function returns the values of the min and max of a time-series
peaks <- function(x) {
	if (min(x)==max(x)) {return(min(x))}
	l <- length(x)
	xm1 <- c(x[-1], x[l])
	xp1 <- c(x[1], x[-l])
	z<-unique(x[x > xm1 & x > xp1 | x < xm1 & x < xp1])
	if (length(z)==0){return(min(x))}
	return(z)
}

# Plot a bifurcation diagram as a function of rd
# This might take a while to run
rd.vals<- seq(1,3,0.001)
plot(0,0, xlim=range(rd.vals), ylim=c(0,15),type="n", xlab=bquote(r[d]), ylab="N*")
for (rd in rd.vals) {
	out <- Dlogis(rd=rd,K=10,N0=0.01,T=2000)
	l <- length(out)%/%10 # use only the last 10% of time-steps to avoid transients
	out <- out[(9*l):(10*l)]
	p <- peaks(out)
	points(rep(rd, length(p)), p, pch=".")
}

# OR using apply statements to do all rd values at once.
# This might take a while to run
rd.vals<- seq(1,3,0.001)
out<-sapply(rd.vals, function(rd){ out<-Dlogis(rd=rd,K=10,N0=0.01,T=2000) })
l<-nrow(out)%/%10 # use only the last 10%
out <- out[(9*l):(10*l),]
p<-apply(out,2,peaks)

plot(0,0, xlim=range(rd.vals), ylim=c(0,15),type="n", xlab=bquote(r[d]), ylab="N*")
lapply(1:length(p),function(y){points(rep(rd.vals[y],length(p[[y]])),p[[y]],pch='.')})

##################################################################
##################################################################
##################################################################
# Stop here in class
##################################################################
##################################################################
##################################################################

#########################
# Contrast to dynamics of continuous logistic
#########################
#########################
require(deSolve)

Clogis<-function(t,y,parms){
	N<-y[1]
	with(as.list(parms),{
		dNdt<-r*N*(1-N/K)
		return(list(dNdt))
	})
}

# Examine the effects of increasing r on the population's dynamics
par(mar=c(5,5,1,1))

N0<-c(N=0.01)
t<-1:30
params<-c(K=10,r=1)
out<-data.frame(ode(y=N0,times=t,func=Clogis,parms=params))
plot(out$time,out$N,type='o',ylim=c(0,1.1*max(out$N)))

################
# Bifurcation plot of continuous logistic
plot(0,0, xlim=c(1,3), ylim=c(0,15),type="n", xlab='r', ylab="N*")
for (r in seq(1,3,0.001)) {
	params<-c(K=10,r=r)
	out <- data.frame(ode(y=N0,times=t,func=Clogis,parms=params))[['N']]
	l <- length(out)%/%10 # use only the last 10% of time-steps to avoid transients
	out <- out[(9*l):(10*l)]
	p <- peaks(out)
	points(rep(r, length(p)), p, pch=".")
}

##################################################################
##################################################################
##################################################################
# Stop here in class
##################################################################
##################################################################
##################################################################

#########################
# Contrast to dynamics of continuous logistic with time-delay
#########################
#########################
require(PBSddesolve)

Claglogis <- function(t,y,parms) {
	N<-y[1]
	with(as.list(parms),{
		if (t < tau) {lag <- N}
		else {lag <- pastvalue(t - tau)}
		dNdt<-r*N*(1-lag[1]/K)
		return(list(dNdt))
	})
}

t=seq(0,2000,0.1)
N0=0.01
params=list(r=1, K=10, tau=2)

# Examine the effects of different delay-times (by varying tau > 0)
# solve the dde system
out <- dde(y=N0,times=t,func=Claglogis,parms=params)
plot(out$time,out$y1,type='o',cex=0.5, pch=23,col='grey',xlim=c(0,40))


# Bifurcation plot as function of Tau
plot(0,0, xlim=c(1,3), ylim=c(0,60),type="n", xlab=expression(tau), ylab="N*")
for (tau in seq(1,3,0.01)) {
	params <- list(r=1, K=10, tau=tau)
	out <- data.frame(dde(y=N0,times=t,func=Claglogis,parms=params))[,2]
	l <- length(out)%/%10 # use only the last 10% of time-steps to avoid transients
	out <- out[(9*l):(10*l)]
	p <- peaks(out)
	points(rep(tau, length(p)), p, pch=21,cex=0.2)
}

# Bifurcation plot as function of r
plot(0,0, xlim=c(1,3), ylim=c(0,60),type="n", xlab='r', ylab="N*")
for (R in seq(1,3,0.01)) {
	params <- list(r=R, K=10, tau=1)
	out <- data.frame(dde(y=N0,times=t,func=Claglogis,parms=params))[,2]
	l <- length(out)%/%10 # use only the last 10% of time-steps to avoid transients
	out <- out[(9*l):(10*l)]
	p <- peaks(out)
	points(rep(R, length(p)), p, pch=21,cex=0.2)
}


##################################################################
##################################################################
##################################################################
# Stop here in class
##################################################################
##################################################################
##################################################################
