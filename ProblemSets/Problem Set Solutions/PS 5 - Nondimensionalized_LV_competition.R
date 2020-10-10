#################################################################
######### Theoretical Ecology - Problem set 5 ###################
#################################################################
# sage: # define parameters and variables
# sage: u, a12, u2, rho, a21 = var('u1 a12 u2 rho a21')
# sage: du1dtau=u1*(1-u1-a12*u2)
# sage: du2dtau=r
# sage: # Solve for zero-isoclines
# sage: isou1=solve(du1dtau==0,u2)
# sage: print(isou1);
# sage: isou2=solve(du2dtau==0,u1)
# sage: print(isou2);
# sage: # rearrange in order to plot 2nd isocline in R
# sage: solve(isou2,u2)
# [
# u2 == -(u1 - 1)/a12
# ]
# [
# u1 == -(u2 - 1)/a21
# ]
# [u2 == -a21*u1 + 1]
# sage: # determine equilibria
# sage: eq1= -(u1 - 1)/a12 

##################################################################################
##################################################################################
rm(list=ls())
	
source('VectorField.R')
library(deSolve)

isou1<-function(x){-a21*x + 1}
isou2<-function(x){-(x - 1)/a12}

rho=1 # species have equal density-independent growth rates

par(mfrow=c(2,2),pty='s',xaxs='i',yaxs='i')
# species 1 dominant
a21=1.5;a12=0.5;
plotVectorField(function(u1,u2){c(u1*(1-u1-a12*u2),rho*u2*(1-u2-a21*u1))},c(0,2),c(0,2)) 
curve(isou2,add=T)
curve(isou1,add=T,lty=2)
legend('topright',legend=c('u1','u2'),lty=c(1,2),bg='white')
title(main='Dominance by u1',xlab='u1',ylab='u2')

# species 2 dominant
a21=0.5;a12=1.5;
plotVectorField(function(u1,u2){c(u1*(1-u1-a12*u2),rho*u2*(1-u2-a21*u1))},c(0,2),c(0,2)) 
curve(isou2,add=T)
curve(isou1,add=T,lty=2)
legend('topright',legend=c('u1','u2'),lty=c(1,2),bg='white')
title(main='Dominance by u2',xlab='u1',ylab='u2')

# coexistence
a21=0.5;a12=0.7;
plotVectorField(function(u1,u2){c(u1*(1-u1-a12*u2),rho*u2*(1-u2-a21*u1))},c(0,2),c(0,2)) 
curve(isou2,add=T)
curve(isou1,add=T,lty=2)
legend('topright',legend=c('u1','u2'),lty=c(1,2),bg='white')
title(main='Coexistence',xlab='u1',ylab='u2')

# priority effect
a21=1.5;a12=1.7;
plotVectorField(function(u1,u2){c(u1*(1-u1-a12*u2),rho*u2*(1-u2-a21*u1))},c(0,2),c(0,2)) 
curve(isou2,add=T)
curve(isou1,add=T,lty=2)
legend('topright',legend=c('u1','u2'),lty=c(1,2),bg='white')
title(main='Priority effect',xlab='u1',ylab='u2')

##################################################################################
model<-function(t, x, params){
  u1<-x[1];
  u2<-x[2];
with(as.list(parameters),{
    du1dt<-u1*(1-u1-a12*u2)
    du2dt<-rho*u2*(1-u2-a21*u1)
    res<-c(du1dt,du2dt)
    list(res)
 })}

T<-20
t<-seq(0, T, by=1)
xstart <- c(u1=0.1,u2=0.1)

par(mfrow=c(3,2))

# species 1 dominant
parameters<-c(a21=1.5,a12=0.5,rho=1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1~time,data=out,type='l',ylim=c(0,1.5))
lines(out$time,out$u2,lty=2)
legend('topleft',legend=c('u1','u2'),lty=c(1,2),bg='white')
title('Dominance by u1')

# species 2 dominant
parameters<-c(a21=0.5,a12=1.5,rho=1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1~time,data=out,type='l',ylim=c(0,1.5))
lines(out$time,out$u2,lty=2)
legend('topleft',legend=c('u1','u2'),lty=c(1,2),bg='white')
title('Dominance by u2')

# coexistence
parameters<-c(a21=0.5,a12=0.7,rho=1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1~time,data=out,type='l',ylim=c(0,1.5))
lines(out$time,out$u2,lty=2)
legend('topleft',legend=c('u1','u2'),lty=c(1,2),bg='white')
title('Coexistence')

# priority effect - starting abundance 1
xstart <- c(u1=0.1,u2=0.1)
parameters<-c(a21=1.5,a12=1.7,rho=1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1~time,data=out,type='l',ylim=c(0,1.5))
lines(out$time,out$u2,lty=2)
legend('topleft',legend=c('u1','u2'),lty=c(1,2),bg='white')
title('Priority effect')

# priority effect - starting abundance 2
xstart <- c(u1=0.5,u2=0.1)
parameters<-c(a21=1.5,a12=1.7,rho=2)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1~time,data=out,type='l',ylim=c(0,1.5))
lines(out$time,out$u2,lty=2)
legend('topleft',legend=c('u1','u2'),lty=c(1,2),bg='white')
title('Priority effect')

##################################################################################
##################################################################################
##################################################################################

