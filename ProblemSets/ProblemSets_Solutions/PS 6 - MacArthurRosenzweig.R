#################################################################
######### Theoretical Ecology - Problem set 6 ###################
#################################################################
# MacArthur-Rosenzweig
rm(list=ls())
source('VectorField.R')
library(deSolve)

model<-function(R,C){c(b*R*(1-ai*R)-a*R*C/(1+a*h*R), e*a*R*C/(1+a*h*R)-d*C)}

MR<-function(t,y,params){
	R<-y[1]
	C<-y[2]
	with(as.list(params),{
		dRdt<-b*R*(1-ai*R)-a*R*C/(1+a*h*R)
		dCdt<-e*a*R*C/(1+a*h*R)-d*C
		return(list(c(dRdt,dCdt)))
})}

b<-0.8;ai<-0.1;a<-0.25;e<-0.1;h<-1.5;d<-0.05

isoR<-function(x){(b+b*x*(a*h-ai-a*ai*h*x))/a}
isoC<-d/(a*e-a*d*h)

quartz(width=10,height=5)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),pty='s',xaxs='i',yaxs='i',mgp=c(1.75,0.4,0), tcl=-0.3)
par(mfrow=c(1,2))

xlims<-c(0,10);ylims<-c(0,8)
plotVectorField(model,xlims,ylims,arrow.length=0.08,grid.points=20) 
box(lwd=4)
title(xlab='Resource',ylab='Consumer')
curve(isoR,lwd=2,add=T)
abline(v=isoC,lty=2,lwd=2)

params<-c(b=b,a=a,e=e,d=d,ai=ai,h=h)
Time<-seq(0,1000,by=0.1)
out<-ode(c(R0=8,C0=6),Time,MR,params)
lines(out[,2],out[,3],lwd=5,col='red')


# limit cycles
d<-0.035

isoR<-function(x){(b+b*x*(a*h-ai-a*ai*h*x))/a}
isoC<-d/(a*e-a*d*h)
plotVectorField(model,xlims,ylims,arrow.length=0.08,grid.points=20) 
box(lwd=4)
title(xlab='Resource',ylab='Consumer')
curve(isoR,lwd=2,add=T)
abline(v=isoC,lty=2,lwd=2)

params<-c(b=b,a=a,e=e,d=d,ai=ai,h=h)
Time<-seq(0,1000,by=0.1)
out<-ode(c(R0=8,C0=6),Time,MR,params)
lines(out[,2],out[,3],lwd=5,col='red')

##############################################


quartz(width=10,height=3)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),pty='s',xaxs='i',yaxs='i', mgp=c(1.75,0.4,0), tcl=-0.3,mfrow=c(1,3))

# MacArthur-Rosenzweig w/ Type III
b<-0.8;ai<-0.1;a<-0.25;e<-0.1;h<-1.5;m<-3;d<-0.05
isoR<-function(x){-(x^2*ai*b + (a*ai*x^(m + 2) - a*x^(m + 1))*b*h - x*b)/(x^m*a)}
xlims<-c(0,10);ylims<-c(0,8)
curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
box(lwd=4)
legend('topright',legend='TypeII\nFunctional\nResponse',bty='n')

# MacArthur-Rosenzweig w/ Refuge
G<-0.2
isoR<-function(x){(x^3*a*ai*b*h - (G*a*ai*b*h + a*b*h - ai*b)*x^2 + (G*a*b*h -b)*x)/(G*a - x*a)}

curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
box(lwd=4)
legend('topright',legend='Prey\nRefuge',bty='n')

# MacArthur-Rosenzweig w/ Immigration
I<-0.2
isox<-function(x){-(x^3*a*ai*b*h - (a*b*h - ai*b)*x^2 - (I*a*h + b)*x - I)/(x*a)
}
xlims<-c(0,10);ylims<-c(0,8)
curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
box(lwd=4)
legend('topright',legend='Prey\nImmigration',bty='n')

#######################################################################
# The mechanism common to all is that the prey experiences reduced predator-induced mortality at low population sizes.  All, in a sense, provide the prey with a "refuge" at low population sizes.  With a Type III functional response this refuge is phenomenological in that the predator "switches" to some unmodeled prey.  With density-independent immigration, the predator cannot suppress the prey at low population sizes because the prey's growth rate is dependent more upon immigration (which the predator can't control) than it is dependent on the number of prey in the population.
#######################################################################
#######################################################################

