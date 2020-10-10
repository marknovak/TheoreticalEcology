source('~/Dropbox/R-Codes/aaaFunctions/VectorField.R')

model<-function(R,C){c(b*R-a*R*C, e*a*R*C-d*C)}

b<-0.5;a<-0.01;e<-0.1;d<-0.2
isoR<-b/a
isoC<-d/(e*a)

quartz(width=7,height=7)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),pty='s',xaxs='i',yaxs='i',mgp=c(1.75,0.4,0), tcl=-0.3)

xlims<-c(0,3*isoC);ylims<-c(0,3*isoR)
plotVectorField(model,xlims,ylims,arrow.length=0.08) 
box(lwd=4)
title(xlab='Resource',ylab='Consumer')
abline(h=isoR,lwd=2)
abline(v=isoC,lty=2,lwd=2)

##################################################################
library(deSolve)

LV.CR<-function(t,y,params){
	R<-y[1]
	C<-y[2]
	with(as.list(params),{
		dRdt<-b*R-a*R*C
		dCdt<-e*a*R*C-d*C
		return(list(c(dRdt,dCdt)))
	})
}

params<-c(b=b,a=a,e=e,d=d)
Time<-seq(0,100,by=0.1)

out<-ode(c(R0=150,C0=40),Time,LV.CR,params)
lines(out[,2],out[,3])

# Different starting popn sizes
out<-ode(c(R0=100,C0=20),Time,LV.CR,params)
lines(out[,2],out[,3])

#####################################################################
#####################################################################
#####################################################################

# Evaluation Jacobian at equilibrium

J11<- 0
J12<- -d/e
J21<- e*b
J22<- 0

J<-matrix(c(J11,J12,J21,J22),byrow=T,nrow=2)
eigen(J)$values

#####################################################################
#####################################################################
#####################################################################
# Lynx-Hare data
options(stringsAsFactors=F,warn=-1)
dat<-read.csv("/Users/marknovak/Dropbox/BIOE-148/LectureMaterial/Class-1/LynxHare.csv", header=TRUE,nrows=59)
x<-dat$Year; y<-dat$Hare; z<-dat$Lynx;n<-length(z)

quartz(width=8,height=4)
par(mar=c(4,4,2,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),pch=21)
plot(x,y,type='n',xlab='Year',ylab='Popn size')
lines(x,y,type='l',col='red',lwd=4)
lines(x,z,type='l',col='blue',lwd=4)
legend('topleft',c('Hare','Lynx'),lwd=4,col=c('red','blue'),inset=0.01)

quartz(width=4,height=4)
par(mar=c(4,4,2,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),pty='s')
plot(y,z,pch=19,type='n',cex=0.4,xlab='Hare',ylab='Lynx');box(lwd=2)
arrows(y[-n],z[-n],y[-1],z[-1],length=0.1)
#####################################################################
#####################################################################
#####################################################################
# MacArthur-Rosenzweig

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

#different starting popn sizes
out<-ode(c(R0=2,C0=2),Time,MR,params)
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

# different starting popn sizes
out<-ode(c(R0=4,C0=3),Time,MR,params)
lines(out[,2],out[,3],lwd=5,col='red')

d<-0.02;a<-0.14
isoR<-function(x){(b+b*x*(a*h-ai-a*ai*h*x))/a}
isoC<-d/(a*e-a*d*h)
curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
abline(v=isoC,lty=2,lwd=2);box(lwd=4)

##############################################
quartz(width=10,height=5)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),pty='s',xaxs='i',yaxs='i',mgp=c(1.75,0.4,0), tcl=-0.3,mfrow=c(1,2))
# MacArthur-Rosenzweig w/ Type III
b<-0.8;ai<-0.1;a<-0.25;e<-0.1;h<-1.5;m<-3;d<-0.05
isoR<-function(x){-(x^2*ai*b + (a*ai*x^(m + 2) - a*x^(m + 1))*b*h - x*b)/(x^m*a)}

xlims<-c(0,10);ylims<-c(0,8)
curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
box(lwd=4)

# MacArthur-Rosenzweig w/ Refuge
G<-0.2
isoR<-function(x){(x^3*a*ai*b*h - (G*a*ai*b*h + a*b*h - ai*b)*x^2 + (G*a*b*h -
b)*x)/(G*a - x*a)}

xlims<-c(0,10);ylims<-c(0,8)
curve(isoR,lwd=2,xlim=xlims,ylim=ylims,xlab='Resource',ylab='Consumer')
box(lwd=4)

#######################################################################
#######################################################################
#######################################################################
