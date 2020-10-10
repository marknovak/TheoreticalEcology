#############################################################################
# Lotka-Volterra Predation
##########################
# Phase portrait
################
source('VectorField.R')

b<-0.5
a<-0.01
e<-0.1
d<-0.2

isoR <- b/a
isoC <- d/(e*a)

model<-function(R,C){
  c(b*R-a*R*C, 
    e*a*R*C-d*C)
  }

xlims<-c(0,3*isoC)
ylims<-c(0,3*isoR)
quartz(width=8,height=8)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),
    pty='s',xaxs='i',yaxs='i',
    mgp=c(1.75,0.4,0), tcl=-0.3)
plotVectorField(model,xlims,ylims,arrow.length=0.08)
	box(lwd=4)
	title(xlab='Resource',ylab='Consumer')
	abline(h=isoR,lwd=2)
	abline(v=isoC,lty=2,lwd=2)

##########################
# Simulate and add to plot
##########################
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

  out<-ode(c(R0=150, C0=40), Time, LV.CR, params)
  lines(out[,2],out[,3])
  
  # Different starting popn sizes
  out<-ode(c(R0=100, C0=20), Time, LV.CR, params)
  lines(out[,2],out[,3])


############################################################
# Evaluate Jacobian at equilibrium and determine eigenvalues
############################################################
D(expression(b*R-a*R*C),'R')
D(expression(b*R-a*R*C),'C')
D(expression(e*a*R*C-d*C),'R')
D(expression(e*a*R*C-d*C),'C')
  
A11<- 0
A12<- -d/e
A21<- e*b
A22<- 0

A<-matrix(c(A11,A12,A21,A22), byrow=T, nrow=2)
eigs=eigen(A)$values
Re(eigs)
Im(eigs)

######################
# MacArthur-Rosenzweig
######################
MR<-function(t,y,params){
	R<-y[1]
	C<-y[2]
	with(as.list(params),{
		dRdt<-b*R*(1-ai*R)-a*R*C/(1+a*h*R)
		dCdt<-e*a*R*C/(1+a*h*R)-d*C
		return(list(c(dRdt,dCdt)))
})}

b<-0.8
ai<-0.1
a<-0.25
e<-0.1
h<-1.5
d<-0.05

model<-function(R,C){
  c(b*R*(1-ai*R)-a*R*C/(1+a*h*R),
    e*a*R*C/(1+a*h*R)-d*C)
  }

isoR<-function(x){ (b + b*x*(a*h-ai-a*ai*h*x))/a}
isoC<-d/(a*e-a*d*h)

xlims<-c(0,10)
ylims<-c(0,8)
quartz(width=8,height=5)
par(mfrow=c(1,2),cex.lab=1.5,cex.axis=0.8,
    mar=c(4,4,1,1),pty='s', xaxs='i',yaxs='i',
    mgp=c(1.75,0.4,0), tcl=-0.3)


plotVectorField(model,xlims,ylims,arrow.length=0.08,grid.points=20)
	box(lwd=4)
	title(xlab='Resource',ylab='Consumer')
	curve(isoR,lwd=2,add=TRUE)
	abline(v=isoC,lty=2,lwd=2)

	params<-c(b=b,a=a,e=e,d=d,ai=ai,h=h)
	Time<-seq(0,1000,by=0.1)
	out<-ode(c(R0=8,C0=6),Time,MR,params)
	lines(out[,2],out[,3],lwd=5,col='red')

	#different starting popn sizes
	out<-ode(c(R0=2,C0=2),Time,MR,params)
	lines(out[,2],out[,3],lwd=3,col='blue')
	
	points(out[nrow(out),2],out[nrow(out),3],pch=21,bg='blue',cex=2)

	################
	# limit cycles
	d<-0.035

	plotVectorField(model,xlims,ylims,arrow.length=0.08,grid.points=20)
	box(lwd=4)
	title(xlab='Resource',ylab='Consumer')
	curve(isoR,lwd=2,add=TRUE)
	abline(v=isoC,lty=2,lwd=2)

	params<-c(b=b,a=a,e=e,d=d,ai=ai,h=h)
	Time<-seq(0,1000,by=0.1)
	out<-ode(c(R0=8,C0=6),Time,MR,params)
	lines(out[,2],out[,3],lwd=5,col='red')

	# different starting popn sizes
	out<-ode(c(R0=4,C0=3),Time,MR,params)
	lines(out[,2],out[,3],lwd=3,col='blue')


#######################################################################
#######################################################################
#######################################################################

