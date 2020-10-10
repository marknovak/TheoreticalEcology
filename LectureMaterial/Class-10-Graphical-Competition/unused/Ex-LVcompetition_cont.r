rm(list=ls())
source('~/Dropbox/Research/R-Codes/aaaFunctions/VectorField.R')

isou1<-function(x){(1-x)/a12}
isou2<-function(x){1-a21*x}

model<-function(u1,u2){c(u1*(1-u1-a12*u2) , rho*u2*(1-u2-a21*u1))}

rho=1

# quartz(width=12,height=4)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,4,1),pty='s',xaxs='i',yaxs='i',mgp=c(1.75,0.4,0), tcl=-0.3)
par(mfrow=c(1,3))

# coexistence
a21=0.5;a12=0.6;
xlims<-ylims<-c(0,2)
plotVectorField(model,xlims,ylims,arrow.length=0.05) 
curve(isou1,add=TRUE,lwd=2)
curve(isou2,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Coexistence',xlab='Sp1',ylab='Sp2')

# species 1 dominant
a21=1.5;a12=0.5;
xlims<-ylims<-c(0,2)
plotVectorField(model,xlims,ylims,arrow.length=0.05) 
curve(isou1,add=TRUE,lwd=2)
curve(isou2,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Dominance by Sp1',xlab='Sp1',ylab='Sp2')

# priority effect
xlims<-ylims<-c(0,1.1)
a21=1.5;a12=1.7;
plotVectorField(model,xlims,ylims,arrow.length=0.05,grid.points=15) 
curve(isou1,add=TRUE,lwd=2)
curve(isou2,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Priority effect',xlab='Sp1',ylab='Sp2')


#######################################################################
#######################################################################
#######################################################################
# Define system of ordinary differential equations (ODEs) describing population-level growth rates of each state variable
du1dt<-expression(u1*(1-u1-a12*u2))
du2dt<-expression(rho*u2*(1-u2+a21*u1))

# Take partial derivatives of ODEs with respect to each state variable
J11<-D(du1dt,"u1")
J12<-D(du1dt,"u2")
J21<-D(du2dt,"u1")
J22<-D(du2dt,"u2")

J11
J12
J21
J22

# Place PDEs into matrix form (Jacobian).  By making it an expression we can evaluate it repeatedly for different parameter values
J<-expression(matrix(c(eval(J11),eval(J12), eval(J21),eval(J22)), byrow=T,nrow=2))

# Define coexistence equilibrium in terms of parameters (from Sage or Mathematica)
u1star<-expression((1-a12)/(a12*a21 + 1))
u2star<-expression((1+a21)/(a12*a21 + 1))

# Specify parameters
rho<-1
a12<-0.6;a21<-0.5 # coexistence parameter values

# Specify which equilibrium and determine its eigenvalues
# Coexistence equilibrium
# Observe that both lambda < 0.  This equilbrium is stable.
u1<-eval(u1star);u2<-eval(u2star)
eigen(eval(J))$values

# Sp1 boundary equlibrium
# Observer that lambda_1 > 0 and lambda_2 < 0.  Thus u2 can invade u1
u1<-1;u2<-0
eigen(eval(J))$values

# Sp2 boundary equlibrium
# Observer that lambda_1 < 0 and lambda_2 > 0.  Thus u1 can invade u2
u1<-0;u2<-1
eigen(eval(J))$values

# Trivial equlibrium
# Observe that both lambda > 0.  This equilbrium is unstable.
u1<-0;u2<-0
eigen(eval(J))$values

#######################################################################
#######################################################################
#######################################################################
