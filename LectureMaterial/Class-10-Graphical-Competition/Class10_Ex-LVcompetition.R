################
# LV competition
################
library(deSolve)

model<-function(t, x, params){
  N1<-x[1];
  N2<-x[2];
with(as.list(params),{
    dN1dt<-r1*N1*(1-N1/K1-a12*N2/K1)
    dN2dt<-r2*N2*(1-N2/K2-a21*N1/K2)
    out<-c(dN1dt,dN2dt)
    list(out)
 })}

##########
T<-30
t<-seq(0, T, by=1)

xstart <- c(N1=0.01,N2=0.0)
parameters<-c(a21=0.5,a12=0.7,r1=1,r2=1,K1=1,K2=1)
###########
out <- as.data.frame(ode(xstart, t, model, parameters))

#########
quartz(height=4,width=10)
par(mfcol=c(1,2))
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),mgp=c(1.75,0.4,0), tcl=-0.3)

plot(N1~time,data=out,type='l',ylim=c(0,1),lwd=2,ylab='Popn size')
lines(out$time,out$N2,lty=2,lwd=2)
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2))
box(lwd=3)

op<-par(pty='s')
lims<-c(0,1)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='n')
abline(0,1,lty=2,col='grey',lwd=2)
abline(v=tail(out,1)[2],col='grey',lty=3)
abline(h=tail(out,1)[3],col='grey',lty=3)
box(lwd=3)

points(out$N1,out$N2,type='o',pch=21,bg='grey')
par(op)


####################################################
####################################################
####################################################
 # Repeat above to demonstrate with N2(0) = 0
####################################################
####################################################
####################################################
 #  				PAUSE
####################################################
####################################################
####################################################
# 		Continue for 'Graphical analysis'
####################################################
####################################################
####################################################

isoSp1<-function(x){(K1-x)/a12}
isoSp2<-function(x){K2-a21*x}
####################################################

xstart <- c(N1=0.1,N2=0.1)

quartz(width=10,height=5)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),mgp=c(1.75,0.4,0), tcl=-0.3)
par(mfcol=c(2,5),yaxs='i',xaxs='i')

####################
# coexistence
####################
parameters<-c(a21=0.7,a12=0.5,r1=1,r2=1,K1=1,K2=1)
K1=1;a12=0.5;
K2=1;a21=0.7;
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(N1~time,data=out,type='l',ylim=c(0,1.5),lwd=2,ylab='Popn size',col='blue',axes=FALSE)
axis(1,labels=NA);axis(2,labels=NA)
lines(out$time,out$N2,lty=2,lwd=2,col='red')
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2),col=c('blue','red'),bg='white')
box(lwd=3)
title('Coexistence')

op<-par(pty='s')
lims<-c(0,2)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='o',axes=FALSE)
axis(1,at=c(K1,K2/a21),labels=c('K1','K2/a21'))
axis(2,at=c(K2,K1/a12),labels=c('K2','K1/a12'))
box(lwd=3)
# abline(0,1,lty=2,col='grey')
curve(isoSp1,seq(0,2,0.1),add=TRUE,col='blue',lwd=2)
curve(isoSp2,seq(0,2,0.1),add=TRUE,col='red',lwd=2,lty=2)
par(op)

####################
# species 1 dominant
####################
parameters<-c(a21=1.5,a12=0.5,r1=1,r2=1,K1=1,K2=1)
K1=1;a12=0.5;
K2=1;a21=1.5;
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(N1~time,data=out,type='l',ylim=c(0,1.5),lwd=2,ylab='Popn size',col='blue',axes=FALSE)
axis(1,labels=NA);axis(2,labels=NA)
lines(out$time,out$N2,lty=2,lwd=2,col='red')
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2),col=c('blue','red'),bg='white')
box(lwd=3)
title('Sp1 Dominant')

op<-par(pty='s')
lims<-c(0,2)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='o',axes=FALSE)
axis(1,at=c(K1,K2/a21),labels=c('K1','K2/a21'))
axis(2,at=c(K2,K1/a12),labels=c('K2','K1/a12'))
box(lwd=3)
# abline(0,1,lty=2,col='grey')
curve(isoSp1,seq(0,2,0.1),add=TRUE,col='blue',lwd=2)
curve(isoSp2,seq(0,2,0.1),add=TRUE,col='red',lwd=2,lty=2)
par(op)

####################
# species 2 dominant
####################
parameters<-c(a21=0.5,a12=1.5,r1=1,r2=1,K1=1,K2=1)
K1=1;a12=1.5;
K2=1;a21=0.5;
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(N1~time,data=out,type='l',ylim=c(0,1.5),lwd=2,ylab='Popn size',col='blue',axes=FALSE)
axis(1,labels=NA);axis(2,labels=NA)
lines(out$time,out$N2,lty=2,lwd=2,col='red')
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2),col=c('blue','red'),bg='white')
box(lwd=3)
title('Sp2 Dominant')

op<-par(pty='s')
lims<-c(0,2)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='o',axes=FALSE)
axis(1,at=c(K1,K2/a21),labels=c('K1','K2/a21'))
axis(2,at=c(K2,K1/a12),labels=c('K2','K1/a12'))
box(lwd=3)
# abline(0,1,lty=2,col='grey')
curve(isoSp1,seq(0,2,0.1),add=TRUE,col='blue',lwd=2)
curve(isoSp2,seq(0,2,0.1),add=TRUE,col='red',lwd=2,lty=2)
par(op)

####################
# priority effect - N1<N2 starting abundance
####################
xstart <- c(N1=0.05,N2=0.15)
parameters<-c(a21=2,a12=2,r1=1,r2=1,K1=1.4,K2=1.4)
K1=1.4;a12=2;
K2=1.4;a21=2;
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(N1~time,data=out,type='l',ylim=c(0,1.5),lwd=2,ylab='Popn size',col='blue',axes=FALSE)
axis(1,labels=NA);axis(2,labels=NA)
lines(out$time,out$N2,lty=2,lwd=2,col='red')
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2),col=c('blue','red'),bg='white')
box(lwd=3)
title('Alternative Stable States')

op<-par(pty='s')
lims<-c(0,2)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='o',axes=FALSE)
axis(1,at=c(K1,K2/a21),labels=c('K1','K2/a21'))
axis(2,at=c(K2,K1/a12),labels=c('K2','K1/a12'))
box(lwd=3)
# abline(0,1,lty=2,col='grey')
curve(isoSp1,seq(0,2,0.1),add=TRUE,col='blue',lwd=2)
curve(isoSp2,seq(0,2,0.1),add=TRUE,col='red',lwd=2,lty=2)
par(op)

####################
# priority effect - N1>N2 starting abundance
####################
xstart <- c(N1=0.15,N2=0.05)
parameters<-c(a21=2,a12=2,r1=1,r2=1,K1=1.4,K2=1.4)
K1=1.4;a12=2;
K2=1.4;a21=2;
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(N1~time,data=out,type='l',ylim=c(0,1.5),lwd=2,ylab='Popn size',col='blue',axes=FALSE)
axis(1,labels=NA);axis(2,labels=NA)
lines(out$time,out$N2,lty=2,lwd=2,col='red')
legend('topleft',legend=c('Sp1','Sp2'),lty=c(1,2),col=c('blue','red'),bg='white')
box(lwd=3)
title('Alternative Stable States')

op<-par(pty='s')
lims<-c(0,2)
plot(N2~N1,data=out,ylim=lims,xlim=lims,ylab='Sp2',xlab='Sp1',type='o',axes=FALSE)
axis(1,at=c(K1,K2/a21),labels=c('K1','K2/a21'))
axis(2,at=c(K2,K1/a12),labels=c('K2','K1/a12'))
box(lwd=3)
# abline(0,1,lty=2,col='grey')
curve(isoSp1,seq(0,2,0.1),add=TRUE,col='blue',lwd=2)
curve(isoSp2,seq(0,2,0.1),add=TRUE,col='red',lwd=2,lty=2)
par(op)

####################################################
####################################################
####################################################

# source vector field function
source('VectorField.R')

isoN1<-function(x){-x*a21+K2}
isoN2<-function(x){(K1-x)/a12}

model<-function(N1,N2){c(r1*N1*(1-N1/K1-a12*N2/K1) , r2*N2*(1-N2/K2-a21*N1/K2))}

r1<-r2<-1
K1<-K2<-1


quartz(width=10,height=3.5)
par(cex.lab=1.5,cex.axis=0.8,mar=c(4,4,1,1),pty='s',xaxs='i',yaxs='i',mgp=c(1.75,0.4,0), tcl=-0.3)
par(mfrow=c(1,4))

# coexistence
a21=0.5;a12=0.7;
xlims<-ylims<-c(0,2)
plotVectorField(model,xlims,ylims)
curve(isoN2,add=TRUE,lwd=2)
curve(isoN1,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Coexistence',xlab='Sp1',ylab='Sp2')

# species 1 dominant
a21=1.5;a12=0.5;
xlims<-ylims<-c(0,2)
plotVectorField(model,xlims,ylims)
curve(isoN2,add=TRUE,lwd=2)
curve(isoN1,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Dominance by Sp1',xlab='Sp1',ylab='Sp2')

# species 2 dominant
a21=0.5;a12=1.5;
xlims<-ylims<-c(0,2)
plotVectorField(model,xlims,ylims)
curve(isoN2,add=TRUE,lwd=2)
curve(isoN1,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Dominance by Sp2',xlab='Sp1',ylab='Sp2')

# priority effect
a21=1.5;a12=1.7;
xlims<-ylims<-c(0,1.1)
plotVectorField(model,xlims,ylims,grid.points=15)
curve(isoN2,add=TRUE,lwd=2)
curve(isoN1,add=TRUE,lty=2,lwd=2)
box(lwd=4)
legend('topright',legend=c('Sp1','Sp2'),lty=c(1,2),bg='white')
title(main='Priority effect',xlab='Sp1',ylab='Sp2')

