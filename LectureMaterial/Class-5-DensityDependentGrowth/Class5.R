##########################################################################################
##########################################################################################
rm(list=ls()) # clears workspace

# Define function for Discrete logistic
Dlogis<-function(N0,T,r,a){
	N<-numeric();N[1]<-N0
	for(t in 1:T){ N[t+1]<-N[t]+r*(1-a*N[t])*N[t] }
	N
}

#########################
# N(t) vs. t
quartz(width=14,height=3.5)
par(mfrow=c(1,4));par(mar=c(4,4,1.5,1),lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)

	N<-Dlogis(1,20,0.5,0)
	Time<-1:(length(N))
	plot(N~Time,xlab='Time',ylim=c(0,120),type='b')
	legend('bottomright',legend='Density-independent',bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,0.01)
	Time<-1:(length(N))
	plot(N~Time,xlab='Time',ylim=c(0,120),type='b')
	legend('bottomright',legend='Negative\ndensity-dependence',xjust=1,bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,-0.2)
	Time<-1:(length(N))
	plot(N~Time,xlab='Time',ylim=c(0,120),type='b')
	legend('bottomright',legend='Positive\ndensity-dependence',xjust=1,bty='n',cex=1.7)

	Na<-Dlogis(1,20,0.5,0.015)
	Nb<-Dlogis(95,20,0.2,0.015)
	Time<-1:(length(Na))
	plot(Na~Time,xlab='Time',ylim=c(0,120),type='b',ylab="N(t)")
	points(Time,Nb,type='b')

############################
# N(t+1) vs N(t)
quartz(width=12,height=4)
par(mfrow=c(1,3));par(mar=c(4,4,1.5,1),lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3,pty='s')

	N<-Dlogis(1,20,0.5,0)
	plot(N[-1]~N[-length(N)],ylab='N(t+1)',xlab='N(t)',type='b',ylim=c(0,120),xlim=c(0,120))
	abline(0,1,col='grey')
	legend('bottomright',legend='Density-independent',bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,0.01)
	Time<-1:(length(N))
	plot(N[-1]~N[-length(N)],ylab='N(t+1)',xlab='N(t)',type='b',ylim=c(0,120),xlim=c(0,120))
	abline(0,1,col='grey')
	legend('bottomright',legend='Negative\ndensity-dependence',xjust=1,bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,-0.2)
	Time<-1:(length(N))
	plot(N[-1]~N[-length(N)],ylab='N(t+1)',xlab='N(t)',type='b',ylim=c(0,120),xlim=c(0,120))
	abline(0,1,col='grey')
	legend('bottomright',legend='Positive\ndensity-dependence',xjust=1,bty='n',cex=1.7)

############################
# (N(t+1)-N(t)) / N(t) = (N(t+1)/N(t))-1
quartz(width=12,height=4)
par(mfrow=c(1,3));par(mar=c(4,4,1.5,1),lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3,pty='s')

	N<-Dlogis(1,20,0.5,0)
	y<-(N[-1]-N[-length(N)])/N[-length(N)]
	plot(y~N[-length(N)],ylab=expression((N[t+1]-N[t])/N[t]),xlab=expression(N[t]),type='b', ylim=c(0,1),xlim=c(0,110))
	legend('bottomright',legend='Density-independent',bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,0.01)
	Time<-1:(length(N))
	y<-(N[-1]-N[-length(N)])/N[-length(N)]
	plot(y~N[-length(N)],ylab=expression((N[t+1]-N[t])/N[t]),xlab=expression(N[t]),type='b', ylim=c(0,.6),xlim=c(0,110))
	legend('bottomleft',legend='Negative\ndensity-dependence',xjust=1,bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,-0.2)
	Time<-1:(length(N))
	y<-(N[-1]-N[-length(N)])/N[-length(N)]
	plot(y~N[-length(N)],ylab=expression((N[t+1]-N[t])/N[t]),xlab=expression(N[t]),type='b', ylim=c(0,10),xlim=c(0,110))
	legend('bottomright',legend='Positive\ndensity-dependence',xjust=1,bty='n',cex=1.7)

############################
# N(t+1)-N(t)
quartz(width=12,height=4)
par(mfrow=c(1,3));par(mar=c(4,4,1.5,1),lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3,pty='s')

	N<-Dlogis(1,20,0.5,0)
	y<-(N[-1]-N[-length(N)])
	plot(y~N[-length(N)],ylab=expression(N[t+1]-N[t]),xlab=expression(N[t]),type='b', ylim=c(0,100),xlim=c(0,110))
	legend('bottomright',legend='Density-independent',bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,0.01)
	Time<-1:(length(N))
	y<-(N[-1]-N[-length(N)])
	plot(y~N[-length(N)],ylab=expression(N[t+1]-N[t]),xlab=expression(N[t]),type='b', ylim=c(0,13),xlim=c(0,110))
	legend('bottom',legend='Negative\ndensity-dependence',xjust=1,bty='n',cex=1.7)

	N<-Dlogis(1,20,0.5,-0.2)
	Time<-1:(length(N))
	y<-(N[-1]-N[-length(N)])
	plot(y~N[-length(N)],ylab=expression(N[t+1]-N[t]),xlab=expression(N[t]),type='b', ylim=c(0,100),xlim=c(0,110))
	legend('bottomright',legend='Positive\ndensity-dependence',xjust=1,bty='n',cex=1.7)

############################

# N<-Dlogis(1,20,0.5,0.01)
# Time<-1:(length(N))
# y<-(N[-1]-N[-length(N)])/N[-length(N)]
# par(mar=c(4,4,1.5,1),cex.axis=1.5,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3,pty='s',xaxs='i',yaxs='i')
# plot(y~N[-length(N)],ylab='Per capita growth rate',xlab=expression(N[t]),type='l',lty=2,lwd=5, ylim=c(0,.6),xlim=c(0,110),axes=F)
# axis(1,at=c(0,100),labels=c(0,expression(1/alpha)))
# axis(2,at=c(0,0.5),labels=c(0,expression(r[d])),las=1)
# box(lwd=2)

# N<-Dlogis(1,20,0.5,0.01)
# Time<-1:(length(N))
# y<-(N[-1]-N[-length(N)])/N[-length(N)]
# par(mar=c(4,4,1.5,1),cex.axis=1.5,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3,pty='s',xaxs='i',yaxs='i')
# plot(y~N[-length(N)],ylab='Per capita growth rate',xlab=expression(N[t]),type='n',lty=2,lwd=5, ylim=c(0,.6),xlim=c(0,110),axes=F)
# abline(h=0.5,lty=2,lwd=5)
# axis(1,at=c(0,100),labels=c(0,expression(1/alpha)))
# axis(2,at=c(0,0.5),labels=c(0,expression(r[d])),las=1)
# box(lwd=2)

par(mfrow=c(1,2),mar=c(4,4,1.5,1),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	N<-Dlogis(1,20,0.5,0.01)
	Time<-1:(length(N))
	plot(N~Time,xlab='Time',ylim=c(0,120),type='l',lwd=4,axes=F);box(lwd=2)
	axis(2,at=c(0,100),labels=c(0,'K'),las=1)

	y<-(N[-1]-N[-length(N)])
	plot(y~N[-length(N)],ylab=expression(Delta*N[Delta*t]),xlab=expression(N[t]),type='l',xlim=c(0,110), lwd=4,axes=F);box(lwd=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)



#########################################################
library(deSolve)
ClogisK<-function(t,y,p){
	N<-y[1]
	with(as.list(p),{
		dNdt<-r*N*(1-N/K)
		return(list(dNdt))
	})
}

r=0.1
K=100
out<-data.frame(ode(y=c(N=0.01),times=c(1:200),func=ClogisK,parms=c(r=r, K=K)))
N<-out$N
dNdt<- N[-1]-N[-length(N)]
rK4<-(r*K)/4

par(mar=c(4,4,1.5,1),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	axis(1,at=c(0,50,100),labels=c(0,'K/2','K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	segments(50,0,50,rK4,lty=2,lwd=2,col='grey')
	segments(0,rK4,50,rK4,lty=2,lwd=2,col='grey')


# Fixed harvest effort
par(mar=c(4,4,1.5,4),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	segments(0,0,120,2.5,lty=1,lwd=2,col='black')
	mtext('Harvest\nrate',side=4,line=2.5,cex=2)
	text(65,0.4,'New equilibrium')
	arrows(65,0.6,77,1.5,length=0.1)

# Fixed harvest rate
par(mar=c(4,4,1.5,4),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	segments(0,1.5,100,1.5,lty=1,lwd=2,col='black')
	axis(4,at=c(0,rK4),labels=c(0,'maxH'),las=1)
	mtext('Harvest\nrate',side=4,line=2.5,cex=2)
	text(30,0.7,'Unstable')
	text(70,0.7,'Stable')
	arrows(c(30,70),c(0.9,0.9),c(20,80),c(1.4,1.4),length=0.1)

# One series of figures with increasing harvest rates
x<-N[-length(N)]
par(mfrow=c(2,2),mar=c(4,4,1.5,4),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	mtext('Harvest rate',side=4,line=1,cex=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	y<-curve(0.02*x,from=0,to=100,add=T,n=length(x))
	polygon(c(x,rev(y$x)),c(dNdt,rev(y$y)),col='grey')
	arrows(c(20,40,60,90),c(0.7,1,1.3,1.5),c(20,40,60,90),c(1.3,1.7,2,1.2),length=0.05,code=3)

plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	mtext('Harvest rate',side=4,line=1,cex=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	y<-curve(0.05*x,from=0,to=100,add=T,n=length(x))
	polygon(c(x,rev(y$x)),c(dNdt,rev(y$y)),col='grey')
	arrows(c(20,40,60,90),c(1.1,2.1,2.7,2),c(20,40,60,90),c(1.4,2.3,2.5,1.2),length=0.05,code=3)

plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	mtext('Harvest rate',side=4,line=1,cex=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	y<-curve(1.5+0*x,from=0,to=100,add=T,n=length(x))
	polygon(c(x,rev(y$x)),c(dNdt,rev(y$y)),col='grey')
	arrows(c(5,40,60,95),c(1.3,1.7,1.7,1.3),c(5,40,60,95),c(0.8,2.1,2.1,0.8),length=0.05,code=3)

plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	mtext('Harvest rate',side=4,line=1,cex=2)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)
	y<-curve(2+0*x,from=0,to=100,add=T,n=length(x))
	polygon(c(x,rev(y$x)),c(dNdt,rev(y$y)),col='grey')
	arrows(c(5,40,60,95),c(1.3,2.1,2.1,1.3),c(5,40,60,95),c(0.8,2.3,2.3,0.8),length=0.05,code=3)



# Stochasticity
par(mar=c(4,4,1.5,4),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(2,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='N',type='n',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	points(jitter(dNdt,4000)~jitter(N[-length(N)],4000),pch=19)
	axis(1,at=c(0,100),labels=c(0,'K'),las=1)
	axis(2,at=rK4,labels='rK/4',las=1)

##################################
par(mar=c(4,4,1.5,4),cex.axis=1.2,lwd=2,cex.lab=2,mgp=c(1,0.5,0),tcl=-0.3)
	plot(dNdt~N[-length(N)],ylab='dN/dt',xlab='Population size (N)',type='l',xlim=c(0,110),ylim=c(0,2.7), lwd=4,axes=F);box(lwd=2)
	mtext('Harvest',side=4,line=1,cex=2)
	segments(0,0,120,2.5,lty=2,lwd=2,col='black')
	text(90,2.2,'B',cex=2)
	segments(0,1.5,120,1.5,lty=3,lwd=2,col='black')
	text(100,1.2,'A',cex=2)

###################################################################################
###################################################################################
###################################################################################

