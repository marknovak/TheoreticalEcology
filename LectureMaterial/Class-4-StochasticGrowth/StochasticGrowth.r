
################################################################################
# Uniform distribution as in Case pg 37
#################################################################################
library(sfsmisc) # for eaxis function
par(mfrow=c(1,2),mar=c(3,4,1,1),mgp=c(1.75,0.3,0),tcl=-0.2,lwd=1,cex.axis=1,cex.lab=1.4)
Col<-as.vector(col2rgb('black')/255);Col<-rgb(Col[1],Col[2],Col[3], alpha=0.1)
gmean<-function(x){ exp(mean(log(na.omit(x)))) }
amean<-function(x){ mean(na.omit(x)) }

# set stochastic parameters
N0<-2
lambda <- 1.1
l.min<-0.1
l.max<-2.1

TimeSteps<-200; t <- c(1:TimeSteps);reps<-5000;
rlambda<-array(NA,dim=c(reps,length(t)));
N<-array(NA,dim=c(reps,(length(t)+1))); N[,1]<-N0;

for (r in 1:reps){
	for (Time in t){
		rlambda[r,Time]<-runif(1,l.min,l.max)
		N[r,Time+1] <- N[r,Time]*rlambda[r,Time]
	}
}

amu.lambda<-mean(rlambda)
gmu.lambda<-gmean(rlambda)

pN.alam<-pN.glam<-rep(NA,TimeSteps)
pN.alam[1]<-pN.glam[1]<-N0
for(Time in t){
	pN.alam[Time+1]<-pN.alam[Time]*amu.lambda
	pN.glam[Time+1]<-pN.glam[Time]*gmu.lambda
}

g.Nmu<-apply(N,2,gmean)
a.Nmu<-apply(N,2,amean)

ylim<-range(N[,1:50])
plot(N[1,],type='n',log='y',lwd=3,ylab='N_t',xlab='Time', ylim=ylim,xlim=c(0,50),axes=FALSE)
	eaxis(2);axis(1);box(lwd=1)
	for(i in 1:nrow(N)){lines(N[i,],lwd=1,col=Col)}
	lines(pN.alam,col='orange',lwd=2)
	lines(a.Nmu,col='green',lwd=2)
	lines(pN.glam,col='blue',lwd=2)
	lines(g.Nmu,col='red',lwd=2)

legend('bottomleft',legend=c('Arithmetic mean predicted','Arithmetic mean observed','Geometric mean predicted','Geometric mean observed'),lwd=2,col=c('orange','green','blue','red'),bty='n',cex=1,inset=0)

ylim<-range(N)
plot(N[1,],type='n',log='y',lwd=3,ylab='N_t',xlab='Time', ylim=ylim,axes=FALSE)
	eaxis(2);axis(1);box(lwd=1)
	for(i in 1:nrow(N)){lines(N[i,],lwd=1,col=Col)}
	lines(pN.alam,col='orange',lwd=2)
	lines(a.Nmu,col='green',lwd=2)
	lines(pN.glam,col='blue',lwd=2)
	lines(g.Nmu,col='red',lwd=2)

legend('bottomleft',legend=c('Arithmetic mean predicted','Arithmetic mean observed','Geometric mean predicted','Geometric mean observed'),lwd=2,col=c('orange','green','blue','red'),bty='n',cex=1,inset=0)


################################################################################
################################################################################
################################################################################

