rm(list=ls()) # clears workspace
options(stringsAsFactors=F,warn=-1)
###########################################################
###########################################################
# Sum of Polynomials - Lynx-Hare Dynamics
###########################################################
###########################################################
rm(list=ls()) # clears workspace
options(stringsAsFactors=F,warn=-1)
dat<-read.csv("LynxHare.csv", header=TRUE,nrows=59)
x<-dat$Year; y<-dat$Hare; z<-dat$Lynx

quartz(width=8,height=4)
par(mar=c(4,4,2,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),pch=21)
plot(x,y,type='n',xlab='Year',ylab='Popn size')
lines(x,y,type='l',col='red',lwd=4)
lines(x,z,type='l',col='blue',lwd=4)
legend('topleft',c('Hare','Lynx'),lwd=4,col=c('red','blue'),inset=0.01)

####

quartz(width=8,height=4)
par(mar=c(4,4,2,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),pch=21,lwd=2)
plot(x,y,xlab='Year',ylab='Popn size',bg='grey')
legend('topleft',c('Hare'),pch=21,pt.bg='grey',inset=0.01)

#######
nrow(dat) # years of data


#######
quartz(width=8,height=5)
par(mfrow=c(3,1),mar=c(4,4,1,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),lwd=2,pch=21)

p=1
plot(x,y,xlab='Time',ylab='Pop. size',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))

p=3
plot(x,y,xlab='Time',ylab='Pop. size',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))

p=12
plot(x,y,xlab='Time',ylab='Pop. size',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))

quartz(width=8,height=5)
par(mfrow=c(3,1),mar=c(4,4,1,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),lwd=2,pch=21)

p=22
plot(x,y,xlab='Time',ylab='Pop. size',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))

p=26
plot(x,y,xlab='Time',ylab='Pop. size',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))

###########################################################
###########################################################
# Sum of Polynomials - Randomly distributed data
###########################################################
###########################################################
n=20
x<-seq(1,n,1)
y<-runif(n)

########
quartz(width=10,height=5)
par(mfrow=c(2,3),mar=c(4,4,1,1),tcl=-0.2,cex.lab=2,cex.axis=1,mgp=c(2,0.3,0),lwd=2,pch=21)

plot(x,y,xlab='Time',ylab='Random y',bg='grey')

for (p in c(2,6,12,n-1)){
	plot(x,y,xlab='Time',ylab='Random y',bg='grey',main=paste('Polynomal degree: ',p),cex=1.5)
	lines(x,predict(lm(y~poly(x,p)),data.frame(x=x)))
	Sys.sleep(2)
}

#######################################################################
#######################################################################
#######################################################################
#######################################################################
