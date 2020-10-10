###############################
# Taylor expansion example
##########################
pfunc<-function(x){a+b*x+c*x^2+d*x^3+e*x^4}
a = -1; b = 3.24; c = 5.23; d = 2.8; e = 0.38


quartz(height=5,width=7)
par(mar=c(3,3,1,1),lwd=2,cex=1.5,las=1,bty="l",mgp=c(1,0.4,0),cex.lab=2)

x<-seq(-3,2,1)
curve(pfunc(x),xlim=c(-3,2),ylim=c(-15,15),lwd=10,xlab='N',ylab="f(N)",axes=FALSE)
box(lwd=3)
axis(1,at=0,label=0,cex=2)
points(0,a,pch=21,lwd=4,cex=3,bg="red")

legend("bottomright", c("f(N)","constant:      f(0)", expression(linear:~~~~~~~~~~f(0)+f*minute(0)*n),         expression(quadratic:~~~f(0)+f*minute(0)*n+(f*second(0)/2)*n^2),         expression(cubic:~~~~~~~~~~f(0)+f*minute(0)*n+(f*second(0)/2)*n^2+(f*minute*second(0)/6)*n^3)),       lty=c(1,2,1,4,5),lwd=c(5,rep(2,4)),cex=0.75,bty="n",col=c('black','grey','red','grey','grey'))

abline(h=a,lty=2,col='grey',lwd=4)
curve(a+b*x,lty=1,add=TRUE,col='red',lwd=4,lend=2)
curve(a+b*x+c*x^2,lty=4,add=TRUE,col='grey',lwd=4)
curve(a+b*x+c*x^2+d*x^3,lty=5,add=TRUE,col='grey',lwd=4)

