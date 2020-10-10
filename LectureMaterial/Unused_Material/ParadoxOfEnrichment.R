
# Paradox of enrichment

library(deSolve)

predpreyRM<-function(t,y,p){
	H<-y[1]
	P<-y[2]
	with(as.list(p),{
		dHdt<-b*H*(1-alpha*H)-w*P*H/(D+H)
		dPdt<-e*w*P*H/(D+H)-s*P
		return(list(c(dHdt,dPdt)))
	})
}


b<-0.8;e<-0.07;s<-0.2;w<-5;D<-400;alpha<-5e-04;H<-0:(1/alpha)
Hiso<-expression(b/w*(D+(1-alpha*D)*H-alpha*H^2))
HisoStable<-eval(Hiso)

p.RM<-c(b=b,alpha=alpha,e=e,s=s,w=w,D=D)
Time<-100
RM<-ode(c(500,110),1:Time,predpreyRM,p.RM)

par(pty='s',yaxs='i',xaxs='i')
plot(H,eval(Hiso),type='l',ylab="P",xlab="H",ylim=c(0,1.1*max(RM[,3])),lwd=3)

abline(v=s*D/(e*w-s),lty=2,lwd=3)
arrows(RM[-Time,2],RM[-Time,3],RM[-1,2],RM[-1,3],length=0.1)
box(lwd=2)

