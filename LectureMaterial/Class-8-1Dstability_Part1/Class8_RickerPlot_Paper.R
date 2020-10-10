
Dlogis<-function(N0=0.01,K=10,T=40,rd){
    N<-numeric();N[1]<-N0
    for (t in 1:T){     N[t+1] <- N[t]+rd*(1-N[t]/K)*N[t]	 }
    N
}
f<-function(x){ x+rd*(1-x/K)*x }

####################
K=10;rd=1.5

# quartz(width=6,height=6)
op<-par(pty='s',yaxs='i',mar=c(5,5,1,1),cex.lab=2)
	curve(f,0,18,ylim=c(0,16),lwd=5)
	abline(0,1,lwd=3);text(13,13,'f(x)=x',srt=45,pos=3);grid(15,15,col='black')
	legend('top',horiz=T,legend=bquote(list(K==.(K),r[d]==.(rd))),bg='white')
par(op)

# quartz(width=5,height=3)
op<-par(mar=c(5,5,1,1))
	out<-Dlogis(rd=1.5)
	plot(out,type='o',xlab='Time',ylab='N',ylim=c(0,1.1*max(out)))
par(op)

####################
K=10;rd=2.3

# quartz(width=6,height=6)
op<-par(pty='s',yaxs='i',mar=c(5,5,1,1),cex.lab=2)
	curve(f,0,18,ylim=c(0,16),lwd=5)
	abline(0,1,lwd=3);text(13,13,'f(x)=x',srt=45,pos=3);grid(15,15,col='black')
	legend('top',horiz=T,legend=bquote(list(K==.(K),r[d]==.(rd))),bg='white')
par(op)

# quartz(width=5,height=3)
op<-par(mar=c(5,5,1,1))
	out<-Dlogis(rd=2.3)
	plot(out,type='o',xlab='Time',ylab='N',ylim=c(0,1.1*max(out)))
par(op)

####################
K=10;rd=2.5

# quartz(width=6,height=6)
op<-par(pty='s',yaxs='i',mar=c(5,5,1,1),cex.lab=2)
	curve(f,0,18,ylim=c(0,16),lwd=5)
	abline(0,1,lwd=3);text(13,13,'f(x)=x',srt=45,pos=3);grid(15,15,col='black')
	legend('top',horiz=T,legend=bquote(list(K==.(K),r[d]==.(rd))),bg='white')
par(op)

# quartz(width=5,height=3)
op<-par(mar=c(5,5,1,1))
	out<-Dlogis(rd=2.5)
	plot(out,type='o',xlab='Time',ylab='N',ylim=c(0,1.1*max(out)))
par(op)
#####################



pdf('RickerPlots.pdf',heigh=12,width=10)
par(mfcol=c(3,2),pty='s',yaxs='i',mar=c(5,5,1,1),cex.lab=2)
for(i in 1:2){for(rd in c(1.5,2.3,2.5)){
	curve(f,0,18,ylim=c(0,16),lwd=5)
	abline(0,1,lwd=3);text(13,13,'f(x)=x',srt=45,pos=3);grid(15,15,col='black')
	legend('top',horiz=T,legend=bquote(list(K==.(K),r[d]==.(rd))),bg='white')
}}
dev.off()
