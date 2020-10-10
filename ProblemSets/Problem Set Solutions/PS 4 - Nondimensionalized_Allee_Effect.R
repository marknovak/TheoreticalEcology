#################################################################
######### Theoretical Ecology - Problem set 4 ###################
#################################################################

rm(list=ls())
library(deSolve)
#####################################
model<-function(t, x, params){
  N<-x[1];
with(as.list(parameters),{
    dN <- r*N*(1-N/K)*(N/K-A/K)
    res<-c(dN)
    list(res)
 })}

T<-100
t<-seq(0, T, by=1)
parameters  <- c(r=0.5,K=2,A=1)

quartz(width=8,height=8)
par(mfrow=c(2,2))

xstart <- c(N=0.4)
out <- as.data.frame(lsoda(xstart, t, model, parameters))
plot(N~time,data=out,type='l',ylim=c(0,2.2))

xstart <- c(N=1.1)
out <- as.data.frame(lsoda(xstart, t, model, parameters))
plot(N~time,data=out,type='l',ylim=c(0,2.2))

#####################################
model2<-function(t, X, params){
  x<-X[1];
with(as.list(parameters),{
    dx <- x*(1-x)*(x-b)
    res<-c(dx)
    list(res)
 })}

Tau<-100*0.5 # r=0.5
t<-seq(0, Tau, by=1)
parameters  <- c(b=0.5) # A/K = 0.5

xstart <- c(x=0.2)
out <- as.data.frame(lsoda(xstart, t, model2, parameters))
plot(x~time,data=out,type='l',ylim=c(0,1.1))

xstart <- c(x=0.55)
out <- as.data.frame(lsoda(xstart, t, model2, parameters))
plot(x~time,data=out,type='l',ylim=c(0,1.1))

#################################################################
#################################################################
#################################################################
