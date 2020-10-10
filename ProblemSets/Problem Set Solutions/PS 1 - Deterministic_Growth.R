#################################################################
######### Theoretical Ecology - Problem set 1 ###################
#################################################################
rm(list=ls()) # Clears workspace

lambda <- 2  # this line sets the variable "lambda" equal to 2
TimeSteps <- 10 # sets the variable "TimeSteps" equal to 10
t <- c(1:TimeSteps) # creates a vector from 1 to the number of timesteps
N <- numeric() # creates an empty numeric vector
N[1] <- 1

for (Time in t){ N[Time+1] <- N[Time]*lambda }

T<-c(t,TimeSteps+1)
plot(N)
plot(T,N)
plot(T,N,log='y')

lambda<-1.8
for (Time in t){ N[Time+1] <- N[Time]*lambda }
lines(T,N)

dN<-N[-1]/N[-length(N)]
plot(dN~N[-length(N)],log='x')
abline(h=lambda, lty=2, col='red')


lambda <- 0.8
TimeSteps <- 30
t <- c(1:TimeSteps)
N <- numeric()
N[1] <- 1000

for (Time in t){ N[Time+1] <- N[Time]*lambda }
qe<-length(N[which(N>N[1]/100)]) # timesteps to quasi-extinction


T<-c(t,TimeSteps+1)

par(mfrow=c(1,2))

plot(T,N,pch=19,cex=0.3,type='o')
abline(h=N[1]*0.01)
abline(v=qe+1) # time-point of quasi-extinction

plot(T,N,log='y',pch=19,cex=0.3,type='o')
abline(h=N[1]*0.01)
abline(v=qe+1)  # time-point of quasi-extinction


#################################################################
#################################################################
#################################################################
#################################################################