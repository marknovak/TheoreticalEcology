#################################################################
######### Theoretical Ecology - Problem set 2 ###################
#################################################################
###############################################
# PART A. STOCHASTIC DENSITY-INDEPENDENT GROWTH
###############################################
lambda <- 1.01  # this line sets the variable "lambda" equal to 2
TimeSteps <- 2500 # sets the variable "TimeSteps" equal to 10
t <- c(1:TimeSteps) # creates a vector from 1 to the number of timesteps
N <- numeric() # creates an empty numeric vector
N[1] <- 1

Col<-as.vector(col2rgb('black')/255) # to convert a named colour into its code
Col<-rgb(Col[1],Col[2],Col[3], alpha=0.2) # to convert the color into a semi-transparent colour; alpha sets the level of transparency

# set stochastic parameters
mu<-0
sigma<-0.2

# Run individual realization
rlambda<-numeric() # to store stochasticlambda values
for (Time in t){
    rlambda[Time]<-lambda*exp(rnorm(1,mu,sigma))
    N[Time+1] <- N[Time]*rlambda[Time]
}
plot(N,log='y',ylim=c(1,10^20))

mean.rlambda<-mean(rlambda);mean.rlambda
# ANSWER: arithmetic mean is > given lambda

#################
# Run and plot all realizations
rlambda<-numeric() # to store stochastic lambda values
reps<-150 # specify number of replicate realizations
plot(N,type='n',log='y',ylim=c(1,10^20)) #  to initialize plot
for (r in 1:reps){
    for (Time in t){
        rlambda[Time]<-lambda*exp(rnorm(1,mu,sigma))
        N[Time+1] <- N[Time]*rlambda[Time]
    }
    lines(N,col=Col)
}

###################
# Run and plot all realizations, saving N and rlambda for all
rlambda<-array(NA,dim=c(reps,length(t))) # to store stochastic lambda values
N<-array(NA,dim=c(reps,(length(t)+1))) # to store popn size values
N[,1]<-1
for (r in 1:reps){
    for (Time in t){
        rlambda[r,Time]<-lambda*exp(rnorm(1,mu,sigma))
        N[r,Time+1] <- N[r,Time]*rlambda[r,Time]
    }
}

gmean<-function(x){ exp(mean(log(x)))}
gmean.rlambda<-gmean(rlambda);gmean.rlambda
mean.rlambda<-mean(rlambda);mean.rlambda

plot(N[1,],log='y',ylim=c(1,10^20),type='n')
for (r in 1:nrow(N)){
    lines(N[r,],col=Col)
}
gN<-numeric()
gN[1]<-1
for (Time in t){
    gN[Time+1] <- gN[Time]*gmean.rlambda
}
lines(gN,col='green',lwd=3)

aN<-numeric()
aN[1]<-1
for (Time in t){
    aN[Time+1] <- aN[Time]*mean.rlambda
}
lines(aN,col='red',lwd=3)

#########################
N.final<-N[,ncol(N)]

hist(log(N.final),col='grey')
abline(v=mean(log(N.final)),col='red',lwd=2)

obs.mu<-mean(log(N.final));obs.mu
obs.var<-var(log(N.final));obs.var

exp.mu<-TimeSteps*log(lambda);exp.mu 
exp.var<-TimeSteps*sigma^2;exp.var    

exp(obs.mu)
exp(exp.mu)
det.N<-N[1]*lambda^TimeSteps;det.N
abline(v=log(det.N),col='green',lwd=2)


#################################################################
# Part B.  STOCHASTIC DENSITY-DEPENDENT GROWTH (Logistic growth)
#################################################################
reps<-4
TimeSteps <- 20
t <- c(1:TimeSteps) # creates a vector from 1 to the number of timesteps
N<-array(NA,dim=c(reps,TimeSteps+1))

# varying N(0)
N[,1]<-c(1,10,50,120)
r.mu<-1
r.sd<-0
K.mu<-100
K.sd<-0

for (rep in 1:reps){
    for (Time in t){
        r<-rnorm(1,r.mu,r.sd)
        K<-rnorm(1,K.mu,K.sd)
        N[rep,Time+1] <- N[rep,Time]+r*(1-(N[rep,Time]/K))*N[rep,Time]
    }
}

plot(N[1,],type='n',ylim=c(0,140))
abline(h=K.mu,col='grey')
for (i in 1:nrow(N)){
    lines(N[i,])
}

pdN<-(N[,-1]-N[,-ncol(N)])/N[,-ncol(N)]
dN<-(N[,-1]-N[,-ncol(N)])
par(mfrow=c(1,2))
plot(dN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(dN)))
for (i in 1:nrow(dN)){
    lines(N[i,-ncol(N)],dN[i,])
}
plot(pdN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(pdN)))
for (i in 1:nrow(pdN)){
    lines(N[i,-ncol(N)],pdN[i,])
}

# varying r
N[,1]<-1
r.mu<-1
r.sd<-0.5
K.mu<-100
K.sd<-0

for (rep in 1:reps){
    for (Time in t){
        r<-rnorm(1,r.mu,r.sd)
        K<-rnorm(1,K.mu,K.sd)
        N[rep,Time+1] <- N[rep,Time]+r*(1-(N[rep,Time]/K))*N[rep,Time]
    }
}

plot(N[1,],type='n',ylim=c(0,140))
abline(h=K.mu,col='grey')
for (i in 1:nrow(N)){
    lines(N[i,])
}

pdN<-(N[,-1]-N[,-ncol(N)])/N[,-ncol(N)]
dN<-(N[,-1]-N[,-ncol(N)])
par(mfrow=c(1,2))
plot(dN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(dN)))
for (i in 1:nrow(dN)){
    lines(N[i,-ncol(N)],dN[i,])
}
plot(pdN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(pdN)))
for (i in 1:nrow(pdN)){
    lines(N[i,-ncol(N)],pdN[i,])
}

# varying K
N[,1]<-1
r.mu<-1
r.sd<-0
K.mu<-100
K.sd<-10

for (rep in 1:reps){
    for (Time in t){
        r<-rnorm(1,r.mu,r.sd)
        K<-rnorm(1,K.mu,K.sd)
        N[rep,Time+1] <- N[rep,Time]+r*(1-(N[rep,Time]/K))*N[rep,Time]
    }
}

plot(N[1,],type='n',ylim=c(0,140))
abline(h=K.mu,col='grey')
for (i in 1:nrow(N)){
    lines(N[i,])
}

pdN<-(N[,-1]-N[,-ncol(N)])/N[,-ncol(N)]
dN<-(N[,-1]-N[,-ncol(N)])
par(mfrow=c(1,2))
plot(dN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(dN)))
for (i in 1:nrow(dN)){
    lines(N[i,-ncol(N)],dN[i,])
}
plot(pdN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(pdN)))
for (i in 1:nrow(pdN)){
    lines(N[i,-ncol(N)],pdN[i,])
}

# varying both r and K
N[,1]<-1
r.mu<-1
r.sd<-0.5
K.mu<-100
K.sd<-10

for (rep in 1:reps){
    for (Time in t){
        r<-rnorm(1,r.mu,r.sd)
        K<-rnorm(1,K.mu,K.sd)
        N[rep,Time+1] <- N[rep,Time]+r*(1-(N[rep,Time]/K))*N[rep,Time]
    }
}

plot(N[1,],type='n',ylim=c(0,140))
abline(h=K.mu,col='grey')
for (i in 1:nrow(N)){
    lines(N[i,])
}

pdN<-(N[,-1]-N[,-ncol(N)])/N[,-ncol(N)]
dN<-(N[,-1]-N[,-ncol(N)])
par(mfrow=c(1,2))
plot(dN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(dN)))
for (i in 1:nrow(dN)){
    lines(N[i,-ncol(N)],dN[i,])
}
plot(pdN[1,]~N[1,-ncol(N)],type='n',ylim=c(0,max(pdN)))
for (i in 1:nrow(pdN)){
    lines(N[i,-ncol(N)],pdN[i,])
}

# IMPLICITLY ASSUMING THAT EFFECT OF ENVIRONMENT ON r IS INDEPENDENT OF IT'S EFFECT ON K

#####################################################################################
#####################################################################################
#####################################################################################


