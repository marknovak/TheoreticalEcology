##################################################################
######### Theoretical Ecology - Problem set 3b ###################
#########       Model-fitting & comparison     ###################
##################################################################
rm(list=ls()) # clears workspace
#################################################################
#################################################################
# Import 'greattit.csv' time series from csv file:
# Can also use file.choose
dat<-read.csv(file.choose())

Year<-dat[,1];N_obs<-dat[,2]
T<-length(Year)

log_L_obs<-log(N_obs[-1]/N_obs[-T])

plot(Year,N_obs,ylab='Number observed',xlab='Year',type='o')
plot(N_obs[1:T-1],log_L_obs,ylab='log(N_t+1 / N_t)',xlab='N_t')

# Create functions for 3 alternative models
# i) Exponential model (constant popn growth)
m1<-function(x){b[1]+0*x}
# ii) Ricker model of density-dependence
m2<-function(x){b[1]*(1-x/b[2])}
# iii) Theta-logistic model of density-dependence
m3<-function(x){b[1]*(1-(x/b[2])^b[3])}

# Use nonlinear least-squares to estimate parameters and sigma of each model and plot the resulting best-fit functions (with the data)
x<-N_obs[1:T-1]; y<-as.vector(log_L_obs);n<-length(x)

plot(x,y,xlab='N_t',ylab='log(N_t+1 / N_t)')

# i) Exponential model (constant population growth)
efit<-nls(y~b1+0*x,start=list(b1=0)) 
s_m1<-sqrt(mean(resid(efit)^2))
b<-beta_m1<-coef(efit)
curve(m1,c(min(x):max(x)),add=TRUE,lty=1)

# ii) Ricker model
rfit<-nls(y~b1*(1-x/b2), start=list(b1=beta_m1[1],b2=1))
s_m2<-sqrt(mean(resid(rfit)^2))
b<-beta_m2<-coef(rfit)
curve(m2,c(min(x):max(x)),add=TRUE,lty=2)

# iii) Theta-logistic
tfit<-nls(y~b1*(1-(x/b2)^b3),start=list(b1=beta_m2[1],b2=beta_m2[2],b3=1))
s_m3<-sqrt(mean(resid(tfit)^2))
b<-beta_m3<-coef(tfit)
curve(m3,c(min(x):max(x)),add=TRUE,lty=3)

legend('topright',legend=c('exp','ricker','theta'),lty=1:3)


MLE_results<-array(0,dim=c(3,6))
for (f in 1:3){ 
	if (f==1){  # Exponential model (constant population growth), # param = 2 (beta_1 and sigma)
		b=beta_m1;sigma=s_m1;np=2
		y_exp<-m1(x)
	}
	if (f==2){  # Ricker model # param = 3 (beta_1 beta_2 and sigma)
		b=beta_m2;sigma=s_m2;np=3
		y_exp<-m2(x)
	}
	if (f==3){  # Theta-logistic model # param = 4 (beta_1 beta_2 and beta_3 and sigma)
		b=beta_m3;sigma=s_m3;np=4
		y_exp<-m3(x)
	}
	nLL<-(n/2)*log(2*pi)+(n/2)*log(sigma^2)+(1/(2*sigma^2))*sum((y-y_exp)^2)
	AIC<-2*nLL + 2*np
	AICc<-2*nLL + 2*np*(n/(n-np-1))
	MLE_results[f,1:4]<-c(f,nLL,AIC,AICc)
}

# Calculate delta AIC and AIC-weights
delta_AICc<-MLE_results[,4]-min(MLE_results[,4])

for (f in 1:3){
	AIC_wt<-exp(-0.5*delta_AICc[f])/sum(exp(-0.5*delta_AICc))
	MLE_results[f,5]<-delta_AICc[f]
	MLE_results[f,6]<-AIC_wt
}
colnames(MLE_results)<-c('p','nLL','AIC','AICc','delta_AICc','AIC_wt')
print(MLE_results)

# Conclusion:  Ricker model performs best.

########################################################################
########################################################################
########################################################################










