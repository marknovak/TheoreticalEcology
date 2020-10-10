
# Define system of ordinary differential equations (ODEs) describing population-level growth rates of each state variable
du1dt<-expression(u1*(1-u1-a12*u2))
du2dt<-expression(rho*u2*(1-u2+a21*u1))

# Take partial derivatives of ODEs with respect to each state variable
Ddu1du1<-D(du1dt,"u1")
Ddu1du2<-D(du1dt,"u2")
Ddu2du1<-D(du2dt,"u1")
Ddu2du2<-D(du2dt,"u2")

Ddu1du1
Ddu1du2
Ddu2du1
Ddu2du2

# Place PDEs into matrix form (Jacobian).  By making it an expression we can evaluate it repeatedly for different parameter values
J<-expression(matrix(c(eval(Ddu1du1),eval(Ddu1du2), eval(Ddu2du1),eval(Ddu2du2)), byrow=T,nrow=2))

# Define coexistence equilibrium in terms of parameters (from Sage)
u1star<-expression((1-a12)/(a12*a21 + 1))
u2star<-expression((1+a21)/(a12*a21 + 1))

# Specify parameters
rho<-1
a12<-0.6;a21<-0.5
# a12<-1.5;a21<-1.7

# Specify which equilibrium and determine its eigenvalues
# Coexistence equilibrium
u1<-eval(u1star);u2<-eval(u2star)
eigen(eval(J))$values

# Sp1 boundary equlibrium
u1<-1;u2<-0
eigen(eval(J))$values

# Trivial equlibrium
u1<-0;u2<-0
eigen(eval(J))$values

