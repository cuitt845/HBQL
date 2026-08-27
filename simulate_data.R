######################################
###  Function for generating simulated data ###
######################################
library(MASS)

simulate_data=function(pi,delta1,delta2,sigma1,sigma2,gamma,n,p1,p2)
{
O1= mvrnorm(n, rep(0, p1), diag(p1))    
A1 = sample(c(-1,1), n, prob=c(pi,1-pi), replace=TRUE)
A2 = sample(c(-1,1), n, prob=c(pi,1-pi), replace=TRUE)
O2 = mvrnorm(n, rep(0, p2), diag(p2))
O2[,1] = delta1*O1[,1] + delta2*A1 + rnorm(n, 0, sigma2)
R1 = rep(0,n)
H_psi = gamma[5]*O1[,1]+ gamma[6]*A1 + gamma[7]*O2[,1] 
R2 = gamma[1] + gamma[2]*O1[,1] + gamma[3]*A1 + gamma[4]*O1[,1]*A1 + H_psi*A2 + rnorm(n, 0, sigma1)     
Y2 = R2

datasets=list(O1=O1,A1=A1,O2=O2,A2=A2,Y2=Y2)
return(datasets)
}