##################################################
##    Function for sampling sigma^2 in the Gibbs sampler     ##
##################################################

gib.sig=function(n_person,len_theta,sig_hat,gamma,X,Y,theta,del = 0.1)
{
 s0=(0.1*sig_hat)/n_person
 s1  =  sig_hat*max(log(n_person), 0.01*len_theta^{2 + del}/n_person)
T1=gamma*s1+(1-gamma)*s0
D=diag(1/T1)
alpha1=10^(-4)
alpha2=10^(-4)
gib.sig=1/rgamma(1,alpha1+(0.5*n_person)+(0.5*len_theta),alpha2+(0.5*t(Y-X%*%theta)%*%(Y-X%*%theta))+(0.5*t(theta)%*%D%*%theta))
return(gib.sig)
}
