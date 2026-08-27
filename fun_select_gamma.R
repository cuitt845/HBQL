###############################################################
#### Function for selecting active covariates at stages 2 and 1 using HBQL ###  
###############################################################
# nsplit: number of coefficient blocks used in block Gibbs sampling; 
# it may be adjusted according to the dimensions p1 and p2.
fun_select_gamma=function(datasets0,KK2,KK1,niter,nburn)
{
##########################
##### read data ############
O1=datasets0$O1
A1=datasets0$A1
O2=datasets0$O2
A2=datasets0$A2
Y2=datasets0$Y2
Y1=0
S1=cbind(rep(1,n),O1) 
H1=cbind(rep(1,n),O1) 
S2=cbind(rep(1,n), O1, A1, O1*A1)  
H2 = cbind(O1, A1, O2) 
X1 = cbind(S1, H1*A1)     
X2 = cbind(S2, H2*A2)    
n=length(O1[,1])
##########################################
len_gamma2=len_theta2=dim(X2)[2]
len_gamma1=len_theta1=dim(X1)[2]
len_beta2=dim(S2)[2]
len_psi2=dim(H2)[2]
len_beta1=dim(S1)[2]
len_psi1=dim(H1)[2]
theta2.sample=matrix(rep(0,len_theta2*niter),nrow=len_theta2,ncol=niter)
gamma2.sample=matrix(rep(5,len_theta2*niter),nrow=len_theta2,ncol=niter)
theta1.sample=matrix(rep(0,len_theta1*niter),nrow=len_theta1,ncol=niter)
gamma1.sample=matrix(rep(5,len_theta1*niter),nrow=len_theta1,ncol=niter)
A2.opt.sample=matrix(rep(0,n*niter),nrow=n,ncol=niter)
Y.opt.sample=matrix(rep(0,n*niter),nrow=n,ncol=niter)
A1.opt.sample=matrix(rep(0,n*niter),nrow=n,ncol=niter)
sig2.sample=rep(0,niter)
sig1.sample=rep(0,niter)
##################### Stage 2 ###############################
choicep2=function(x)
{   
  return(x - KK2 + qnorm(0.9)*sqrt(x*(1- x/len_theta2)))
}
cp2 = uniroot(choicep2, c(1,KK2))$root
q_n2=cp2/len_theta2
## initial values at stage 2 ##
theta2.initial=theta.initial(X2,Y2,len_theta2)
sig_hat2=var(Y2)
gamma2=rep(0,len_gamma2)
## HBIG at stage 2 ##
gib.sig2=gib.sig(n,len_theta2,sig_hat2,gamma2,X2,Y2,theta2.initial,0.1)
gib.theta2=gib.theta(X2, Y2, theta2.initial, gib.sig2, gamma2,del = 0.1, nsplit =50)
gib.gamma2=gib.gamma(n,len_theta2,gib.sig2,q_n2,gib.theta2)
gib.psi2=gib.theta2[(len_beta2+1):(len_theta2)]
A2.opt.sample[,1]=sign(H2%*%gib.psi2)
for(i in 1:n)
{
	if(A2.opt.sample[i,1]==A2[i])
       {
       	  Y.opt.sample[i,1]=Y2[i]
       }else{
       	  Y.opt.sample[i,1]=rnorm(1,(c(S2[i,],H2[i,]*A2.opt.sample[i,1]))%*%gib.theta2,sqrt(gib.sig2))
      }
}
sig2.sample[1]=gib.sig2
theta2.sample[,1]=gib.theta2
gamma2.sample[,1]=gib.gamma2
for(brn in 2:nburn)
{
    sig2.sample[brn]=gib.sig(n,len_theta2,sig2.sample[brn-1],gamma2.sample[,brn-1],X2,Y2,theta2.sample[,brn-1],del = 0.1)
    theta2.sample[,brn]=gib.theta(X2, Y2, theta2.sample[,brn-1], sig2.sample[brn], gamma2.sample[,brn-1],del = 0.1,nsplit =50)
    gamma2.sample[,brn]=gib.gamma(n,len_theta2,sig2.sample[brn],q_n2,theta2.sample[,brn])
    gib.psi2=theta2.sample[(len_beta2+1):(len_theta2),brn]
    A2.opt.sample[,brn]=sign(H2%*%gib.psi2)
    X2.opt=cbind(S2,H2*A2.opt.sample[,brn])
    for(i in 1:n)
    {
	     if(A2.opt.sample[i,brn]==A2[i])
          {
       	        Y.opt.sample[i,brn]=Y2[i]
          }else{
       	        Y.opt.sample[i,brn]=rnorm(1,X2.opt[i,]%*%theta2.sample[,brn],sqrt(sig2.sample[brn]))
         }
    }
}
######################### Stage 1  ################################################
choicep1=function(xx)
{   
  return(xx - KK1 + qnorm(0.9)*sqrt(xx*(1- xx/len_theta1)))
}
cp1 = uniroot(choicep1, c(1,KK1))$root
q_n1=cp1/len_theta1
## initial values at stage 1 ##
theta1.initial=theta.initial(X1,(Y1+Y.opt.sample[,1]),len_theta1)
sig_hat1=var(Y1+Y.opt.sample[,1])
gamma1=rep(0,len_gamma1)
## HBIG at stage 1 ##
gib.sig1=gib.sig(n,len_theta1,sig_hat1,gamma1,X1,(Y1+Y.opt.sample[,1]),theta1.initial,0.1)
gib.theta1=gib.theta(X1,(Y1+Y.opt.sample[,1]), theta1.initial, gib.sig1, gamma1,del = 0.1,nsplit =50)
gib.gamma1=gib.gamma(n,len_theta1,gib.sig1,q_n1,gib.theta1)
gib.psi1=gib.theta1[(len_beta1+1):(len_theta1)]
A1.opt.sample[,1]=sign(H1%*%gib.psi1)
sig1.sample[1]=gib.sig1
theta1.sample[,1]=gib.theta1
gamma1.sample[,1]=gib.gamma1
for(brn in 2:nburn)
{
        sig1.sample[brn]=gib.sig(n,len_theta1,sig1.sample[brn-1],gamma1.sample[,brn-1],X1,(Y1+Y.opt.sample[,brn-1]),theta1.sample[,brn-1],del = 0.1)
        theta1.sample[,brn]=gib.theta(X1, (Y1+Y.opt.sample[,brn-1]), theta1.sample[,brn-1], sig1.sample[brn], gamma1.sample[,brn-1],del = 0.1,nsplit =50)
        gamma1.sample[,brn]=gib.gamma(n,len_theta1,sig1.sample[brn],q_n1,theta1.sample[,brn])
        gib.psi1=theta1.sample[(len_beta1+1):(len_theta1),brn]
        A1.opt.sample[,brn]=sign(H1%*%gib.psi1)
        
}
#################### Summary ###############################
#####################################################
##        Select active covariates at stage 2 in the k-th replication ##
#####################################################
gamma2.final=rep(0,len_gamma2)
for(i in 1:length(gamma2.sample[,1]))
{
    gamma2.sample0=length(which(gamma2.sample[i,1001:nburn]==0))
    gamma2.sample1=length(which(gamma2.sample[i,1001:nburn]==1))
    if(gamma2.sample1>=gamma2.sample0){
    	    gamma2.final[i]=1
    	  }else{
    	    gamma2.final[i]=0
    	  }
}
gamma2.final1=which(gamma2.final==1)
####################################################
##     Select active covariates at stage 1 in the k-th replication ##
####################################################
gamma1.final=rep(0,len_gamma1)
for(i in 1:length(gamma1.sample[,1]))                       
{
    gamma1.sample0=length(which(gamma1.sample[i,1001:nburn]==0))
    gamma1.sample1=length(which(gamma1.sample[i,1001:nburn]==1))
    if(gamma1.sample1>=gamma1.sample0){
    	   gamma1.final[i]=1
    	  }else{
    	   gamma1.final[i]=0
    	  }
}
gamma1.final1=which(gamma1.final==1)

sig1.hat=mean(sig1.sample[1001:nburn])
sig2.hat=mean(sig2.sample[1001:nburn])
return(list(gamma2.final1=gamma2.final1,gamma1.final1=gamma1.final1,gamma1.sample=gamma1.sample,gamma2.sample=gamma2.sample,theta1.sample=theta1.sample,theta2.sample=theta2.sample,sig1.sample=sig1.sample,sig2.sample=sig2.sample,sig1.hat=sig1.hat,sig2.hat=sig2.hat,Y.opt.sample=Y.opt.sample))
}

