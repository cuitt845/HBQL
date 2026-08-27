################################################################################
# Simulation Example for the HBQL Method
#
# This script illustrates Setting 1 with p1 = 302 and p2 = 603.
# Settings 2 and 3 can be reproduced by modifying the corresponding
# model parameters specified below.
###############################################################################################
library(glmnet)
library(ncvreg)
library(mvtnorm)
library(MASS)
source("simulate_data.R")
source("theta.initial.R")
source("gib.theta.R")
source("gib.sig.R")
source("gib.gamma.R")
source("fun_select_gamma.R")
set.seed(2019)

###########################################
pi=0.5
delta1=0.5   ###can be modified for different Settings
delta2=0.5  ###can be modified for different Settings
sigma1=1    
sigma2=1    
gamma=c(1,1,1,1,1,1,1)   ###can be modified for different Settings
n=100         ### sample size n=100
p1=150       ###m1=150; can be modified for different dimension
p2=150      ###m2=150; can be modified for different dimension
KK1=KK2=max(40,log(n))
niter=5000
nburn=5000
times=1000   ###can be modified for different replication
cut=1001    ### burn in: 1000

##stage 1 dim##
len_o1=p1
len_H1=p1+1
len_beta1=p1+1
len_psi1=p1+1
len_theta1=2*p1+2
len_gamma1=2*p1+2
long_stage1=len_theta1

##stage 2 dim##
len_o2=p2
len_H2=2*p1+2
len_beta2=2*p1+2
len_psi2=p1+p2+1
len_theta2=3*p1+p2+3
len_gamma2=3*p1+p2+3
long_stage2=len_theta2

datasets2_200=matrix(rep(0,(times*n)*(p2+3*p1+4)),nrow=(times*n),ncol=(p2+3*p1+4))
datasets1_200=matrix(rep(0,(times*n)*(2*p1+2)),nrow=(times*n),ncol=(2*p1+2))
##########################################
sum_gamma2.final1=matrix(rep(0,times*50),nrow=times,ncol=50)
sum_gamma1.final1=matrix(rep(0,times*50),nrow=times,ncol=50)
#PIP.stage2.sum=matrix(rep(0,times*dim(X2)[2]),nrow=times,ncol=dim(X2)[2])
#PIP.stage1.sum=matrix(rep(0,times*dim(X1)[2]),nrow=times,ncol=dim(X1)[2])
poster.mean.theta2.sum=matrix(rep(0,times*long_stage2),nrow=times,ncol=long_stage2)
poster.mean.theta1.sum=matrix(rep(0,times*long_stage1),nrow=times,ncol=long_stage1)
#############################################################################
start=Sys.time()
for(k in 1:times)
{
    datasets0=simulate_data(pi,delta1,delta2,sigma1,sigma2,gamma,n,p1,p2)
    O1=datasets0$O1
    A1=datasets0$A1
    O2=datasets0$O2
    A2=datasets0$A2
    Y2=datasets0$Y2
    S1=cbind(rep(1,n),O1) 
    H1=cbind(rep(1,n),O1) 
    S2=cbind(rep(1,n), O1, A1, O1*A1)  
    H2 = cbind(O1, A1, O2) 
    X1 = cbind(S1, H1*A1)     
    X2 = cbind(S2, H2*A2)    
    datasets1_200[(1+(k-1)*n):(n+(k-1)*n),]=X1
    datasets2_200[(1+(k-1)*n):(n+(k-1)*n),1:(p2+3*p1+3)]=X2
    datasets2_200[(1+(k-1)*n):(n+(k-1)*n),p2+3*p1+4]=Y2

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
########### HBQL ###############################
    select_gamma_1_2=fun_select_gamma(datasets0,KK2,KK1,niter,nburn)
    gamma2.final1=select_gamma_1_2$gamma2.final1
    gamma1.final1=select_gamma_1_2$gamma1.final1
    
    gamma2.sample.all=select_gamma_1_2$gamma2.sample[,cut:nburn]
    gamma1.sample.all=select_gamma_1_2$gamma1.sample[,cut:nburn]
## stage 2 ##
    if(length(gamma2.final1)==0 )
    {  
         sum_gamma2.final1[k,]=0
         ##PIP.stage2.sum[k,]=0
	     poster.mean.theta2.sum[k,]=0
    }else{
          if(length(gamma2.final1)==1)
          {      
   	          sum_gamma2.final1[k,1]=gamma2.final1
   	          ##PIP_stage2=mean(gamma2.sample.all[gamma2.final1,])
              theta2.sample.all=select_gamma_1_2$theta2.sample[,cut:nburn]
              poster_mean_theta2=mean(theta2.sample.all[gamma2.final1,])         
               ##PIP.stage2.sum[k,1]=PIP_stage2
	           poster.mean.theta2.sum[k,1]=poster_mean_theta2
           }else{
               sum_gamma2.final1[k,1:length(gamma2.final1)]=gamma2.final1
               ##PIP_stage2=rowMeans(gamma2.sample.all[gamma2.final1,])
               theta2.sample.all=select_gamma_1_2$theta2.sample[,cut:nburn]
               poster_mean_theta2=rowMeans(theta2.sample.all[gamma2.final1,])
               for(j in 1:length(gamma2.final1))
               {
          	          ##PIP.stage2.sum[k,gamma2.final1[j]]=PIP_stage2[j]
	                  poster.mean.theta2.sum[k,gamma2.final1[j]]=poster_mean_theta2[j]
               }
           }
    }

## stage 1 ##
    if(length(gamma1.final1)==0 )
    {  
       sum_gamma1.final1[k,]=0
       ##PIP.stage1.sum[k,]=0
       poster.mean.theta1.sum[k,]=0
     }else{
	        if(length(gamma1.final1)==1)
	        {
		       sum_gamma1.final1[k,1]=gamma1.final1
		       ##PIP_stage1=mean(gamma1.sample.all[gamma1.final1,])
               theta1.sample.all=select_gamma_1_2$theta1.sample[,cut:nburn]
               poster_mean_theta1=mean(theta1.sample.all[gamma1.final1,])    
               ##PIP.stage1.sum[k,1]=PIP_stage1
               poster.mean.theta1.sum[k,1]=poster_mean_theta1	
	        }else{
                sum_gamma1.final1[k,1:length(gamma1.final1)]=gamma1.final1
                ##PIP_stage1=rowMeans(gamma1.sample.all[gamma1.final1,])
                theta1.sample.all=select_gamma_1_2$theta1.sample[,cut:nburn]
                poster_mean_theta1=rowMeans(theta1.sample.all[gamma1.final1,])
                for(jj in 1:length(gamma1.final1))
                { 
                     ##PIP.stage1.sum[k,gamma1.final1[jj]]=PIP_stage1[jj]
                     poster.mean.theta1.sum[k,gamma1.final1[jj]]=poster_mean_theta1[jj]
                 }
           }
    }
}
end <- Sys.time()
cat("Running time:", difftime(end, start, units = "secs"), "seconds\n")

write.table(datasets2_200,"datasets2_200.txt")
write.table(datasets1_200,"datasets1_200.txt")

write.table(sum_gamma2.final1,"sum_gamma2.final1.txt")
##write.table(PIP.stage2.sum,"PIP.stage2.sum.txt")
write.table(poster.mean.theta2.sum,"poster.mean.theta2.sum.txt")


write.table(sum_gamma1.final1,"sum_gamma1.final1.txt")
##write.table(PIP.stage1.sum,"PIP.stage1.sum.txt")
write.table(poster.mean.theta1.sum,"poster.mean.theta1.sum.txt")



