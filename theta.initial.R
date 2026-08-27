#######################################
#### Function for initializing theta             #####
#######################################

theta.initial=function(X,Y,len_theta)   
{
   theta_initial=rep(0,len_theta)
   for(i in 1:length(X[1,]))
   {
       XtX_i.inv=solve(t(X[,i])%*%X[,i])
       theta_initial[i]=XtX_i.inv%*%t(X[,i])%*%Y
   }
   return(theta_initial)
}
