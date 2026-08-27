#########################################################
##     Function for the conditional posterior distribution of  gamma     ##
#########################################################

gib.gamma=function(n_person,len_theta,sig_posterior,q_n,gib_theta)
{
        s1  =  sig_posterior*max(0.01*len_theta^{2.1}/n_person,log(n_person))      
        s0=(0.1*sig_posterior)/n_person
        s = seq(1:len_theta)
        prob = sapply(s, function(j) q_n* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s1))/ (q_n* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s1))+ (1-q_n)* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s0))))
       nan_prob=which(is.na(prob)==TRUE)
       if(length(nan_prob)==0)
       {
            tmp = (runif(len_theta) - prob)
            gamma=rep(3,length(prob))
        for(i in 1:length(prob))
        {
            if(tmp[i]>0){ 
            	    gamma[i]=0 
            	}else{ 
            		gamma[i]=1 
            	}
            
        }
        if(sum(gamma) > n_person/2)    # if the selected model is too large 
        {
            indz = which(gamma == 1)
            gamma[indz[which(prob[gamma] < rev(sort(prob[gamma]))[round(n_person/2)])]] = 0 
        }
        return(gamma)
      }else{
        prob = sapply(s, function(j) q_n* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s1))/ (q_n* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s1))+ (1-q_n)* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s0))))
        prob[nan_prob]= sapply(nan_prob, function(j) q_n* dnorm(gib_theta[j], 0, sqrt(sig_posterior*s1))/ (1e-8))
        tmp = (runif(len_theta) - prob)
        gamma=rep(3,length(prob))
        for(i in 1:length(prob))
        {
            if(tmp[i]>0){ 
            	gamma[i]=0 
            	}else{ 
            	gamma[i]=1 
            	}    
        }
        if(sum(gamma) > n_person/2)    # if the selected model is too large 
        {
            indz = which(gamma == 1)
            gamma[indz[which(prob[gamma] < rev(sort(prob[gamma]))[round(n_person/2)])]] = 0 
        }
        return(gamma)
      } 
}