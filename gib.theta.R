###############################################
##    Function for sampling theta in the Gibbs sampler     ##
###############################################
# ##nsplit: number of coefficient blocks used in block Gibbs sampling;
# ##it may be adjusted according to the dimensions p1 and p2.
    gib.theta <- function(X, Y, B, sig, Z, del = 0.1, nsplit )
    {   p = dim(X)[2]; n = dim(X)[1]
        s0=(0.1*sig)/n
        s1  = sig* max(log(n), 0.01*p^{2 + del}/n)   
        T1 = Z*s1 + (1-Z)*s0      
        vsize = p %/%nsplit 
        G = t(X)%*%X
        vec = seq(1:vsize)
        Xtmp = X
        for(s in 1:nsplit)
        {
            svec = (s-1)*vsize +vec
            COV = (G[svec,svec] + diag(1/T1[svec]))
            ec = eigen(COV)
            COVsq = ec$vectors %*% diag(1/sqrt(ec$values)) %*% t(ec$vectors)
            B[svec] = COVsq %*% ( COVsq %*% (t(X[,svec]) %*% Y - G[svec, -svec] %*% B[-svec]) + sqrt(sig)*rnorm(vsize)  )    #block updating 
        }
        if( p > nsplit*vsize)
        {
            svec = (nsplit*vsize +1):p
            COV = (G[svec,svec] + diag(1/T1[svec]))
            ec = eigen(COV)
            COVsq = ec$vectors %*% diag(1/sqrt(ec$values)) %*% t(ec$vectors)
            B[svec] = COVsq %*% (COVsq %*% (t(X[,svec]) %*% Y - G[svec, -svec] %*% B[-svec]) + sqrt(sig)*rnorm(length(svec)))
        }
        return(B)
    }

