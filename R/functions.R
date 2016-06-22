
                                      # SUBFUNCTIONS

paramGammaMoment <- function(x){ ## gamma approximation by the moment method

    m1 <- mean(x)
    m2 <- mean(x^2)

    k <- m1^2 / ((m2 - m1^2))

    Theta <-  (m2 - m1^2)/ m1

    return(list(k = k, theta =Theta))
}

tr <- function(x){ ## trace function
    return(sum(diag(x)))
}

scale.ref<- function(x,scale=FALSE,center=TRUE,ref=NULL){
    # TRUE/FALSE : scale/center or not on x himself
    # ref : object issued from scale or scale.ref
    if(is.vector(x)){
        x <- matrix(x,ncol=1)
    }
    if(is.null(ref) || length(attributes(ref)) == 1){
        return(scale(x,scale=scale,center=center))
    } else {
        if(scale == TRUE)
            scale <- attributes(ref)[["scaled:scale"]]
        if(center == TRUE)
            center <- attributes(ref)[["scaled:center"]]
        return(scale(x,center=center,scale=scale))
    }
}

# SUBFUNCTIONS FOR RIDGE

ridge.permut<-function(X,y,lambda1=0,lambda2=0,alpha=1,beta=rep(1,ncol(X)),iter=500,group=1:ncol(X),penalty="quadrupen",lawWithOne=FALSE){
    y <- as.vector(y)
    n <- length(y)


    ## Calculation of adaptive penalty
    if(penalty == "glmnet"){
        pen <- penal.enet.glm(lambda=lambda1,alpha=alpha,beta=beta,n=n)
        group <- group[pen$NonZero]
    }
    if(penalty == "quadrupen"){
        pen <- penal.enet.quad(lambda1=lambda1,lambda2=lambda2,beta=beta)
        pen$NonZero
        group <- group[pen$NonZero]
    }
    if(penalty == "SGL"){
        pen <- penal.sgl(lambda=lambda1,alpha=alpha,beta=beta,n=n,group=group)
        group <- pen$group
    }
    if(penalty == "grplasso"){
        pen <- penal.grplasso(lambda=lambda1,beta=beta,group=group)
        group <- pen$group
    }
    L <- pen$L
    beta <- pen$beta
    NonZero<-pen$NonZero


    p <- length(beta)
    perm <- matrix(1:n,nrow=n,ncol=iter)
    perm<-apply(perm,2,sample)
    if(p == 1){
        X <- matrix(X,ncol=1)
    }
    ## to avoid problemto invert t(X)%*%X + L we use the woodbury formula
    ## where (A+B)^-1 = A^-1 - A^-1B(I_p + A^-1B)^-1A^-1
    Ainv <- diag(1/diag(L),p,p)
    B <- crossprod(X)
    T <- Ainv - Ainv%*%B%*%solve(diag(1,p) + Ainv%*%B,Ainv)
#    T <- solve(t(X)%*%X + L)
    H.Omega <- X%*%T%*%t(X)
    Xty <- crossprod(X,y)
    beta.orig <- T%*%Xty
    yhat<-matrix(0,ncol=iter+1,nrow=n)
    yhat[,1] <- X%*%beta.orig ## y.hat estimate without permutation, the same for all groups

    ngroup <- length(unique(group))
    Fstat <- matrix(0,ncol=ngroup,nrow=iter+1)

    index <- unique(group)
    pval <- rep(0,ngroup)
    pval.gamma <- rep(0,ngroup)
    H.omega <- list()
    for(j in 1:ngroup){ #@ For each we need to estim y.hat to build our F-statistic (current group is noticed by j)
        sizeG <- sum(group==index[j])
        isJ <- which(group==index[j]) ### Index of variables contained in the jth group
        isNotJ <- which(group!=index[j]) ### Index of variables contained out the jth group
        xj<- matrix(X[,isJ],ncol=sizeG)
        x<- X[,isNotJ]

        Xj <- matrix(0,nrow=n,ncol=iter*sizeG)


        for(i in 1:iter){  # Generation of all permuted version of variables contained in the jth group
          Xj[,((i-1)*sizeG+1):(sizeG+ (i-1)*sizeG)] <- xj[perm[,i],]
        }
        a<-crossprod(xj) + L[isJ,isJ] ## initialisation of A
        Xjty<- crossprod(Xj,y) ## crossprod of Y with all permuted variables and groups
        if(ngroup == 1){ ## Special case if only one group exists
            yhat.omega <- matrix(0,nrow=n)
            A <- solve(matrix(a,ncol=sizeG)) ## A is the same for each permutations for this special case
            for(i in 1:iter){
                beta <- A %*% matrix(Xjty[(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG)],ncol=1)
                yhat[,(i+1)] <- Xj[,(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG)] %*% beta
            }
        } else {  ## General case for multiple groups
            if(sizeG == 1) ## calcul of D for group at size 1
                D <- T[isNotJ,isNotJ]-  (matrix(T[isNotJ,isJ],ncol=sizeG)%*%matrix(1/T[isJ,isJ],ncol=length(isJ))%*%matrix(T[isJ,isNotJ],nrow=sizeG))
            else ## calcul of D for group at size upper to 1
                D <- T[isNotJ,isNotJ] -  (T[isNotJ,isJ]%*%solve(T[isJ,isJ],T[isJ,isNotJ]))
            beta.omega <- D%*%Xty[isNotJ] ## coefficient of reduce model
            H.omega[[j]] <- X[,isNotJ]%*%D%*%t(X[,isNotJ]) ##
            beta<-matrix(c(rep(0,length(isJ)),beta.omega),ncol=iter,nrow=p)


            if(!lawWithOne || j == 1){ ## test to know if it needs to do permutation for the jth group, an approximation could be used when all groups have the same size
                Vj <- t(Xj)%*%x
                U <- - D%*% t(Vj)
                for(i in 1:iter){ ## Correction of beta.omega fo each permutation
                    Uj <- rbind(diag(rep(1,length(isJ))),matrix(U[,(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG)],ncol=sizeG))
                    A <- matrix(solve((a + matrix((Vj[(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG),]),nrow=sizeG)%*%Uj[-(1:length(isJ)),])),ncol=sizeG)
                    beta[,i] <- beta[,i]+Uj%*%A%*%t(Uj)%*%matrix(c(Xjty[(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG)],Xty[isNotJ]),ncol=1)
                    yhat[,(i+1)] <- matrix(X[,isNotJ],ncol=(ncol(X)-sizeG))%*%matrix(beta[-(1:sizeG),i],ncol=1) + (matrix(Xj[,(1+ (i-1)*sizeG):(sizeG+ (i-1)*sizeG)],ncol=length(isJ)) %*% beta[1:sizeG,i]) ## yhat of full (permuted) model
                }
            }
            yhat.omega <- x%*%beta.omega ## y.hat of reduce model
        }

        rss.omega <- RSS.permut(y,yhat.omega) ## RSS for reduce model
        rss.Omega<-RSS.permut(y,yhat) ## RSS for full (unpermuted and permuted model)
        Fstat.Omega <- Ftest.permut(rss.omega,rss.Omega)$F ## Fstat calculcation for (unpermuted and permuted model if exist)
        if(!lawWithOne || j == 1){
            Fn <- ecdf(Fstat.Omega[-1]) ## Distribution of Fstat under null hypothesis
            gam<-paramGammaMoment(Fstat.Omega[-1]) ## Distribution of Fstat under null hypothesis by gamma law approximation
        }
        pval[j] <- Fn(Fstat.Omega[1])  ## pval for j th group
        pval.gamma[j] <- pgamma(Fstat.Omega[1],shape=gam$k,scale=gam$theta) ## pval for j the group based on gamma law approximation
        Fstat[,j] <- Fstat.Omega
    }
    return(list(pval=1 - pval,pval.gamma= 1 - pval.gamma,Fstat=Fstat,beta=beta.orig,nameGr=index,H.Omega=H.Omega,H.omega=H.omega))
}


RSS.permut<- function(y,yhat){ ### RSS calculation
        return(colSums((y-yhat)^2))
    }

Ftest.permut<- function(RSS1,RSS2){ ### F-statistic calculation
        return(list(F = (RSS1-RSS2)/(RSS2),full = RSS2,reduce=RSS1-RSS2))
    }



