ridge<-function(X,y,lambda1=0,lambda2=0,alpha=1,beta=rep(1,ncol(X)),group=1:ncol(X),penalty="quadrupen"){
    y <- as.vector(y)
    n <- length(y)
    p <- ncol(X)


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


    if(length(NonZero) == 0){
        return(list(beta.ridge=(rep(0,p)),pval=NULL))
    }

    X <- X[,NonZero]
    if(length(NonZero) == 1){
        X <- matrix(X,ncol=1)
    }
    XtX <- crossprod(X)

                                        #    Calcul Ridge

    Xy<-crossprod(X,y)


    Ainv <- diag(1/diag(L),p,p)
    B <- XtX
    T <- Ainv - Ainv%*%B%*%solve(diag(1,p) + Ainv%*%B,Ainv)

    beta.ridge<- T%*%Xy
    beta <- rep(0,p)
    beta[NonZero] <- beta.ridge

    out <- list()
    out$group <- group
    out$beta <- beta
    out$y.hat <- (X%*%beta.ridge)
    return(out)
}

ridge.svd<-function(X,y,lambda=0,beta=rep(1,ncol(X))){
    y <- as.vector(y)
    n <- length(y)
    if(length(beta) == 0){
        print("ERROR : no variables")
        return()
    }
    X<-matrix(X,ncol=length(beta))
    p <- ncol(X)
    n.lambda <- length(lambda)
    #y.tilde <- y - mean(y)
    out <- list()
    NonZero <- which(beta!=0)
    out$beta <- matrix(0,p,n.lambda)

    if(length(NonZero) == 0){
        out$y.hat <- (X%*%out$beta)
        return(out)
    }

    y.tilde <- y
    beta0 <- mean(y)

    X <- matrix(X[,NonZero],ncol=length(NonZero))
    beta <- as.vector(1/abs(beta[NonZero]))
    D<- diag(beta^(-1/2),ncol=length(beta),nrow=length(beta))
    x <- X%*%D
    SVD <- svd(x)
    v <- SVD$v
    v <- cbind(v,matrix(0,nrow(v),max(0,nrow(v)-ncol(v))))
    u <- SVD$u
    u2<-cbind(u,matrix(0,nrow(u),max(0,ncol(x)-nrow(u))))
    d <- SVD$d
    d<-c(d,rep(0,max(0,ncol(x)-length(d))))

    Delta <- rep(d,n.lambda)/(rep(d^2,n.lambda)+rep(lambda,each=ncol(x)))
    beta.ridge  <- v %*% matrix((rep(t(u) %*% y.tilde, n.lambda) * Delta),nrow=ncol(x))
    beta.ridge <- D%*%beta.ridge


    out$beta[NonZero,] <- beta.ridge
    out$y.hat <- (X%*%beta.ridge)
    return(out)
}


