cleaning.ridge <- function(X,y,screening=NULL,lambda1=0,lambda2=0,alpha=1,beta=rep(1,ncol(X)),nbIter=500,group=1:ncol(X),penalty="quadrupen",scale=FALSE,center=TRUE,lawWithOne=FALSE){
    refX <- NULL
    refY <- NULL

    if(!is.null(screening) && class(screening) == "screening"){
        lambda1 <- screening$lambda1.min
        lambda2 <- screening$lambda2.min
        beta <- screening$beta.min
        alpha <- NULL
        penalty <- "quadrupen"
        refY <- screening$ref$y
        refX <- screening$ref$X
    }
    if(!is.null(screening) && class(screening) == "screening.group"){
        lambda1 <- screening$lambda1.min
        lambda2 <- NULL
        beta <- screening$beta.min
        alpha <- screening$alpha.min
        penalty <- "SGL"
        refY <- screening$ref$y
        refX <- screening$ref$X
    }

    y <- scale.ref(y,scale=scale,center=center,ref=refY)
    X <- scale.ref(X,scale=scale,center=center,ref=refX)

    beta.ridge <- rep(0,ncol(X))
    Nzero <- which(beta != 0)
    beta <- beta[Nzero]
    group<-group[Nzero]
    X <- X[,Nzero]
    Pvalue <- rep(0,length(unique(group)))
    Pvalue.gamma <- rep(0,length(unique(group)))
    out <- list()

    Fstat <- rep(0,length(unique(group)))
    Fstat.permut <- matrix(0,ncol=length(unique(group)),nrow=nbIter)
    ##comput ditrib t-stat
    res <- ridge.permut(X,y,lambda1=lambda1,beta=beta,iter=nbIter,alpha=alpha,lambda2=lambda2,penalty=penalty,group=group,lawWithOne=lawWithOne)
    Pvalue<-res$pval
    Pvalue.gamma<-res$pval.gamma
    Fstat<-res$Fstat[1,]
    Fstat.permut<-matrix(res$Fstat[-1,],ncol=length(unique(group)))
    beta.ridge[Nzero] <- res$beta
    names(Fstat) <- res$nameGr
    colnames(Fstat.permut) <- res$nameGr
    out$Fstat <- Fstat
    out$Fstat.permut <- Fstat.permut



    names(Pvalue) <- unique(group)
    names(Pvalue.gamma) <- unique(group)

    out$pval <- Pvalue
    out$pval.gamma <- Pvalue.gamma
    out$beta <- beta.ridge
    out$group<-group
    out$H.Omega <- res$H.Omega
    out$H.omega <- res$H.omega
    class(out) <- "cleaning"
    return(out)
}
