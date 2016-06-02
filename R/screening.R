screening <- function(X,y,lambda1 = NULL,lambda2=0,nfolds=10,scale=TRUE,center=TRUE,se1=FALSE){

    if(nfolds <= 2) {
        print("ERROR : nfolds must greather than 2")
        return(NULL)
    } else {
        X<-scale(X,scale=scale,center=center)
        y<-scale(y,scale=scale,center=center)
    }

    CV <- crossval(X,y,intercept=FALSE,normalize=FALSE,naive=TRUE,lambda1=lambda1,lambda2=lambda2,K=nfolds,verbose=FALSE,mc.cores=1)

    if (se1 == FALSE) {
        beta.min <- CV@beta.min
        lambda1.min <- CV@lambda1.min
        lambda2.min <- CV@lambda2.min
    } else {
        beta.min <- CV@beta.1se
        lambda1.min <- CV@lambda1.1se
        lambda2.min <- CV@lambda2.min
    }

    Mscreen <- list()
    Mscreen$ref <- list()
    Mscreen$ref$X <- X
    Mscreen$ref$y <- y
    Mscreen$beta.min <- beta.min
    Mscreen$S <- which(Mscreen$beta.min != 0)
    Mscreen$lambda1.min <- lambda1.min
    Mscreen$lambda2.min <- lambda2.min
    class(Mscreen) <- "screening"
    return(Mscreen)
}
