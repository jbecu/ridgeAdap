screening.group <- function(X,y,lambda = NULL,alpha=0,nlam=20,nfolds=10,group=1:ncol(X),thresh=0.001,scale=TRUE,center=TRUE){
    X <- scale(X,scale=scale,center=center)
    y <- scale(y,scale=scale,center=center)
    data = list(x = X, y = y)
    cvm <- -1
    for(a in alpha){
        CV<-cvSGL(data, index=group, maxit = 1000, thresh = thresh, min.frac = 0.05, nlam = nlam, nfold = nfolds, standardize = FALSE, verbose = FALSE, alpha = a, lambdas = lambda)

        if(cvm == -1 || cvm > min(CV$lldiff) ){
            cvm <- min(CV$lldiff)
            index.lambda.min <- which.min(CV$lldiff)
            lambda.min <- CV$lambdas[index.lambda.min]
            beta.min <- CV$fit$beta[,index.lambda.min]
            alpha.min<-a
        }
    }
    Mscreen <- list()
    Mscreen$ref <- list()
    Mscreen$ref$X <- X
    Mscreen$ref$y <- y
    Mscreen$beta.min <- beta.min
    Mscreen$S <- which(Mscreen$beta.min != 0)
    Mscreen$lambda1.min <- lambda.min
    Mscreen$alpha.min <- alpha.min
    class(Mscreen) <- "screening.group"
    return(Mscreen)
}
