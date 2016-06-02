cleaning.lm <- function(X,y,screening=NULL,beta=rep(1,ncol(X)),test="student",group=1:ncol(X),scale=FALSE,center=TRUE){

    ##test peux etre student ou fisher

    refX <- NULL
    refY <- NULL
    if(!is.null(screening)){
        refX <- screening$ref$X
        refY <- screening$ref$y
        beta <- screening$beta.min
    }
    y <- scale.ref(y,scale=scale,center=center,ref=refY)
    X <- scale.ref(X,scale=scale,center=center,ref=refX)

    S <- which(beta != 0)
    n <- length(y)
    if(length(S) > n){
        S <- as.vector(order(abs(beta),decreasing=TRUE)[1:(n-1)])
    }
    lm.res <- lm(y ~ X[,S]-1)

    out <- list()
    coef <- coefficients(summary(lm.res))
    if(test == "student"){

        pVal <- coef[,4]

        out$pval <- rep(1,length(S))
        out$T <- rep(0,length(S))
        out$sd <-rep(0,length(S))

        out$pval[is.element(S,S)] <- pVal
        out$T[is.element(S,S)] <- coef[,3]
        out$sd[is.element(S,S)] <-coef[,2]

        out$beta.min<-rep(0,ncol(X))
        out$beta.min[S] <- as.vector(coef[,1])


        names(out$pval) <- as.vector(S)
        names(out$T) <- as.vector(S)
        names(out$sd) <- as.vector(S)
        out$df <- lm.res$df.residual
        out$test <- test
    } else {
        ngroup <- length(unique(group[S]))
        index <- unique(group[S])
        pVal <- rep(0,ngroup)
        F<-pVal
        for(j in 1:ngroup){
            isG <- S[is.element(S,which(group != index[j]))]
            df2 <- lm.res$df.residual
            if(ngroup != 1){
                lm.res.little <- lm(y ~ X[,isG]-1)
                df1 <- lm.res.little$df.residual - lm.res$df.residual

                F[j] <- ((sum(lm.res.little$residuals^2)-sum(lm.res$residuals^2)) / df1)/(sum(lm.res$residuals^2)/df2)
            } else {
                df1 <- length(y)-lm.res$df.residual
                F[j] <- ((sum(y^2)-sum(lm.res$residuals^2)) / df1)/(sum(lm.res$residuals^2)/df2)
            }

            pVal[j] <- 1-pf(F[j],df1,df2)
        }

        out$pval <- pVal
        names(out$pval) <- index
        out$F <- F
        names(out$F) <- index
        out$beta.min<-rep(0,ncol(X))
        out$beta.min[S] <- as.vector(coef[,1])
        out$test <- "fisher"
    }


    class(out) <- "cleaning"
    return(out)
}
