## Sous fonctions automatisant les

penal.enet.quad <- function(lambda1=0,lambda2=0,beta=NULL){
    NonZero <- which(beta!=0)

    beta <- beta[NonZero]
    penal <- (lambda1 / abs(beta)) + lambda2 ## quadrupen naive enet
    if(length(beta) == 1){
        L <-matrix(penal)
    } else {
        L <- diag(penal)
    }
    return(list(beta=beta,L=L,NonZero=NonZero))
}

penal.enet.glm <- function(lambda=0,alpha=0,beta=NULL,n=0){
    NonZero <- which(beta!=0)

    beta <- beta[NonZero]
    penal <- lambda *n * ( ( alpha /( abs(beta))) + ((1 - alpha)))
    if(length(beta) == 1){
        L <-matrix(penal)
    } else {
        L <- diag(penal)
    }
    return(list(beta=beta,L=L,NonZero=NonZero))
}

penal.sgl <- function(lambda=0,alpha=0,beta=NULL,group=NULL,n=0){
    ##group vecteur de taille p
    #  penal <- lambda*n* ( alpha / abs(beta) + (1 - alpha) / (sqrt(rowSums(group) * rowSums(t(beta*group)^2))) )
    G <- matrix(0,length(beta),length(beta))
    for(i in 1:length(group))
        G[i,which(group == group[i])] <- 1
    NonZero <- which(beta!=0)
    wk <- sqrt(rowSums(G))[NonZero]
    beta <- beta[NonZero]
    G <- G[NonZero,NonZero]
    print(alpha)
    penal <- lambda*n*( alpha / abs(beta) + (1 - alpha) * wk / sqrt(rowSums(t(beta^2*G))))
    if(length(beta) == 1){
        L <-matrix(penal)
    } else {
        L <- diag(penal)
    }
    return(list(beta=beta,group=group[NonZero],L=L,NonZero=NonZero))
}

penal.grplasso <- function(lambda=0,beta=NULL,group=NULL,wk=NULL){
    ##group vecteur de taille p
    #  penal <- lambda*n* ( alpha / abs(beta) + (1 - alpha) / (sqrt(rowSums(group) * rowSums(t(beta*group)^2))) )
    G <- matrix(0,length(beta),length(beta))
    for(i in 1:length(group))
        G[i,which(group == group[i])] <- 1
    NonZero <- which(beta!=0)
    if(is.null(wk)){
        wk <- sqrt(rowSums(G))[NonZero]
    }
    beta <- beta[NonZero]
    G <- G[NonZero,NonZero]

    penal <- lambda/2*( wk / sqrt(rowSums(t(beta^2*G))))
    if(length(beta) == 1){
        L <-matrix(penal)
    } else {
        L <- diag(penal)
    }
    return(list(beta=beta,group=group[NonZero],L=L,NonZero=NonZero))
}

