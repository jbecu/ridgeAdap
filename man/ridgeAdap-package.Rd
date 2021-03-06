\name{ridgeAdap-package}
\Rdversion{1.1}
\alias{ridgeAdap-package}
\alias{ridgeAdap}
\docType{package}
\title{
Adaptive ridge
}
\description{
  This package is dedicated to feature selection based on the adaptive ridge method.
}
\details{

  The key idea of the adaptive ridge regression is the use of specific
  penalty into the ridge regression that could mimic lasso and it's variation.
  The equivalence between ridge regression with specific penalty and
  lasso show that this penalty is the better way to transfert signal
  from lasso to ridge.

  The adaptive ridge constitute same at the adaptive lasso, to learn
  elastic-net or sparseg group lasso estimates in a first step, and use the estimates of lasso to
  define a specific penalty to the ridge regression. This penalty definition differs 
  if elastic-net or sparse group lasso was used in the first step, that
  explain the name of "adaptive" ridge regression.

Variable selection is done with a step protocol based on Larry Wasserman (2009) works.

Af first step (\code{\link{screening},\link{screening.group}}), to select variables by a sparse approach on the half of
the dataset. This package use the Elastic-net regression, which is write
as follows :
\if{latex}{\deqn{%
      \hat{\beta}_{\lambda_1,\lambda_2} = \arg \min_{\beta}
     RSS(beta) + \lambda_1 \|
      \beta \|_1 + \lambda_2  \|\beta\|_2. }}
 \if{text}{\deqn{beta.hat(lambda1, lambda2) =
      argmin_beta RSS(beta) + lambda1 || beta||_1 + lambda2
      || beta||_2.}}

A second step (\code{\link{cleaning.ridge}}), to test variables in the subset
of varibles selected in the first step. This step is applied on the
second half of the dataset with the adaptive ridge regression, which is
write
as follows :
 \if{latex}{\deqn{%
      \hat{\beta}_{\Omega} = \arg \min_{\beta}
     RSS(beta) + \sum{\omega_j \beta_j},}}
 \if{text}{\deqn{beta.hat(Omega) =     
      argmin_beta RSS(beta) + sum(omega_j beta_j),}}
with a specific penalty for the elastic-net
  \if{latex}{\deqn{%
      \omega_j = \frac{\lambda_1}{|\hat{\beta}_j(\lambda_1)|} + \lambda_2.}}
 \if{text}{\deqn{
      omega_j = lambda1/|beta_j(lambda1)| + lambda2.}}
The significance test for the ridge regressor is based on the
fisher's staistic and permutations.

Parameters \eqn{\lambda_1,\lambda} are chosen by k-folds
Cross-Validation on the entire process.
  
  This package is under development, all feedback are appreciate.
}
\author{
  Jean-Michel Becu.

  Maintener : Jean-Michel Becu (jean-michel.becu@inria.fr)
}

\references{
  Becu J.M., Ambroise C., Dalmasso C. and Grandvalet Y.
  Beyond Support in Two-Stage Variable Selection. \emph{Statistics and
  Computing}, \bold{26}:1-11,2016.

  For further information :

  L. Wasserman and K. Roeder. High-dimensional variable selection. \emph{The
  Annals of Statistics}, \bold{37}(5A):2178–2201, 2009.
  Wasserman

  Y. Grandvalet and S. Canu. Outcomes of the equivalence of adaptive
  ridge with least absolute shrinkage. In M. S. Kearns, S. A. Solla, and
  D. A. Cohn, editors, \emph{Advances in Neural Information Processing Systems
  11 (NIPS 1998)}, pages 445–451. MIT Press, 1999.

  }
\keyword{package}
\examples{

library("ridgeAdap")
#generate dataset with block's type
n <- 250
p <- 500
rho <- 0.5

nRel <- 25
SNR <- 4
group<-sort(rep(1:20,25))

S <-matrix(unlist(lapply(1:p,FUN=function(x){group == group[x]})),ncol=p)
sigma<-matrix(0,p,p)
sigma[S==TRUE]<-rho
diag(sigma)<-1

beta<-rep(0,p)
beta[sample(1:p,nRel)]<-runif(nRel,0.1,1)
U <- t(chol(sigma))
random.normal <- matrix(rnorm(p*n,0,1), nrow=p,ncol=n)
X <- U\%*\% random.normal
X<-t(X)

y<-as.vector(t(matrix(beta,nrow=p))\%*\%t(X))
var<-sqrt(t(matrix(beta,ncol=1))\%*\%sigma\%*\%matrix(beta,ncol=1))*(1/(SNR))
y<-y+rnorm(n,0,var)

##selection procedure
subsets <-split(1:n,1:2)
Screen<-screening(X[subsets[[1]],],y[subsets[[1]]])
Clean<-cleaning.ridge(X[subsets[[2]],],y[subsets[[2]]],screening=Screen)
selectVar<-as.integer(names(which(p.adjust(Clean$pval,"BH") <= 0.05)))

Screen.group<-screening.group(X[subsets[[1]],],y[subsets[[1]]],group=group)
Clean.group<-cleaning.ridge(X[subsets[[2]],],y[subsets[[2]]],screening=Screen.group,group=group)

selectGroup<-as.integer(names(which(p.adjust(Clean.group$pval,"BH") <= 0.05)))

real <- which(beta!=0) ## index for variables with beta != 0
realInGroup <- table(group[beta!=0]) ## Number of variables with beta!=
0 for explanatory groups

##estimation procedure

n <- 250
p <- 500
rho <- 0.5

nRel <- 25
SNR <- 4
group<-sort(rep(1:20,25))

S <-matrix(unlist(lapply(1:p,FUN=function(x){group ==
group[x]})),ncol=p)
sigma<-matrix(0,p,p)
sigma[S==TRUE]<-rho
diag(sigma)<-1

beta<-rep(0,p)
beta[sample(1:p,nRel)]<-runif(nRel,0.1,1)
U <- t(chol(sigma))
random.normal <- matrix(rnorm(p*n,0,1), nrow=p,ncol=n)
X <- U\%*\% random.normal
X<-t(X)

random.normal <- matrix(rnorm(p*n,0,1), nrow=p,ncol=n)
X.test <- U\%*\% random.normal
X.test<-t(X.test)

y<-as.vector(t(matrix(beta,nrow=p))\%*\%t(X))
var<-sqrt(t(matrix(beta,ncol=1))\%*\%sigma\%*\%matrix(beta,ncol=1))*(1/(SNR))
y<-y+rnorm(n,0,var)

esti<-estimation.AR(X,y)

esti$beta.min

}