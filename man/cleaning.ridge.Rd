\name{cleaning.ridge}
\Rdversion{1.1}
\alias{cleaning.ridge}
\title{
Cleaning.ridge}
\description{
Cleaning.ridge}
\usage{
cleaning.ridge(X,y,screening=NULL,lambda1=0,lambda2=0,alpha,beta=rep(1,ncol(X)),
nbIter=1000,group=1:ncol(X),penalty="quadrupen",scale=FALSE,center=TRUE)
}
\arguments{
  \item{X}{
    a matrix object with n*p size representing the data
  }
  \item{y}{
    a vector of size p representing the output variable
  }
  \item{screening}{
    object issued from the screening or screening.group functions. If is
    not NULL then the beta.min, lambda1.min, lambda2.min, penalty type informations contained in this object will be used instead the beta main parameter. Default is NULL.
  }
  \item{lambda1}{
    The l1 penalty to used.
  }
  \item{lambda2}{
    The l2 penalty to use in case of elatic-net for screening.
  }
  \item{alpha}{
    The \eqn{\alpha} mixing parameter to use in case of sparse group
    lasso for screening or if glmnet is used for screening.
  }
  \item{beta}{
    a vector of size p representing the weight of variable. A variable with weight
    equal to zero signify will not appear in the ordinary least square
    regression. No constraint to the number of variables with a non-zero weight.
  }
  \item{nbIter}{
    Number of iterations used to simulate the null hypothesis
    law. Bigger is better. Default is 1000. Could be less with the gamma
    law approximation.
  }
  \item{group}{
    a vector of size p representing the group index for each
    variables. Default is 1:ncol(X) which represent the particular case
    whit no group. Warning : this vector must be in
    ascending group order and variables ordered in this sense.
  }
  \item{penalty}{
    a string of characters to determine the adaptive ridge specific
    penalty  to be used. This should be one of "quadrupen" to
    elastic-net by quadrupen package, "glmnet" to elastic-net by
    glmnet package, "grplasso" to group lasso by grplasso package and "SGL"
    to sparse group lasso by SGL package at screening. By default is "quadrupen".
  }
  \item{scale}{
    a boolean to indicate if X and y will be scaled. Default is
    FALSE.
  }
  \item{center}{
    a boolean to indicate if X and y will be centered. Default
    is TRUE.
  }
  \item{lawWithOne}{
    a boolean to indicate if one distribution under the null hypothesis
    for all groups of variables could be used instead of one distribution for
    each of them. Warning : It is an approximation and that work only
    if all groups have the same size. Default is FALSE. 
    }
}
\value{
  An object of class 'cleaning' with pvalue for each variables or group
  variables. Two pvalues are given, the first is estimates on the real
  distribution of F-statistic under the null hypothsesis, the second
  pval.gamma is estimates on the gamma low fitting on the real
    distribution of F-statistic. This approximation could be used to
    minimize the number of needed iteration. Contains equally estimates of adaptive ridge, F statistics
  on estimates and F statistics simulated under the H0 model.
}
\author{
JM BECU
}
%\note{
%  \code{screening} requires \code{glmnet} package.
%}
\examples{
 # see ridgeAdap help
}
\keyword{cleaning}

