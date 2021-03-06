\name{ridge.permut}
\Rdversion{1.1}
\alias{ridge.permut}
\title{
Ridge.permut}
\description{
Ridge.permut}
\usage{
ridge.permut(X,y,lambda1=0,lambda2=0,alpha=1,beta=rep(1,ncol(X)),iter=1000,
group=1:ncol(X),penalty="quadrupen")
}
\arguments{
 \item{X}{
    an matrix object with n*p size
  }
  \item{y}{
    an vector of size p
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
    a vector of size p representing the weight of variable. A
    variable with weight
    equal to zero signify will not appear in the ordinary
    least square
    regression. No constraint to the number of variables
    with a non-zero weight.
  }
  \item{iter}{
    Number of iterations used to simulate the null hypothesis
    law. Bigger is better. Default is 1000.
  }
  \item{group}{
    a vector of size p representing the group index for each
        variables. Default is 1:ncol(X) which represent the particular
	case
	whit no group. Warning : this vector must be in
	ascending group order and variables ordered in this sense.
      }
      \item{penalty}{
	a string of characters to determine the adaptive ridge
	specific
	penalty  to be used. This should be one of "quadrupen"
	to
	elastic-net by quadrupen package, "glmnet" to elastic-net by
	glmnet package, "grplasso" to group lasso by grplasso package and
	"SGL"
        to sparse group lasso by SGL package at screening. By default is
	"quadrupen".
      }
      \item{lawWithOne}{
	    a boolean to indicate if one distribution under the null
	    hypothesis
	        for all groups of variables could be used instead of one
		distribution for
		    each of them. Warning : It is an approximation and
		    that work only
		        if all groups have the same size. Default is
			FALSE.
			    }
    }
\value{
An list with ridge estimates, F statistic an F statistics under H0 and
p-value for each variables or group significancy. Contains equally
pval.gamma obtaines with the gamma law approximation.
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
\keyword{ridge}

