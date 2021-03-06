\name{ridge}
\Rdversion{1.1}
\alias{ridge}
\title{
Ridge}
\description{
Ridge}
\usage{
ridge(X,y,lambda1=0,lambda2=0,alpha=1,beta=rep(1,ncol(X)),group=1:ncol(X),
penalty="quadrupen")
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
  \item{group}{
    a vector of size p representing the group index for each
        variables. Default is 1:ncol(X) which represent the particular
	case
	whit no group. Warning : this vector must be in
	ascending group order and variables ordered in this sense.
	  }
  \item{penalty}{
    a string of characters to determine the adaptive ridge specific
    penalty  to be used. This should be one of "quadrupen" to
    elastic-net by quadrupen package, "glmnet" to elastic-net by
        glmnet package, "grplasso" to group lasso by grplasso package
	and "SGL"
	    to sparse group lasso by SGL package at screening. By
	    default is "quadrupen".
  }
}
\value{
An list with ridge estimates, y estimates.
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

