\name{screening.group}
\Rdversion{1.1}
\alias{screening.group}
\title{
Screening.group}
\description{
Screening.group}
\usage{
screening.group(X,y,lambda = NULL,alpha=0,nlam=20,nfolds=10,group=1:ncol(X),
thresh=0.001,scale=FALSE,center=TRUE)
}
\arguments{
 \item{X}{
    a matrix object with n*p size representing the data
  }
  \item{y}{
    a vector of size p representing the output variable
  }
  \item{lambda}{
    a vector representing the l1 penalty. Default is NULL and
    automatically defined by SGL.
  }
  \item{alpha}{
    a vector representing the \eqn{\alpha} parameter for the sparse group lasso. If alpha = 1 that correspond
    to the lasso solution, if alpha = 0 that correspond to the group
    lasso and otherwise to the sparse group lasso solution. Default is 0.
  }
  \item{nlam}{
    an integer representing the size of the lambda grid if lambda is
    NULL. Default is 20.
    }
  \item{nfolds}{
    an integer representing the number of folds for the k-folds
    cross-validation used to choose l1 and l2 penatly. Default is 10 and
    must be greater than 2.
  }
  \item{group}{
    a vector of size p representing the group index for each variables.
    Default is 1:ncol(X) which represent the particular case whit no
    group. Warning this vector must be in ascending group order and
    variables in X needs to be ordered in this sense.
  }
  \item{thresh}{
    convergence threshold for change in beta.
    }
  \item{scale}{
    a boolean to indicate if X and y will be scaled. Default is TRUE.
  }
  \item{center}{
    a boolean to indicate if X and y will be centered. Default is TRUE.
  }
}

\value{
  A object of class 'screening' with 5 elements
  \item{beta.min}{
    a vector of size p with the elastic.net estimates
  }
  \item{S}{
    a vector contaning the index of variable which $B_j != 0$
  }
  \item{lambda1.min}{
    l1 penalty choose by CV
  }
  \item{alpha.min}{
    alpha mixing parameter choose by CV
  }
  \item{ref}{
    a list with two component contains respectively the X and y variance
  and mean used to scaling and/or centering the data.
}
  
}
\author{
JM BECU
}
\note{
  \code{screening.group} requires \code{SGL} package.
}
\examples{
 # see ridgeAdap help
}
\keyword{screening}

