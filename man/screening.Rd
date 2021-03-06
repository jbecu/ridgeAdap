\name{screening}
\Rdversion{1.1}
\alias{screening}
\title{
Screening}
\description{
Screening}
\usage{
screening(X,y,lambda1 = NULL,lambda2=0,nfolds=10,scale=FALSE,center=TRUE,se1=FALSE)
}
\arguments{
 \item{X}{
    a matrix object with n*p size representing the data
  }
  \item{y}{
    a vector of size p representing the output variable
  }
  \item{lambda1}{
    a vector representing the l1 penalty. Default is NULL and
    automatically defined by quadrupen.
  }
  \item{lambda2}{
    a vector representing the l2 penalty. If lambda2 = 0 that correspond
    to the lasso solution, otherwise to the elastic-net solution. Default is 0.
  }
  \item{nfolds}{
    An integer representing the number of folds for the k-folds
    cross-validation used to choose l1 and l2 penatly. Default is 10 and
    must be greater than 2.
  }
  \item{scale}{
    a boolean to indicate if X and y will be scaled. Default is FALSE.
  }
  \item{center}{
    a boolean to indicate if X and y will be centered. Default is TRUE.
  }
  \item{se1}{
    a boolean to indicate if using the se1 rule for the cross validation
    with nfolds upper than 2.
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
    l1 penalty chosen by CV
  }
  \item{lambda2.min}{
    l2 penalty chosen by CV
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
  \code{screening} requires \code{quadrupen} package.
}
\examples{
 # see ridgeAdap help
}
\keyword{screening}

