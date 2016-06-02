\name{cleaning.lm}
\Rdversion{1.1}
\alias{cleaning.lm}
\title{
Cleaning.lm}
\description{
Cleaning.lm}
\usage{
cleaning.lm(X,y,screening=NULL,beta=rep(1,ncol(X)),test="student",
group=1:ncol(X),scale=FALSE,center=TRUE)
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
    not NULL then the beta.min vector contained in this object will be used instead the
    beta main parameter.
  }
  \item{beta}{
    a vector of size p representing the weight of variable. A variable with weight
    equal to zero signify will not appear in the ordinary least square
    regression. The number of variables with a non-zero weight must be
    inferious to the sample size of X.
  }
  \item{test}{
    a string of characters to determine the statistical test to be
    used. This should be one of "student" or "fisher". Student or fisher
    test are the same when variables are not grouped.  If variables are
    grouped then the student t-test measure the significancy of each
    variables and the fisher F-test measure the significancy of each
    groups of variables. Default is "student".
  }
  \item{group}{
    a vector of size p representing the group index for each
    variables. Default is 1:ncol(X) which represent the particular case
    whit no group. Warning : variables and this vector must be in
    ascending group order.
  }
  \item{scale}{
    a boolean to indicate if X and y will be scaled. Default is
    FALSE.
  }
  \item{center}{
    a boolean to indicate if X and y will be centered. Default
    is TRUE.
  }
}
\value{
  An object of class 'cleaning' with pvalue for each variables or group
  variables corresponding to the satistical test choice. Contain
  equally estimates of ols and t or F statistics.
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

