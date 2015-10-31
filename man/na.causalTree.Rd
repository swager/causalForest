\name{na.causalTree}
\alias{na.causalTree}
\title{
  Handles Missing Values in an Rpart Object 
}
\usage{
na.causalTree(x)
}
\description{
  Handles missing values in an \code{"causalTree"} object.
}
\arguments{
\item{x}{
  a model frame.
}}
\details{
  Default function that handles missing values when calling the
  function \code{causalTree}.

  It omits cases where part of the response is missing or all the
  explanatory variables are missing.
}
\keyword{tree}
