\name{causalTree-internal}
\alias{pred.causalTree}
\alias{causalTree.matrix}
\title{
  Internal Functions 
}
\description{
  Internal functions, only used by packages \pkg{causalTree} and \pkg{ipred}.
}
\usage{
pred.causalTree(fit, x)
causalTree.anova(y, offset, parms, wt)
causalTree.class(y, offset, parms, wt)
causalTree.matrix(frame)
causalTree.poisson(y, offset, parms, wt)
causalTreeco(tree, parms)
}
\arguments{
  \item{fit}{a tree fitted by \code{causalTree}.}
  \item{x}{a matrix of predictors.}
  \item{y}{the responses.}
  \item{offset}{an offset, or \code{NULL}.}
  \item{parms}{a list of parameters, usually empty.}
  \item{wt}{case weights.}
  \item{frame}{model frame (from call to \code{causalTree})}
  \item{tree}{a tree fitted by \code{causalTree}.}
}
\value{
  For \code{causalTreeco} the x,y plotting coordinates of the nodes.
  
  \code{causalTree.anova}, \code{causalTree.class} and
  \code{causalTree.poisson} return a list with components
  \item{y}{(adjusting for \code{offset} if necessary),}
  \item{parms}{as input,}
  \item{numresp}{the number of responses,}
  \item{summary}{a function to be invoked by \code{\link{summary.causalTree}},}
  \item{text}{a function to be invoked by \code{\link{text.causalTree}}.}
}
\keyword{internal}
