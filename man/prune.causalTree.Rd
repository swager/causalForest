\name{prune.causalTree}
\alias{prune.causalTree}
\alias{prune}
\title{
  Cost-complexity Pruning of an Rpart Object 
}
\description{
Determines a nested sequence of subtrees of the supplied \code{causalTree} object 
by recursively \code{snipping} off the least important splits, based on the 
complexity parameter (\code{cp}). 
}
\usage{
prune(tree, \dots)

\method{prune}{causalTree}(tree, cp, \dots)
}
\arguments{
\item{tree}{
fitted model object of class \code{"causalTree"}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{causalTree} function. 
}
\item{cp}{
  Complexity parameter to which the \code{causalTree} object will be trimmed.
}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{
A new \code{causalTree} object that is trimmed to the value \code{cp}. 
}
\seealso{
\code{\link{causalTree}} 
}
\examples{
z.auto <- causalTree(Mileage ~ Weight, car.test.frame)
zp <- prune(z.auto, cp = 0.1)
plot(zp) #plot smaller causalTree object
}
\keyword{tree}
