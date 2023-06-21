\name{nlmeans}
\alias{nlmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NLMeans filter in 1D/2D/3D
}
\description{
Implements the Non-Local-Means Filter of Buades et al 2005
}
\usage{
nlmeans(x, lambda, sigma, patchhw = 1, searchhw = 7, pd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  1, 2 or 3-dimensional array of obseved response (image intensity) data.
}
  \item{lambda}{
  scale factor for kernel in image space.
}
  \item{sigma}{
  error standard deviation (for additive Gaussian errors).
}
  \item{patchhw}{
  Half width of patches in each dimension (patchsize is \code{(2*patchhw+1)^d}
  for d-dimensional array).
}
  \item{searchhw}{
  Half width of search area (size of search area is \code{(2searchhw+1)^d}
  for d-dimensional array)).
}
  \item{pd}{
  If \code{pd < (2*patchhw+1)^d} use \code{pd} principal components instead
  of complete patches.
}
}
\details{
The implementation follows the description of the Non-Local-Means Filter of
Buades et al 2005 on \url{http://www.numerical-tours.com/matlab/denoisingadv_6_nl_means/#biblio}
that incorporates dimension reduction for patch comparisons by PCA.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list of class \code{"nlmeans"} with components
\item{theta}{Denoised array}
\item{lambda}{Scale parameter used}
\item{sigma}{The error standard deviation}
\item{patchhw}{Half width of patches}
\item{pd}{Effective patchsize used}
\item{searchhw}{Half width of search area}
}
\references{
J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06 .

A. Buades, B. Coll and J. M. Morel (2006). A review of image denoising algorithms, with a new one.
Simulation, 4, 490-530. DOI:10.1137/040616024.

\url{http://www.numerical-tours.com/matlab/denoisingadv_6_nl_means/#biblio}
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}
}
\note{
use \code{setCores='number of threads'} to enable parallel execution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
