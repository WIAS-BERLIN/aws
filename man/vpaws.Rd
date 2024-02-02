\name{vpaws}
\alias{vpaws}
\alias{vpawscov}
\alias{vpawscov2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
vector valued version of function \code{\link{paws}} with homogeneous covariance structure
}
\description{The function implements a vector-valued version the propagation separation approach that
uses patches instead of individuel voxels for comparisons in parameter space. Functionality is analog to function \code{\link{vaws}}. Using patches allows for an improved
handling of locally smooth functions and in 2D and 3D for improved smoothness of
discontinuities at the expense of increased computing time.

}
\usage{
vpaws(y, kstar = 16, sigma2 = 1, invcov = NULL, mask = NULL, scorr = 0, spmin = 0.25,
      ladjust = 1, wghts = NULL, u = NULL, patchsize = 1)
vpawscov(y, kstar = 16, invcov = NULL, mask = NULL, scorr = 0, spmin = 0.25, ladjust = 1, 
      wghts = NULL, maxni = FALSE, patchsize = 1)
vpawscov2(y, kstar = 16, invcov = NULL, mask = NULL, scorr = 0, spmin = 0.25,
      lambda = NULL, ladjust = 1, wghts = NULL, patchsize = 1,
      data = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{\code{y}  can be a full array of vector valued data, or, if \code{mask} is provided,
   be a matrix with columns corresponding to points/pixel/voxel within the mask. In the first case
  \code{dim(y)} determines the dimensionality and extend of the grid design, in the second
  case tis information is obtained from the dimensions of \code{mask}.
  the first component varies over components of the response vector.
}
  \item{kstar}{
maximal number of steps to employ. Determines maximal bandwidth.
}
  \item{sigma2}{
specifies a homogeneous error variance.
}
\item{invcov}{
 array (or matrix) of voxelwise inverse covariance matrixes, first index corresponds to
 upper diagonal inverse covariance matrix.
}
  \item{mask}{
logical mask. All computations are restrikted to design poins within the mask.
}
  \item{scorr}{
The vector \code{scorr} allows to specify a first order correlations of the noise for each coordinate direction,
    defaults to 0 (no correlation).
}
  \item{spmin}{
determines the form (size of the plateau) in the adaptation kernel.
    Not to be changed by the user.
}
  \item{ladjust}{
factor to increase the default value of lambda
}
  \item{wghts}{\code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space.
}
  \item{u}{
a "true" value of the regression function, may be provided to
   report risks at each iteration. This can be used to test the propagation condition with \code{u=0}
}
  \item{patchsize}{
positive integer defining the size of patches. Number of grid points within the patch is \code{(2*patchsize+1)^d} with \code{d} denoting the dimensionality of the design.
}
\item{maxni}{require growing sum of weights}
  \item{lambda}{explicit value of lambda}
  \item{data}{optional vector-valued images to be smoothed using the weighting scheme of the last step}
  \item{verbose}{logical: provide information on progress.}
}
\details{see \code{\link{vaws}}.
   Parameter \code{y}
   The procedure is supposed to produce superior results if the assumption of a
   local constant image is violated or if smooothness of discontinuities is desired.

   Function \code{vpawscov2} is intended for internal use in package \code{qMRI} only.
}
\value{ function \code{vpaws} returns
returns an object of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{numeric(0)}
                        \item{ni = "integer"}{integer(0)}
                        \item{mask = "logical"}{logical(0)}
                        \item{theta = "numeric"}{Estimates of regression function, \code{length: length(y)}}
                        \item{hseq = "numeric"}{sequence of bandwidths employed}
                        \item{mae = "numeric"}{Mean absolute error for each iteration step if u was specified,  numeric(0) else}
                        \item{psnr = "numeric"}{Peak signal-to-noise ratio for each iteration step if u was specified,  numeric(0) else}
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.Currently also uses factor \code{1/ni} instead of the correct
                        \code{sum(wij^2)/ni^2}}
                        \item{xmin = "numeric"}{numeric(0)}
                        \item{xmax = "numeric"}{numeric(0)}
                        \item{wghts = "numeric"}{numeric(0), ratio of distances \code{wghts[-1]/wghts[1]}}
                        \item{degree = "integer"}{0}
                        \item{hmax  = "numeric"}{effective hmax}
                        \item{sigma2  = "numeric"}{provided or estimated error variance}
                        \item{scorr = "numeric"}{scorr}
                        \item{family = "character"}{family}
                        \item{shape = "numeric"}{shape}
                        \item{lkern  = "integer"}{integer code for lkern,
                        1="Plateau", 2="Triangle", 3="Quadratic", 4="Cubic", 5="Gaussian"}
                        \item{lambda = "numeric"}{effective value of lambda}
                        \item{ladjust = "numeric"}{effective value of ladjust}
                        \item{aws = "logical"}{aws}
                        \item{memory = "logical"}{memory}
                        \item{homogen = "logical"}{homogen}
                        \item{earlystop = "logical"}{FALSE}
                        \item{varmodel = "character"}{"Constant"}
                        \item{vcoef = "numeric"}{numeric(0)}
                        \item{call = "function"}{the arguments of the call to \code{aws}}
  If \code{y} contained only information (condensed data) for positions within a mask, then the returned object only contains
  results for these positions.
}
\references{J. Polzehl, K. Tabelow (2019). Magnetic Resonance Brain Imaging:
             Modeling and Data Analysis Using R. Springer, Use R! series. Appendix A.
             Doi:10.1007/978-3-030-29184-6.

             J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06 .
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}
}
\note{
use \code{setCores='number of threads'} to enable parallel execution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also \code{\link{vaws}}, \code{\link{lpaws}}, \code{\link{vawscov}},\code{link{awsdata}}
}
\examples{
\dontrun{
setCores(2)
y <- array(rnorm(4*64^3),c(4,64,64,64))
yhat <- vpaws(y,kstar=20)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
