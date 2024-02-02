\name{paws}
\alias{paws}
\alias{pawsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Adaptive weigths smoothing using patches
}
\description{The function implements a version the propagation separation approach that
uses patches instead of individuel voxels for comparisons in parameter space. Functionality is analog to function \code{\link{aws}}. Using patches allows for an improved
handling of locally smooth functions and in 2D and 3D for improved smoothness of
discontinuities at the expense of increased computing time.
}
\usage{
paws(y, hmax = NULL, mask=NULL, onestep = FALSE, aws = TRUE, family = "Gaussian",
     lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, shape = NULL,
     scorr = 0, spmin = 0.25, ladjust = 1, wghts = NULL, u = NULL,
     graph = FALSE, demo = FALSE, patchsize = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{array \code{y} containing the observe response (image intensity) data.
  \code{dim(y)} determines the dimensionality and extend of the grid design.
}
\item{mask}{logical array defining a mask. All computations are restricted to the mask.
}
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5}
  for 1D, 2D, 3D images, respectively.
  In case of \code{lkern="Gaussian"} the bandwidth is assumed to be given in full width half maximum (FWHM) units, i.e., \code{0.42466} times gridsize.
}
  \item{onestep}{
apply the last step only (use for test purposes only)
}
  \item{aws}{
logical: if TRUE structural adaptation (AWS) is used.
}
  \item{family}{\code{family} specifies the probability distribution. Default is \code{family="Gaussian"}, also implemented
   are "Bernoulli", "Poisson", "Exponential", "Volatility", "Variance" and "NCchi". \code{family="Volatility"} specifies a Gaussian distribution with
   expectation 0 and unknown variance. \code{family="Volatility"} specifies that \code{p*y/theta} is distributed as \eqn{\chi^2} with \code{p=shape}
   degrees of freedom. \code{family="NCchi"} uses a noncentral Chi distribution with \code{p=shape} degrees of freedom and noncentrality parameter \code{theta}
}
  \item{lkern}{
character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian".
The default "Triangle" is equivalent to using an Epanechnikov kernel, "Quadratic" and  "Cubic" refer to a Bi-weight and Tri-weight
 kernel, see Fan and Gijbels (1996). "Gaussian" is a truncated (compact support) Gaussian kernel.
 This is included for comparisons only and should be avoided due to its large computational costs.
}

  \item{aggkern}{
character: kernel used in stagewise aggregation, either "Triangle" or "Uniform"
}
  \item{sigma2}{
\code{sigma2} allows to specify the variance in case of \code{family="Gaussian"}. Not used if \code{family!="Gaussian"}.
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==length(y)} then \code{sigma2}
   is assumed to contain the pointwise variance of \code{y} and a heteroscedastic variance model is used.
}
  \item{shape}{Allows to specify an additional shape parameter for certain family models. Currently only used for family="Variance", that is \eqn{\chi}-Square distributed observations
  with \code{shape} degrees of freedom.
}
  \item{scorr}{
The vector \code{scorr} allows to specify a first order correlations of the noise for each coordinate direction,
    defaults to 0 (no correlation).
}
  \item{spmin}{
Determines the form (size of the plateau) in the adaptation kernel.
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
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}.
}
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}.
}
  \item{patchsize}{
  positive integer defining the size of patches. Number of grid points within the patch is \code{(2*patchsize+1)^d} with \code{d} denoting the dimensionality of the design.
}
}
\details{ see \code{\link{aws}. The procedure is supposed to produce superior results if the assumption of a local constant image is violated or if smooothness of discontinuities is desired.}
}
\value{
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
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.}
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

\seealso{See also \code{\link{aws}}, \code{\link{lpaws}}, \code{\link{vpaws}},\code{link{awsdata}}
}
\examples{\dontrun{
setCores(2)
y <- array(rnorm(64^3),c(64,64,64))
yhat <- paws(y,hmax=6)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
