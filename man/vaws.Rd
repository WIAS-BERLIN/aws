\name{vaws}
\alias{vaws}
\alias{vawscov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
vector valued version of function \code{\link{aws}}
The function implements the propagation separation approach to
nonparametric smoothing (formerly introduced as Adaptive weights smoothing)
for varying coefficient likelihood models with vector valued response on a 1D, 2D or 3D grid.
}
\description{The function implements a version the propagation separation approach that
uses vector valued instead of scalar responses.
}
\usage{
vaws(y, kstar = 16, sigma2 = 1, mask = NULL, scorr = 0, spmin = 0.25,
     ladjust = 1, wghts = NULL, u = NULL, maxni = FALSE)
vawscov(y, kstar = 16, invcov = NULL, mask = NULL, scorr = 0, spmin = 0.25,
          ladjust = 1, wghts = NULL, u = NULL, maxni = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{\code{y} contains the observed response data. \code{dim(y)} determines the dimensionality and extend of the grid design. First component varies over components of the response vector.
}
  \item{kstar}{
maximal number of steps to employ. Determines maximal bandwidth.
}
  \item{sigma2}{
specifies a homogeneous error variance.
}
\item{invcov}{
 array of voxelwise inverse covariance matrixes, first index corresponds to
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
  \item{maxni}{
If TRUE use \eqn{max_{l<=k}(N_i^{(l)}} instead of \eqn{(N_i^{(k)}} in the definition of the statistical penalty.
}
}
\details{ see \code{\link{aws}}. Expets vector valued responses. Currently only implements the case of additive Gaussian errors.

}
\value{
returns anobject of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{numeric(0)}
                        \item{ni = "integer"}{integer(0)}
                        \item{mask = "logical"}{logical(0)}
                        \item{theta = "numeric"}{Estimates of regression function, \code{length: length(y)}}
                        \item{hseq = "numeric"}{sequence of bandwidths employed}
                        \item{mae = "numeric"}{Mean absolute error for each iteration step if u was specified,  numeric(0) else}
                        \item{psnr = "numeric"}{Peak signal-to-noise ratio for each iteration step if u was specified,  numeric(0) else}
                        \item{var = "numeric"}{approx. variance of the estimates
                            of the regression function.
                            Please note that this does not reflect variability due to
                            randomness of weights.}
                        \item{xmin = "numeric"}{numeric(0)}
                        \item{xmax = "numeric"}{numeric(0)}
                        \item{wghts = "numeric"}{numeric(0), ratio of distances \code{wghts[-1]/wghts[1]}}
                        \item{degree = "integer"}{0}
                        \item{hmax  = "numeric"}{effective hmax}
                        \item{sigma2  = "numeric"}{provided or estimated (inverse) error variance}
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

J. Polzehl, V. Spokoiny,  Adaptive Weights Smoothing with applications to image restoration,
J. R. Stat. Soc. Ser. B Stat. Methodol.  62 ,  (2000)  , pp. 335--354. DOI:10.1111/1467-9868.00235.

J. Polzehl, V. Spokoiny,  Propagation-separation approach for local likelihood estimation,
Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362. DOI:10.1007/s00440-005-0464-1.
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}
}
\note{
use \code{setCores='number of threads'} to enable parallel execution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also \code{\link{aws}}, \code{\link{vpaws}},\code{link{awsdata}}
}
\examples{\dontrun{
setCores(2)
y <- array(rnorm(4*64^3),c(4,64,64,64))
yhat <- vaws(y,kstar=20)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
