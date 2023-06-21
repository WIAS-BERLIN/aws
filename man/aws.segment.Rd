\name{aws.segment}
\alias{aws.segment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Segmentation by adaptive weights for Gaussian models. }
\description{
The function implements a modification of the adaptive weights smoothing algorithm for
segmentation into three classes. The
}
\usage{
aws.segment(y, level, delta = 0, hmax = NULL, hpre = NULL, mask =NULL,
            varmodel = "Constant", lkern = "Triangle", scorr = 0, ladjust = 1,
            wghts = NULL, u = NULL, varprop = 0.1, ext = 0, graph = FALSE,
            demo = FALSE, fov=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed response data. \code{dim(y)} determines the dimensionality and extend of the grid design.}
  \item{level}{ center of second class }
  \item{delta}{ half width of second class }
  \item{hmax}{ \code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax=250, 12, 5} for \code{dd=1, 2, 3}, respectively.}
  \item{hpre}{ Describe \code{hpre} Bandwidth used for an initial nonadaptive estimate. The first estimate
   of variance parameters  is obtained from residuals with respect to this estimate. }
  \item{mask}{ optional logical mask, same dimensionality as \code{y}}
  \item{varmodel}{Implemented are "Constant", "Linear" and "Quadratic" refering to a polynomial model of degree 0 to 2.}
  \item{lkern}{ character: location kernel, either "Triangle", "Plateau", "Quadratic", "Cubic" or "Gaussian".
  The default "Triangle" is equivalent to using an Epanechnikov kernel, "Quadratic" and  "Cubic" refer to a Bi-weight and Tri-weight
   kernel, see Fan and Gijbels (1996). "Gaussian" is a truncated (compact support) Gaussian kernel.
   This is included for comparisons only and should be avoided due to its large computational costs.}
  \item{scorr}{The vector \code{scorr} allows to specify a first order correlations of the noise for each coordinate direction,
    defaults to 0 (no correlation). }
  \item{ladjust}{ factor to increase the default value of lambda }
  \item{wghts}{ \code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
  \item{u}{ a "true" value of the regression function, may be provided to
   report risks at each iteration. This can be used to test the propagation condition with \code{u=0} }
  \item{varprop}{ Small variance estimates are replaced by \code{varprop} times the mean variance.}
  \item{ext}{ Intermediate results are fixed if the test statistics exceeds the critical value by \code{ext}.}
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}. }
  \item{demo}{ If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{fov}{Field of view. Size of region (sample size) to adjust for in multiscale testing.}
}
\details{
  The image is segmented into three parts by performing multiscale tests
of the hypotheses H1
\code{value >= level - delta} and H2 \code{value <= level + delta}.
Pixel where the first hypotesis is rejected are classified as \code{-1} (segment 1)
while rejection of H2 results in classification  \code{1} (segment 3).
Pixel where neither H1 or H2 are rejected ar assigned to a value \code{0} (segment 2). Critical values for the tests are adjusted for smoothness at the different scales inspected in the iteration process using results from multiscale testing,
see e.g. Duembgen and Spokoiny (2001). Critical values also depend on the
size of the region of interest specified in parameter \code{fov}.

Within segment 2 structural adaptive smoothing is performed while if a pair of pixel belongs to segment 1 or segment 3 the corresponding weight will be nonadaptive.
}
\value{
returns anobject of class \code{aws} with slots
                        \item{y = "numeric"}{y}
                        \item{dy = "numeric"}{dim(y)}
                        \item{x = "numeric"}{numeric(0)}
                        \item{ni = "integer"}{integer(0)}
                        \item{mask = "logical"}{logical(0)}
                        \item{segment = "integer"}{Segmentation results, class numbers 1-3 }
                        \item{theta = "numeric"}{Estimates of regression function, \code{length: length(y)}}
                        \item{mae = "numeric"}{Mean absolute error for each iteration step if u was specified,  numeric(0) else}
                        \item{var = "numeric"}{approx. variance of the estimates of the regression function. Please note that this does not reflect variability due to randomness of weights.}
                        \item{xmin = "numeric"}{numeric(0)}
                        \item{xmax = "numeric"}{numeric(0)}
                        \item{wghts = "numeric"}{numeric(0), ratio of distances \code{wghts[-1]/wghts[1]}}
                        \item{degree = "integer"}{0}
                        \item{hmax  = "numeric"}{effective hmax}
                        \item{sigma2  = "numeric"}{provided or estimated error variance}
                        \item{scorr = "numeric"}{scorr}
                        \item{family = "character"}{"Gaussian"}
                        \item{shape = "numeric"}{NULL}
                        \item{lkern  = "integer"}{integer code for lkern,
                        1="Plateau", 2="Triangle", 3="Quadratic", 4="Cubic", 5="Gaussian"}
                        \item{lambda = "numeric"}{effective value of lambda}
                        \item{ladjust = "numeric"}{effective value of ladjust}
                        \item{aws = "logical"}{aws}
                        \item{memory = "logical"}{memory}
                        \item{homogen = "logical"}{FALSE}
                        \item{earlystop = "logical"}{FALSE}
                        \item{varmodel = "character"}{varmodel}
                        \item{vcoef = "numeric"}{estimated parameters of the variance model}
                        \item{call = "function"}{the arguments of the call to \code{aws.gaussian}}
}
\note{ This function is still experimental and may be changes considerably in future.}
\references{J. Polzehl, H.U. Voss, K. Tabelow (2010). Structural adaptive segmentation for
   statistical parametric mapping, NeuroImage, 52, pp. 515--523. DOI:10.1016/j.neuroimage.2010.04.241

  Duembgen, L. and Spokoiny, V. (2001). Multiscale testing of qualitative hypoteses. Ann. Stat. 29, 124--152.

  Polzehl, J. and Spokoiny, V. (2006). Propagation-Separation Approach for Local
  Likelihood Estimation. Probability Theory and Related Fields. 3 (135)
  335 - 362. DOI:10.1007/s00440-005-0464-1
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}  }
\seealso{\code{\link{aws}}, \code{\link{aws.gaussian}} }
\examples{
require(aws)
}
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ regression }
