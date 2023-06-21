\name{TV_denoising}
\alias{TV_denoising}
\alias{TGV_denoising}
\alias{TV_denoising_colour}
\alias{TGV_denoising_colour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   TV/TGV denoising of image data
}
\description{
   Total variation and total generalized variation are classical energy minimizing methods
   for image denoising.
}
\usage{
TV_denoising(datanoisy, alpha, iter = 1000, tolmean = 1e-06,
             tolsup = 1e-04, scale = 1, verbose=FALSE)
TGV_denoising(datanoisy, alpha, beta, iter = 1000, tolmean = 1e-06,
              tolsup = 1e-04, scale = 1, verbose=FALSE)
TV_denoising_colour(datanoisy, alpha, iter = 1000, tolmean = 1e-06,
                    tolsup = 1e-04, scale = 1, verbose=FALSE)
TGV_denoising_colour(datanoisy, alpha, beta, iter = 1000, tolmean = 1e-06,
                     tolsup = 1e-04, scale = 1, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datanoisy}{
    matrix of noisy 2D image data. In case of \code{TV_denoising_colour}
    and \code{TGV_denoising_colour} and array with third dimension refering to
    RGB channels.
}
  \item{alpha}{
    TV regularization parameter.
}
  \item{beta}{
    additional TGV regularization parameter.
}
  \item{iter}{
    max. number of iterations
}
  \item{tolmean}{
    requested accuracy for mean image correction
}
  \item{tolsup}{
    requested accuracy for max (over pixel) image correction
}
  \item{scale}{
    image scale
}
\item{verbose}{
  report convergence diagnostics.
}
}
\details{
   Reimplementation of original matlab code by Kostas Papafitsoros (WIAS).
}
\value{
   TV/TGV reconstructed image data (2D array)
}
\references{
             J. Polzehl, K. Papafitsoros, K. Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06.

   Rudin, L.I., Osher, S. and Fatemi, E. (1992). Nonlinear total variation based
   noise removal algorithms. Phys. D, 60, 259-268. DOI: 10.1016/0167-2789(92)90242-F.

   Bredies, K., Kunisch, K. and Pock, T. (2010). Total Generalized Variation.
   SIAM J. Imaging Sci., 3, 492-526. DOI:10.1137/090769521.
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% use one of  RShowDoc("KEYWORDS")
