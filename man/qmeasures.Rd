\name{qmeasures}
\alias{qmeasures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Quality assessment for image reconstructions.
}
\description{
   Computes selected criteria for quality assessments of
}
\usage{
qmeasures(img, ref,
  which = c("PSNR", "MAE", "MSE", "RMSE", "SSIM", "MAGE", "RMSGE"),
  mask = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{
   2D/3D image, object of class \code{"aws", "ICIsmooth", "kernsm", "nlmeans"}
   or array.
}
  \item{ref}{
   Reference image (array, matrix or vector) for comparison.
}
  \item{which}{
   Criterion to use for Quality assessment. Please specify a subset of
   \code{"PSNR"} (Peak Signal to Noise Ratio),
   \code{"MAE"} (Mean Absolute Error),
   \code{"MSE"} (Mean Squared Error),
   \code{"RMSE"} (Root Mean Squared Error),
   \code{"SSIM"} (Structural SIMilarity),
   \code{"MAGE"} (Mean Absolute Gradient Error),
   \code{"RMSGE"} (Root Mean Squared Gradient Error).
}
  \item{mask}{
   Logical of same dimension as img/ref. Calculation can be restricted to mask.
}
}
\details{
 Calculates specified quality indices.
}
\value{
 A vector with names as specified in \code{which}.
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{https://www.wias-berlin.de/people/polzehl/}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
