\name{awsweights}
\alias{awsweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate weight scheme that would be used in an additional aws step
}
\description{
Utility function to create a weighting scheme for an additional aws step.
Inteded to be used for illustrations only.
}
\usage{
awsweights(awsobj, spmin = 0.25, inx = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{awsobj}{
    object obtained by a call to function \code{aws}
}
  \item{spmin}{
    Size of the plateau in the adaptation kernel. 
}
  \item{inx}{
    either a matrix of dimension length(awsobj@dy) x number of points
    containing the integer coordinates of points of interest
    or 
    NULL. In the latter case the weight scheme for all points is
    generated.
}
}
\value{
   an array of either dimension
   awsobj@dy x number of points or awsobj@dy x awsobj@dy
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Joerg Polzehl, Vladimir Spokoiny,  Adaptive Weights Smoothing with applications to image restoration,  J. R. Stat. Soc. Ser. B Stat. Methodol.  62 ,  (2000)  , pp. 335--354

Joerg Polzehl, Vladimir Spokoiny,  Propagation-separation approach for local likelihood estimation,  Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362.

Joerg Polzehl, Kostas Papafitsoros, Karsten Tabelow (2020). Patch-Wise Adaptive Weights Smoothing in R, Journal of
             Statistical Software, 95(6), 1-27. doi:10.18637/jss.v095.i06.


%% ~put references to the literature/web site here ~
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{https://www.wias-berlin.de/people/polzehl/}
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also  \code{\link{aws}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% use one of  RShowDoc("KEYWORDS")
