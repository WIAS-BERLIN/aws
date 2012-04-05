########## extract()

extract <- function(x, ...) cat("Data extraction not defined for this class:",class(x),"\n")

setGeneric("extract", function(x, ...) standardGeneric("extract"))

setMethod("extract","aws",function(x, what="y"){
  what <- tolower(what) 
  z <- list(NULL)
  if("y" %in% what) z$y <- x@y
  if("yhat" %in% what) z$yhat <- if(x@degree==0) drop(x@theta) else x@theta
  if("x" %in% what) z$x <- x@x
  if("sigma2" %in% what) z$sigma2 <- x@sigma2
  if("ni" %in% what) z$ni <- x@ni
  if("mask" %in% what) z$mask <- x@mask
  invisible(z)
})
setMethod("extract","awssegment",function(x, what="y"){
  what <- tolower(what) 
  z <- list(NULL)
  if("y" %in% what) z$y <- x@y
  if("yhat" %in% what) z$yhat <- drop(x@theta) 
  if("segment" %in% what) z$segment <- x@segment
  if("x" %in% what) z$x <- x@x
  if("sigma2" %in% what) z$sigma2 <- x@sigma2
  if("ni" %in% what) z$ni <- x@ni
  if("mask" %in% what) z$mask <- x@mask
  invisible(z)
})
setMethod("extract","kernsm",function(x, what="y"){
  what <- tolower(what) 
  z <- list(NULL)
  if("y" %in% what) z$y <- x@y
  if("yhat" %in% what) z$yhat <- x@yhat
  if("vred" %in% what) z$vred <- x@vred
  if("vhat" %in% what) z$vhat <- (median(abs(x@y[-1]-x@y[-n]))/.9538)^2/x@vred
  invisible(z)
})
setMethod("extract","ICIsmooth",function(x, what="y"){
  what <- tolower(what) 
  z <- list(NULL)
  if("y" %in% what) z$y <- x@y
  if("yhat" %in% what) z$yhat <- x@yhat
  if("vhat" %in% what) z$vhat <- x@vhat
  if("vred" %in% what) z$vred <- x@sigma^2/x@vhat
  if("hbest" %in% what) z$hbest <- x@hbest
  invisible(z)
})
################################################################
#                                                              #
# Section for summary(), print(), show() functions (generic)   #
#                                                              #
################################################################

setMethod("print", "aws",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls   :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@dy, collapse="x"), "\n")
    cat("  Max. Bandwidth       :", x@hmax, "  degree=",x@degree,"\n")
    cat("  Lambda               :", x@lambda, "  (ladj=",x@ladj,")\n")
    cat("  Slots", slotNames(x), "\n")
    invisible(NULL)
})
setMethod("summary", "aws",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls   :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@dy, collapse="x"), "\n")
    cat("  Max. Bandwidth       :", object@hmax, "  degree=",object@degree,"\n")
    cat("  Lambda               :", object@lambda, "  (ladj=",object@ladj,")\n")
    cat("  mean sum of weights  :", mean(object$ni), "\n")
    cat("\n")
    invisible(NULL)
})
setMethod("show", "aws",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls   :\n")
    print(object@call)
    d <- length(object@dy)
    if(d==1) if(prod(par()$mfrow<2)) par(mfrow=c(1,2), ...)
    invisible(NULL)
})