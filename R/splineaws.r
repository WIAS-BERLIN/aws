#
#    R - functions  for  Adaptive Weights Smoothing (AWS)
#    in (generalized) local polynomial regression models in 1D and 2D         
#
#    Copyright (C) 2007 Weierstrass-Institut fuer                          
#                       Angewandte Analysis und Stochastik (WIAS)         
#
#    Author:  Joerg Polzehl                                                
#
#  This program is free software; you can redistribute it and/or modify   
#  it under the terms of the GNU General Public License as published by   
#  the Free Software Foundation; either version 2 of the License, or      
#  (at your option) any later version.                                    
#
#  This program is distributed in the hope that it will be useful,        
#  but WITHOUT ANY WARRANTY; without even the implied warranty of         
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
#  GNU General Public License for more details.                           
#
#  You should have received a copy of the GNU General Public License      
#  along with this program; if not, write to the Free Software            
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,  
#  USA.
#
##############################################################################
#
#   Hierarchical local polynomial aws
#
#   order up to cubic     
#
##############################################################################
hiraraws <- function(y,degree=3,hmax=NULL,qlambda=NULL,lkern="Triangle",
                     graph=FALSE,spmin=0.25)
{ 
#
#   first fill numerical derivatives
#
  n <- length(y)
  n1 <- n%/%3
  n2 <- n%/%3
  n3 <- n%/%5
  yd1 <- y[3*(1:n1)]-y[3*(1:n1)-2]
  yd2 <- y[3*(1:n2)]-2*y[3*(1:n2)-1]+y[3*(1:n2)-2]
  yd3 <- y[5*(1:n3)]-2*y[4*(1:n3)-1]+2*y[4*(1:n3)-3]-y[4*(1:n3)-4]
#
# functional values
#
  ergsD3 <- aws1(yd3,NULL,NULL,NULL,0,hmax=hmax,glambda=qlambda,lkern=lkern,spmin=spmin,graph=graph)
  fw3 <- ergsD3$fw
  ni3 <- ergsD3$ni
  sigma2D3 <- ergsD3$sigma2
  fw3all <- as.vector(outer(c(1,.8,.6,.4,.2),fw3[-n3]),"*")+outer(1-c(1,.8,.6,.4,.2),fw3[-1],"*"))
  ni3all <- as.vector(outer(c(1,.8,.6,.4,.2),ni3[-n3],"*")+outer(1-c(1,.8,.6,.4,.2),ni3[-1],"*"))
  fw3all <- c(rep(fw3[1],2),fw3)
  ni3all <- c(rep(ni3[1],2),ni3)
  nall <- length(fw3all)
  if(nall<n) {
     fw3all <- c(fw3all,rep(fw3[n3],n-nall))
     ni3all <- c(ni3all,rep(ni3[n3],n-nall))
  }
  fw32 <- fw3all[3*(1:n2)-1]
  ni32 <- ni3all[3*(1:n2)-1]
  ergsD2 <- aws1(yd2,fw32,ni32,sigma2D3,1,hmax=hmax,glambda=qlambda,lkern=lkern,spmin=spmin,graph=graph)
  fw2 <- ergsD2$fw
  ni2 <- ergsD2$ni
  sigma2D2 <- ergsD2$sigma2
  fw2all <- as.vector(outer((3:1)/3,fw2[-n2]),"*")+outer((0:2)/3,fw2[-1],"*"))
  ni2all <- as.vector(outer((3:1)/3,ni2[-n2],"*")+outer((0:2)/3,ni2[-1],"*"))
  fw2all <- c(fw2[1],fw2all)
  ni2all <- c(ni2[1],ni2all)
  nall <- length(fw2all)
  if(nall<n) {
     fw2all <- c(fw2all,rep(fw2[n2],n-nall))
     ni2all <- c(ni2all,rep(ni2[n2],n-nall))
  }
  fw31 <- rbind(fw2all[2*(1:n1)-1],fw3all[2*(1:n1)-1])
  ni31 <- rbind(ni2all[2*(1:n1)-1],ni3all[2*(1:n1)-1])
  ergsD1 <- aws1(yd1,fw31,ni31,c(sigma2D2,sigma2D3),2,hmax=hmax,glambda=qlambda,lkern=lkern,spmin=spmin,graph=graph)
  fw1 <- ergsD1$fw
  ni1 <- ergsD1$ni
  sigma2D1 <- ergsD1$sigma2
  fw1all <- as.vector(outer((3:1)/3,fw1[-n1]),"*")+outer((0:2)/3,fw1[-1],"*"))
  ni1all <- as.vector(outer((3:1)/3,ni1[-n1],"*")+outer((0:2)/3,ni1[-1],"*"))
  fw1all <- c(fw1[1],fw1all)
  ni1all <- c(ni1[1],ni1all)
  nall <- length(fw1all)
  if(nall<n) {
     fw1all <- c(fw1all,rep(fw1[n1],n-nall))
     ni1all <- c(ni1all,rep(ni1[n1],n-nall))
  }
  fw30 <- rbind(fw1all,fw2all,fw3all)
  ni30 <- rbind(ni1all,ni2all,ni3all)
  ergsD0 <- aws1(y,fw31,ni31,c(sigma2D1,sigma2D2,sigma2D3),3,hmax=hmax,glambda=qlambda,lkern=lkern,spmin=spmin,graph=graph)
  z <- list(y=y,yhat=ergsD0$fw,yhat1D=fw1all,yhat2D=fw2all,yhat3D=fw3all,
                ni=ergsD0$ni,ni1D=ni1all,ni2D=ni2all,ni3D=ni3all)
  class(z) <- "awshir.gaussian"
  z
}

aws1 <- function(y,fw,ni,sigma2D,degree,hmax,qlambda=.965,lkern="Triangle",spmin=.25,graph=FALSE){
   n <- length(y)
   if(is.null(fw)||is.null(ni)){
      fw <- rep(0,n)
      ni <- rep(1,n)
      degree <- 0
      sigma2D <- 1
   }
   lambda <- qchisq(qlambda,1)
   lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Uniform=1,
	            Gaussian=5,2)
   sigma2 <- IQRdiff(y)^2
   hakt <- 1.25
   hincr <- 1.25
   zobj <- list(yhat=y,bi=rep(1,n))
   while(hakt<=hmax){
   zobj <- .Fortran("hiraws",as.double(y),
                             as.integer(n),
                             as.double(sigma2),
                             as.integer(degree),
                             as.double(fw),
                             as.double(ni),
                             as.double(sigma2D),
                             as.double(hakt),
                             as.double(lambda),
                             as.integer(lkern),
                             as.double(spmin),
                             as.double(zobj$yhat)
                             as.double(zobj$bi),
                             yhat=double(n),
                             bi=double(n),
                             DUP=FALSE,
                             PACKAGE="aws")[c("bi","yhat")]
  hakt <- hakt*hincr
}
list(fw=zobj$yhat,ni=zobj$bi,sigma2=sigma2)
}

