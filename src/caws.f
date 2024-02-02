C
C    Copyright (C) 2004 Weierstrass-Institut fuer
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the aws package and contain
C  FORTRAN 77 code needed in R functions aws, vaws,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid
C
C   this is a reimplementation of the original aws procedure
C
C   should be slightly slower for non-Gaussian models (see function kldist)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Model=1    Gaussian
C          Model=2    Bernoulli
C          Model=3    Poisson
C          Model=4    Exponential
C          Model=5    Variance
C          Model=6    Noncentral Chi (Gaussian approximation,
C                     variance mean dependence is introduces via factor bii)
C
C     computing dlog(theta) and dlog(1.d0-theta) outside the AWS-loops
C     will reduces computational costs at the price of readability
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function kldist(model,thi,thj)
      implicit none
      integer model
      double precision thi,thj,z,tthi
      IF (model.eq.1) THEN
C        Gaussian
         z=thi-thj
         kldist=z*z
      ELSE IF (model.eq.2) THEN
C        Bernoulli
         kldist=0.d0
         tthi=(1.d0-thi)
         IF (thi.gt.1.d-10) kldist=kldist+thi*log(thi/thj)
         IF (tthi.gt.1.d-10) kldist=kldist+tthi*log(tthi/(1.d0-thj))
      ELSE IF (model.eq.3) THEN
C        Poisson
         kldist=0.d0
         IF (thi.gt.1.d-10) kldist=thi*log(thi/thj)-thi+thj
      ELSE IF (model.eq.4) THEN
C        Exponential
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE IF (model.eq.5) THEN
C        Variance
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE IF (model.eq.6) THEN
C        Noncentral Chi with Gaussian approximation
         z=thi-thj
         kldist=z*z
      ELSE
C        use Gaussian
         z=thi-thj
         kldist=z*z
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function lkern(kern,xsq)
      implicit none
      integer kern
      double precision xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=exp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws(y,pos,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1                bi0,ai,model,kern,spmin,lwght,wght)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,iindp,jindp
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1       hakt2,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,n1,n2,n3,hakt2,theta,
C$OMP& ih3,lwght,wght,y,pos)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,thetai,bii,swj,swj2,iindp,jindp,
C$OMP& swj0,swjy,sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         thetai=theta(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj
                  z1=jw1
                  z1=z2+z1*z1
                  IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(jindp))
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                     END IF
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jindp)
               END DO
            END DO
         END DO
         ai(iindp)=swjy
         bi(iindp)=swj
         bi2(iindp)=swj2
         bi0(iindp)=swj0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws6(y,pos,n1,n2,n3,hakt,lambda,theta,fnc,bi,
     1                 bi2,bi0,ai,kern,spmin,lwght,wght)
C
C  aws for nc-chi differs from caws through arguments model (missing) and fnc (add)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,kern,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),spmin,spf,fnc(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,iindp,jindp
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z,z1,z2,z3,wj,
     1       hakt2,w1,w2,fnci
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,n1,n2,n3,hakt2,theta,fnc,
C$OMP& ih3,lwght,wght,y,pos)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,thetai,bii,swj,swj2,
C$OMP& swj0,swjy,sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,fnci,z,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         thetai=theta(iindp)
         bii=bi(iindp)/lambda
         fnci=fnc(iindp)
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj
                  z1=jw1
                  z1=z2+z1*z1
                  IF (aws) THEN
                     z=thetai-theta(jindp)
                     sij=bii*z*z/(fnci+fnc(jindp))
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                     END IF
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jindp)
               END DO
            END DO
         END DO
         ai(iindp)=swjy
         bi(iindp)=swj
         bi2(iindp)=swj2
         bi0(iindp)=swj0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws(y,si2,pos,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1           bi0,ai,model,kern,spmin,lwght,wght)
C
C differs from caws by arguments si2 (inverse var)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none
      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),si2(*),spmin
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,iindp,jindp
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1        hakt2,sv1,sv2,spf,w1,w2,wjsi2
      w1=wght(1)
      w2=wght(2)
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,n1,n2,n3,hakt2,hakt,theta,
C$OMP& lwght,wght,y,pos)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,dlw12,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2)
C$OMP& PRIVATE(iind,thetai,bii,swj,iindp,jindp,
C$OMP& swj2,swj0,swjy,sij,sv1,sv2,i1,i2,i3,wj,
C$OMP& j3,jw3,jind3,z3,jwind3,
C$OMP& j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,wjsi2)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         thetai=theta(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jindp)
                  IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(jindp))
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  wjsi2=wj*si2(jindp)
                  swj=swj+wjsi2
                  swj2=swj2+wj*wjsi2
                  swjy=swjy+wjsi2*y(jindp)
               END DO
            END DO
         END DO
         ai(iindp)=swjy
         bi(iindp)=swj
         bi2(iindp)=swj2
         bi0(iindp)=swj0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded) with variance - mean model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgaws(y,pos,si2,n1,n2,n3,hakt,lambda,
     1        theta,bi,bi2,bi0,gi,gi2,ai,kern,spmin,lwght,wght)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none
      external lkern
      double precision lkern
      integer n1,n2,n3,kern,pos(*)
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),si2(*),gi2(*),spmin,gi(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,iindp,jindp
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1        hakt2,sv1,sv2,spf,z,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      w1=wght(1)
      w2=wght(2)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,n1,n2,n3,hakt2,theta,
C$OMP& lwght,wght,y,pos,gi2,gi)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,dlw12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2)
C$OMP& PRIVATE(iind,thetai,bii,swj,iindp,jindp,
C$OMP& swj2,swj0,swjy,sij,sv1,sv2,i1,i2,i3,wj,
C$OMP& j3,jw3,jind3,z3,jwind3,
C$OMP& j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         thetai=theta(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jindp)
                  z1=-jw1
                  z1=z2+z1*z1
                  IF (aws) THEN
C
C      gaussian case only
C
                     z=(thetai-theta(jindp))
                     sij=bii*z*z
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                     END IF
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  swj=swj+wj*si2(jindp)
                  swj2=swj2+wj*wj*si2(jindp)
                  swjy=swjy+wj*si2(jindp)*y(jindp)
               END DO
            END DO
         END DO
         ai(iindp)=swjy
         bi(iindp)=swj
         bi2(iindp)=swj2
         bi0(iindp)=swj0
         gi(iindp)=sv1
         gi2(iindp)=sv2
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,gi,gi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vaws(y,pos,nv,n1,n2,n3,hakt,lambda,theta,bi,vred,
     1                thnew,ncores,spmin,lwght,wght,swjy)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      use omp_lib
      implicit none

      integer nv,n1,n2,n3,ncores,pos(*)
      logical aws
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),vred(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr,iindp,jindp
      double precision bii,biinv,sij,swj,z,z1,z2,z3,wj,hakt2,
     1        w1,w2,spmb,swj2
      external lkern
      double precision lkern
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,vred,nv,n1,n2,n3,hakt2,theta,
C$OMP& ih3,lwght,wght,y,swjy,pos)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,swj2,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         bii=bi(iindp)/lambda
         biinv=1.d0/bii
         spmb=spmin/bii
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=0.d0
                     DO k=1,nv
                        z=theta(k,iindp)-theta(k,jindp)
                        sij=sij+z*z
                        IF(sij.ge.biinv) CYCLE
                     END DO
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-spmin))
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jindp)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iindp)=swjy(k,thrednr)/swj
         END DO
         bi(iindp)=swj
         vred(iindp)=swj2/swj/swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi)
      RETURN
      END
      subroutine vaws2cov(y,pos,nv,nvd,n1,n2,n3,hakt,lambda,theta,bi,
     1                vred,thnew,invcov,ncores,spmin,lwght,wght,swjy,
     2                thi,invcovi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      use omp_lib
      implicit none

      integer nv,n1,n2,n3,ncores,nvd,pos(*)
      logical aws
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),invcov(nvd,*),
     2  vred(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr,iindp,jindp
      double precision bii,sij,swj,z1,z2,z3,wj,hakt2,w1,w2,swj2
      double precision thi(nv,*),invcovi(nvd,*)
      external lkern, KLdistsi
      double precision lkern, KLdistsi
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,nvd,n1,n2,n3,hakt2,theta,invcov,
C$OMP& ih3,lwght,wght,y,swjy,pos,thi,
C$OMP& invcovi,vred)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,swj,swj2,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,thrednr,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         DO k=1,nv
            thi(k,thrednr)=theta(k,iindp)
         END DO
         DO k=1,nvd
            invcovi(k,thrednr)=invcov(k,iindp)
         END DO
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj2=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=bii*KLdistsi(thi(1,thrednr),theta(1,jindp),
     1                   invcovi(1,thrednr),nv)
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jindp)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iindp)=swjy(k,thrednr)/swj
         END DO
         bi(iindp)=swj
         vred(iindp)=swj2/swj/swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi,vred)
      RETURN
      END
