CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cplxawss(y,mask,nv,n1,n2,n3,hakt,lambda,theta,s2,bi,
     1                thnew,s2new,ncores,lwght,wght,swjy)
C
C   y        observed values of regression function
C   nv       number of vector components
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   kritical value
C   theta    estimates from last step   (input)
C   si2      assumes same variance in components and independence
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   ncores   number of cores
C
C   wght     scaling factor for second and third dimension (larger values shrink)
C
!$      use omp_lib
      implicit none

      integer nv,n1,n2,n3,ncores
      logical aws
      integer mask(*)
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),s2(*),
     1  lambda,wght(2),hakt,lwght(*),swjy(nv,ncores),s2new(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision bii,biinv,sij,swj,swj2,z1,z2,z3,wj,hakt2,
     1       w1,w2,spmb,spf
      external lkern,KLdist2
      double precision lkern,KLdist2
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=4.d0/3.d0
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
C$OMP& SHARED(thnew,bi,nv,n1,n2,n3,hakt2,theta,
C$OMP& ih3,lwght,wght,y,swjy,mask,s2new,s2)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,thrednr,swj2)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(mask(iind).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C    nothing to do, final estimate is already fixed by control
         bii=bi(iind)/lambda
         biinv=1.d0/bii
         spmb=0.25d0/bii
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
                  if(mask(jind).eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=KLdist2(theta(1,iind),theta(1,jind),s2(iind))
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-0.25d0))
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bi(iind)=swj
         z1=y(1,iind)-thnew(1,iind)
         z2=y(2,iind)-thnew(2,iind)
         if(swj.gt.1.1d0) THEN
            s2new(iind) = (z1*z1+z2*z2)/2.d0*swj2/(swj2-1.d0)
         ELSE
            s2new(iind) = s2(iind)
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi)
      RETURN
      END

      double precision function KLdist2(thi,thj,s2)
      implicit none
      double precision thi(2), thj(2), s2, z1, z2
      z1 = thi(1)-thj(1)
      z2 = thi(2)-thj(2)
      KLdist2 = (z1*z1+z2*z2)/s2
      RETURN
      END
