CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded) !D
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsp1b(y,fix,nfix,n,degr,hw,hakt,hhom,lambda,theta,
     1        bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   kern     specifies the location kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C
C   temporary arrays set for maximum degree 2
C
      use omp_lib
      implicit none
      external kldistp,lkern
      double precision kldistp,lkern
      integer n,kern,degr,ind(*),nfix,fix(*)
      logical aws,lfix
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,spmin,
     1       bi2(*),hakt,lw(*),w(*),hw,sw(*),slw(*),hhom(n,2)
      integer ih,j1,k,iind,jind,dlw,clw,jw1,
     2        dp1,dp2,ihs,csw,dsw,l,thrednr,trl,trs
      double precision bii(5),sij,swj(5),swj2(5),swj0(5),swjy(5),z1,wj,
     1       hakt2,thij(3),thi(3),zz(5),lwj,yj,hs2,hs,z,cc,spf,
     2       hhommax,hhommin,az1,hfixmax,hnfix,ssij,spmax,
     3       hhomimin,hhomimax
C!$    integer omp_get_thread_num
C!$    external omp_get_thread_num
C   arrays with variable length are organized as
C   theta(n,dp1)
C   bi(n,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      lfix=nfix.gt.(degr+1)
      hnfix=nfix
      hnfix=max(hnfix,.2d0*hakt)
      spmax=degr+1
      spf=1.d0/(spmax-spmin)
      if(degr.eq.0) THEN
         dp1=1
         dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=2
         dp2=3
      ELSE
         dp1=3
         dp2=5
      END IF
      hakt2=hakt*hakt
      ih=FLOOR(hakt)
      dlw=2*ih+1
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=FLOOR(hs)
      dsw=2*ihs+1
      csw=ihs+1
C   compute location weights first
      DO j1=1,dlw
         z1=clw-j1
         lw(j1)=lkern(kern,z1*z1/hakt2)
      END DO
      cc=0.0d0
      call smwghts1(lw,hakt,hw,slw,dlw,dsw,cc)
      thrednr = 0
C  now stochastic term

C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(y,fix,nfix,n,degr,hw,hakt,hhom,lambda,theta,
C$OMP&        bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C$OMP& FIRSTPRIVATE(aws,lfix,hnfix,ih,spf,spmax,dp1,dp2,hakt2,dlw,clw,
C$OMP&        hs2,hs,ihs,csw,dsw)
C$OMP& PRIVATE(j1,k,iind,jind,jw1,l,bii,sij,swj,swj2,swj0,swjy,z1,
C$OMP&         wj,thij,thi,zz,lwj,yj,z,cc,hhommax,hhommin,az1,hfixmax,
C$OMP&         ssij,hhomimin,hhomimax,thrednr,trl,trs)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n

!$       thrednr = omp_get_thread_num()
         trl = thrednr*dlw
         trs = thrednr*dsw
         zz(1)=1.d0
         IF (fix(iind).ne.0) CYCLE
C    nothing to do, final estimate is already fixed by control
         hhomimin=hhom(iind,1)
         hhomimax=hhom(iind,2)
         hhommax=hakt
         hhommin=hakt
         hfixmax=max(hhomimin,hhomimax)
         DO k=1,dp1
            thi(k)=theta(iind+(k-1)*n)
         END DO
         DO k=1,dp2
            bii(k)=bi(iind+(k-1)*n)/lambda
         END DO
C   scaling of sij outside the loop
         DO jw1=1,dlw
            w(jw1+trl)=0.d0
            jind=jw1-clw+iind
            if(jind.lt.1.or.jind.gt.n) CYCLE
            wj=lw(jw1)
            z1=jw1-clw
            zz(2)=z1
            zz(3)=z1*z1
            az1=abs(z1)
            IF (aws.and.(z1.le.-hhomimin.or.z1.ge.hhomimax)) THEN
               DO k=1,dp1
                  thij(k)=theta(jind+(k-1)*n)
               END DO
               thij(1)=thij(1)-thij(2)*z1
               IF (dp1.gt.2) THEN
                  thij(1)=thij(1)+thij(3)*zz(3)
                  thij(2)=thij(2)-2.d0*thij(3)*z1
               END IF
C
C           get difference of thetas
C
               DO k=1,dp1
                  thij(k)=thi(k)-thij(k)
               END DO
               sij=kldistp(dp1,thij,bii,ind)
               IF (sij.le.spmax) THEN
                  hfixmax=max(hfixmax,az1)
                  IF (sij.gt.spmin) THEN
                      ssij=1.d0-spf*(sij-spmin)
                      w(jw1+trl)=wj*ssij
                      if(z1.gt.0) THEN
                         hhommax=min(hhommax,az1)
                      ELSE
                         hhommin=min(hhommin,az1)
                      END IF
                  ELSE
                     w(jw1+trl)=wj
                  END IF
               ELSE
                  w(jw1+trl)=0.d0
                  if(z1.gt.0) THEN
                     hhommax=min(hhommax,az1)
                  ELSE
                     hhommin=min(hhommin,az1)
                  END IF
               END IF
            ELSE
               w(jw1+trl)=wj
            END IF
         END DO
C
C      Smooth the weights
C
         z=0.d0
         DO jw1=1,dlw
            if(jw1.eq.clw) CYCLE
               z=z+w(jw1+trl)
         END DO
         z=(2.d0-z/2.d0)*hw-1+z/2.d0
         z=max(.1d0,min(z,hw))
         cc=min(z-1.d0,1.d0/hakt)
         call smwghts1(w(1+trl),hakt,z,sw(1+trs),dlw,dsw,cc)
         DO k=1,dp2
            swj(k)=0.d0
            swj2(k)=0.d0
            swj0(k)=0.d0
         END DO
         DO k=1,dp1
               swjy(k)=0.d0
         END DO
         DO jw1=csw,dsw
            j1=jw1-csw+iind
            if(j1.lt.1.or.j1.gt.n) CYCLE
            z1=jw1-csw
            lwj=slw(jw1)
            wj=sw(jw1+trs)
            if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
            zz(2)=z1
            zz(3)=z1*z1
            IF(dp1.gt.2) THEN
               zz(4)=z1*zz(3)
               zz(5)=z1*zz(4)
            END IF
            DO k=1,dp2
               swj0(k)=swj0(k)+lwj*zz(k)
            END DO
            if(wj.le.0.d0) CYCLE
            DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
               swj2(k)=swj2(k)+wj*wj*zz(k)
            END DO
            yj=y(j1)
            DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
            END DO
         END DO
         DO jw1=1,csw-1
            j1=jw1-csw+iind
            if(j1.lt.1.or.j1.gt.n) CYCLE
            z1=jw1-csw
            lwj=slw(jw1)
            wj=sw(jw1+trs)
            if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
            zz(2)=z1
            zz(3)=z1*z1
            IF(dp1.gt.2) THEN
               zz(4)=z1*zz(3)
               zz(5)=z1*zz(4)
            END IF
            DO k=1,dp2
               swj0(k)=swj0(k)+lwj*zz(k)
            END DO
            if(wj.le.0.d0) CYCLE
            DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
               swj2(k)=swj2(k)+wj*wj*zz(k)
            END DO
            yj=y(j1)
            DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
            END DO
         END DO
         DO k=1,dp1
            ai(iind+(k-1)*n)=swjy(k)
         END DO
         DO k=1,dp2
            bi(iind+(k-1)*n)=swj(k)
            bi2(iind+(k-1)*n)=swj2(k)
            bi0(iind+(k-1)*n)=swj0(k)
         END DO
         hhom(iind,1)=hhommin
         hhom(iind,2)=hhommax
         if(lfix.and.hakt-hfixmax.ge.hnfix) THEN
            fix(iind)=1
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded), 1D, heteroskedastic
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsph1(y,si,fix,nfix,n,degr,hw,hakt,hhom,lambda,
     1        theta,bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   kern     specifies the location kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C
C   temporary arrays set for maximum degree 2
C
      use omp_lib
      implicit none
      external kldistp,lkern
      double precision kldistp,lkern
      integer n,kern,degr,ind(*),nfix,fix(*)
      logical aws,lfix
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,hhom(*),
     1       bi2(*),hakt,lw(*),w(*),hw,sw(*),slw(*),si(*),spmin
      integer ih,j1,k,iind,jind,dlw,clw,jw1,
     2        dp1,dp2,ihs,csw,dsw,l,thrednr,trl,trs
      double precision bii(5),sij,swj(5),swj2(5),swj0(5),swjy(5),z1,wj,
     1       hakt2,thij(3),thi(3),zz(5),lwj,yj,hs2,hs,z,cc,spf,hhomi,
     2       hhommax,az1,hfixmax,hnfix

C!$    integer omp_get_thread_num
C!$    external omp_get_thread_num

C   arrays with variable length are organized as
C   theta(n,dp1)
C   bi(n,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      lfix=nfix.gt.(degr+1)
      hnfix=nfix
      hnfix=max(hnfix,.2d0*hakt)
      spf=1.d0/(1.d0-spmin)
      if(degr.eq.0) THEN
         dp1=1
         dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=2
         dp2=3
      ELSE
         dp1=3
         dp2=5
      END IF
      hakt2=hakt*hakt
      ih=FLOOR(hakt)
      dlw=2*ih+1
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=FLOOR(hs)
      dsw=2*ihs+1
      csw=ihs+1
C   compute location weights first
      DO j1=1,dlw
         z1=clw-j1
         lw(j1)=lkern(kern,z1*z1/hakt2)
      END DO
      cc=0.0d0
      call smwghts1(lw,hakt,hw,slw,dlw,dsw,cc)
      thrednr = 0
C  now stochastic term
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(y,si,fix,nfix,n,degr,hw,hakt,hhom,lambda,theta,
C$OMP&        bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C$OMP& FIRSTPRIVATE(aws,lfix,hnfix,ih,spf,dp1,dp2,hakt2,dlw,clw,
C$OMP&        hs2,hs,ihs,csw,dsw)
C$OMP& PRIVATE(j1,k,iind,jind,jw1,l,bii,sij,swj,swj2,swj0,swjy,z1,
C$OMP&         wj,thij,thi,zz,lwj,yj,z,cc,hhommax,az1,hfixmax,
C$OMP&         hhomi,thrednr,trl,trs)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
         IF (fix(iind).ne.0) CYCLE
C    nothing to do, final estimate is already fixed by control
!$       thrednr = omp_get_thread_num()
         trl = thrednr*dlw
         trs = thrednr*dsw
         zz(1)=1.d0
         hhomi=hhom(iind)
         hhomi=hhomi
         hhommax=hakt
         hfixmax=hhomi
         DO k=1,dp1
            thi(k)=theta(iind+(k-1)*n)
         END DO
         DO k=1,dp2
            bii(k)=bi(iind+(k-1)*n)/lambda
         END DO
C   scaling of sij outside the loop
         DO jw1=1,dlw
            w(jw1+trl)=0.d0
            jind=jw1-clw+iind
            if(jind.lt.1.or.jind.gt.n) CYCLE
            wj=lw(jw1)*si(jind)
            z1=jw1-clw
            zz(2)=z1
            zz(3)=z1*z1
            az1=abs(z1)
            IF (aws.and.az1.ge.hhomi) THEN
               DO k=1,dp1
                  thij(k)=theta(jind+(k-1)*n)
               END DO
               thij(1)=thij(1)-thij(2)*z1
               IF (dp1.gt.2) THEN
                  thij(1)=thij(1)+thij(3)*zz(3)
                  thij(2)=thij(2)-2.d0*thij(3)*z1
               END IF
C
C           get difference of thetas
C
               DO k=1,dp1
                  thij(k)=thi(k)-thij(k)
               END DO
               sij=kldistp(dp1,thij,bii,ind)
               IF (sij.le.1.d0) THEN
                  hfixmax=max(hfixmax,az1)
                  IF (sij.gt.spmin) THEN
                     w(jw1+trl)=wj*(1.d0-spf*(sij-spmin))
                      hhommax=min(hhommax,az1)
                  ELSE
                     w(jw1+trl)=wj
                  END IF
               ELSE
                  w(jw1+trl)=0.d0
                  hhommax=min(hhommax,az1)
               END IF
            ELSE
               w(jw1+trl)=wj
            END IF
         END DO
C
C      Smooth the weights
C
         z=0.d0
         DO jw1=1,dlw
            if(jw1.eq.clw) CYCLE
               z=z+w(jw1+trl)
         END DO
         z=(2.d0-z/2.d0)*hw-1+z/2.d0
         z=max(.1d0,min(z,hw))
         cc=min(z-1.d0,1.d0/hakt)
         call smwghts1(w(1+trl),hakt,z,sw(1+trs),dlw,dsw,cc)
         DO k=1,dp2
            swj(k)=0.d0
            swj2(k)=0.d0
            swj0(k)=0.d0
         END DO
         DO k=1,dp1
               swjy(k)=0.d0
         END DO
         DO jw1=1,dsw
            j1=jw1-csw+iind
            if(j1.lt.1.or.j1.gt.n) CYCLE
            z1=jw1-csw
            lwj=slw(jw1)
            wj=sw(jw1+trs)
            if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
            zz(2)=z1
            zz(3)=z1*z1
            IF(dp1.gt.2) THEN
               zz(4)=z1*zz(3)
               zz(5)=z1*zz(4)
            END IF
            DO k=1,dp2
               swj0(k)=swj0(k)+lwj*zz(k)
            END DO
            if(wj.le.0.d0) CYCLE
            DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
               swj2(k)=swj2(k)+wj*wj*zz(k)
            END DO
            yj=y(j1)
            DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
            END DO
         END DO
         DO k=1,dp1
            ai(iind+(k-1)*n)=swjy(k)
         END DO
         DO k=1,dp2
            bi(iind+(k-1)*n)=swj(k)
            bi2(iind+(k-1)*n)=swj2(k)
            bi0(iind+(k-1)*n)=swj0(k)
         END DO
         hhom(iind)=hhommax
         if(lfix.and.hakt-hfixmax.ge.hnfix) THEN
            fix(iind)=1
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded), 2D
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsp2(y,fix,nfix,n1,n2,degr,hw,hakt,hhom,lambda,
     1              theta,bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   kern     specifies the location kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C
C   temporary arrays set for maximum degree 2
C
      use omp_lib
      implicit none
      external kldistp,lkern
      double precision kldistp,lkern
      integer n1,n2,kern,degr,ind(*),nfix,fix(*)
      logical aws,lfix
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,spmin,
     1       bi2(*),hakt,lw(*),w(*),hw,sw(*),slw(*),hhom(*)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        iind,jind,jind2,jwind,jwind2,dlw,clw,jw1,jw2,
     2        dp1,dp2,ihs,csw,dsw,l,dlw2,thrednr,trl,trs
      double precision bii(15),sij,swj(15),swj2(15),swj0(15),swjy(6),
     1       z1,z2,wj,hakt2,thij(6),thi(6),zz(15),lwj,hs2,hs,z,cc,
     2       wjy,spf,hhomi,hhommax,az1,hfixmax,hnfix
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
C   arrays with variable length are organized as
C   theta(n1,n2,dp1)
C   bi(n1,n2,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      lfix=nfix.gt.1
      hnfix=max(1.5d0,6.d0+degr-hakt)
      spf=1.d0/(1.d0-spmin)
      if(degr.eq.0) THEN
         dp1=1
         dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=3
         dp2=6
      ELSE
         dp1=6
         dp2=15
      END IF
      hakt2=hakt*hakt
      ih=FLOOR(hakt)
      dlw=2*ih+1
      clw=ih+1
      dlw2=dlw*dlw
      hs=hakt+hw
      hs2=hs*hs
      ihs=FLOOR(hs)
      dsw=2*ihs+1
      csw=ihs+1
      n=n1*n2
C   compute location weights first  sum in slw
      DO j2=1,dlw
         z2=j2-clw
         z2=z2*z2
         ih1=FLOOR(sqrt(hakt2-z2))
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=j1-clw
            lw(jind)=lkern(kern,(z1*z1+z2)/hakt2)
         END DO
      END DO
      cc=0.0d0
      call smwghts2(lw,hakt,hw,slw,dlw,dsw,cc)
C  now stochastic term
      call rchkusr()
      thrednr = 0
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(y,fix,nfix,n1,n2,degr,hw,hakt,hhom,lambda,theta,
C$OMP&        bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C$OMP& FIRSTPRIVATE(aws,lfix,hnfix,ih,spf,dp1,dp2,hakt2,dlw,clw,
C$OMP&        hs2,hs,ihs,csw,dsw,dlw2,n)
C$OMP& PRIVATE(j1,j2,k,iind,jind,jind2,jw1,l,bii,sij,swj,swj2,swj0,
C$OMP&         swjy,z1,wj,thij,thi,zz,lwj,z,cc,hhommax,
C$OMP&         az1,hfixmax,thrednr,trl,trs,
C$OMP&         hhomi,jwind,jwind2,ih1,i1,i2,z2,wjy)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=(iind-i1)/n1+1
            IF (fix(iind).ne.0) CYCLE
C    nothing to do, final estimate is already fixed by control
            zz(1)=1.d0
!$          thrednr = omp_get_thread_num()
            trl = thrednr*dlw2
            trs = thrednr*dsw*dsw
            hhomi=hhom(iind)
            hhomi=hhomi*hhomi-1
            hhommax=hakt2
            hfixmax=hhomi
            DO k=1,dp2
               bii(k)=bi(iind+(k-1)*n)/lambda
            END DO
            DO k=1,dp1
               thi(k)=theta(iind+(k-1)*n)
            END DO
C   scaling of sij outside the loop
            DO jw1=1,dlw2
               w(jw1+trl)=0.d0
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=jw2-clw
C  get directional differences that only depend on i2-j2
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               ih1=FLOOR(sqrt(hakt2-z2*z2))
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jwind=jw1+jwind2
                  wj=lw(jwind)
                  z1=jw1-clw
                  az1=z1*z1+z2*z2
C  get rest of directional differences
                  IF (aws.and.az1.gt.hhomi) THEN
                     IF(dp1.gt.1) THEN
                        zz(2)=z1
                        zz(4)=z1*z1
                        zz(5)=z1*z2
                     END IF
                     DO k=1,dp1
                        thij(k)=theta(jind+(k-1)*n)
                     END DO
                     IF(dp1.gt.1) THEN
                        thij(1)=thij(1)-thij(2)*z1-thij(3)*z2
                        IF (dp1.gt.3) THEN
                           thij(1)=thij(1)+thij(4)*zz(4)+thij(5)*
     1                             zz(5)+thij(6)*zz(6)
                           thij(2)=thij(2)-thij(5)*z2-2.d0*thij(4)*z1
                           thij(3)=thij(3)-thij(5)*z1-2.d0*thij(6)*z2
                        END IF
                     END IF
C
C           get difference of thetas
C
                     DO k=1,dp1
                        thij(k)=thi(k)-thij(k)
                     END DO
                     sij=kldistp(dp1,thij,bii,ind)
                     IF (sij.le.1.d0) THEN
                        hfixmax=max(hfixmax,az1)
                        IF (sij.gt.spmin) THEN
                           w(jwind+trl)=wj*(1.d0-spf*(sij-spmin))
                           hhommax=min(hhommax,az1)
                        ELSE
                           w(jwind+trl)=wj
                        END IF
                     ELSE
                        w(jwind+trl)=0.d0
                        hhommax=min(hhommax,az1)
                     END IF
                  ELSE
                     w(jwind+trl)=wj
                  END IF
               END DO
            END DO
C
C      Smooth the weights
C
            call testwgh2(w(1+trl),dlw,dp1,hw,z)
            z=max(.1d0,min(z,hw))
            cc=min(z-1.d0,1.d0/hakt2)
            call smwghts2(w(1+trl),hakt,z,sw(1+trs),dlw,dsw,cc)
            DO k=1,dp2
               swj(k)=0.d0
               swj2(k)=0.d0
               swj0(k)=0.d0
            END DO
            DO k=1,dp1
               swjy(k)=0.d0
            END DO
            DO jw2=1,dsw
               j2=jw2-csw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dsw
               z2=jw2-csw
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               IF(dp1.gt.3) THEN
                  zz(10)=z2*zz(6)
                  zz(15)=z2*zz(10)
               END IF
               ih1=FLOOR(sqrt(hs2-z2*z2))
               DO jw1=csw-ih1,csw+ih1
                  j1=jw1-csw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  z1=jw1-csw
                  jwind=jw1+jwind2
                  jind=j1+jind2
                  lwj=slw(jwind)
                  wj=sw(jwind+trs)
                  IF(wj.lt.0.d0.or.wj.gt.1.d0) THEN
                     wj=max(0.d0,min(1.d0,wj))
                  END IF
                  if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
                  IF(dp1.gt.1) THEN
                     zz(2)=z1
                     zz(4)=z1*z1
                     zz(5)=z1*z2
                  END IF
                  IF(dp1.gt.3) THEN
                     zz(7)=z1*zz(4)
                     zz(8)=z1*zz(5)
                     zz(9)=z1*zz(6)
                     zz(11)=z1*zz(7)
                     zz(12)=z1*zz(8)
                     zz(13)=z1*zz(9)
                     zz(14)=z1*zz(10)
                  END IF
                  DO k=1,dp2
                     swj0(k)=swj0(k)+lwj*zz(k)
                  END DO
                  if(wj.le.0.d0) CYCLE
                  DO k=1,dp2
                     swj(k)=swj(k)+wj*zz(k)
                     swj2(k)=swj2(k)+wj*wj*zz(k)
                  END DO
                  wjy=wj*y(jind)
                  DO l=1,dp1
                     swjy(l)=swjy(l)+wjy*zz(l)
                  END DO
               END DO
            END DO
            DO k=1,dp1
               ai(iind+(k-1)*n)=swjy(k)
            END DO
            DO k=1,dp2
               bi(iind+(k-1)*n)=swj(k)
               bi2(iind+(k-1)*n)=swj2(k)
               bi0(iind+(k-1)*n)=swj0(k)
            END DO
            hhom(iind)=sqrt(hhommax)
            IF(lfix.and.hakt-sqrt(hfixmax).ge.hnfix) THEN
               fix(iind)=1
            END IF
C         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom,fix)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded), 2D, heteroskedastic
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsph2(y,si,fix,nfix,n1,n2,degr,hw,hakt,hhom,lambda,
     1              theta,bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   kern     specifies the location kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C
C   temporary arrays set for maximum degree 2
C
      use omp_lib
      implicit none
      external kldistp,lkern
      double precision kldistp,lkern
      integer n1,n2,kern,degr,ind(*),nfix,fix(*)
      logical aws,lfix
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,spmin,
     1       bi2(*),hakt,lw(*),w(*),hw,sw(*),slw(*),si(*),hhom(*)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        iind,jind,jind2,jwind,jwind2,dlw,clw,jw1,jw2,
     2        dp1,dp2,ihs,csw,dsw,l,dlw2,thrednr,trl,trs
      double precision bii(15),sij,swj(15),swj2(15),swj0(15),swjy(6),
     1       z1,z2,wj,hakt2,thij(6),thi(6),zz(15),lwj,hs2,hs,z,cc,
     2       wjy,spf,hhomi,hhommax,az1,hfixmax,hnfix
C!$    integer omp_get_thread_num
C!$    external omp_get_thread_num
C   arrays with variable length are organized as
C   theta(n1,n2,dp1)
C   bi(n1,n2,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      lfix=nfix.gt.1
      hnfix=max(1.5d0,6.d0+degr-hakt)
      spf=1.d0/(1.d0-spmin)
      if(degr.eq.0) THEN
         dp1=1
         dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=3
         dp2=6
      ELSE
         dp1=6
         dp2=15
      END IF
      hakt2=hakt*hakt
      ih=FLOOR(hakt)
      dlw=2*ih+1
      clw=ih+1
      dlw2=dlw*dlw
      hs=hakt+hw
      hs2=hs*hs
      ihs=FLOOR(hs)
      dsw=2*ihs+1
      csw=ihs+1
      n=n1*n2
C   compute location weights first  sum in slw
      DO j2=1,dlw
         z2=j2-clw
         z2=z2*z2
         ih1=FLOOR(sqrt(hakt2-z2))
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=j1-clw
            lw(jind)=lkern(kern,(z1*z1+z2)/hakt2)
         END DO
      END DO
      cc=0.0d0
      call smwghts2(lw,hakt,hw,slw,dlw,dsw,cc)
C  now stochastic term
      call rchkusr()
      thrednr = 0
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(y,si,fix,nfix,n1,n2,degr,hw,hakt,hhom,lambda,theta,
C$OMP&        bi,bi2,bi0,ai,kern,spmin,lw,w,slw,sw,ind)
C$OMP& FIRSTPRIVATE(aws,lfix,hnfix,ih,spf,dp1,dp2,hakt2,dlw,clw,
C$OMP&        hs2,hs,ihs,csw,dsw,dlw2,n)
C$OMP& PRIVATE(j1,j2,k,iind,jind,jind2,jw1,l,bii,sij,swj,swj2,swj0,
C$OMP&         swjy,z1,wj,thij,thi,zz,lwj,z,cc,hhommax,
C$OMP&         az1,hfixmax,thrednr,trl,trs,
C$OMP&         hhomi,jwind,jwind2,ih1,i1,i2,z2,wjy)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=(iind-i1)/n1+1
            IF (fix(iind).ne.0) CYCLE
C    nothing to do, final estimate is already fixed by control
            zz(1)=1.d0
!$            thrednr = omp_get_thread_num()
            trl = thrednr*dlw2
            trs = thrednr*dsw*dsw
            hhomi=hhom(iind)
            hhomi=hhomi*hhomi
            hhommax=hakt2
            hfixmax=hhomi
            DO k=1,dp2
               bii(k)=bi(iind+(k-1)*n)/lambda
            END DO
            DO k=1,dp1
               thi(k)=theta(iind+(k-1)*n)
            END DO
C   scaling of sij outside the loop
            DO jw1=1,dlw2
               w(jw1+trl)=0.d0
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=jw2-clw
C  get directional differences that only depend on i2-j2
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               ih1=FLOOR(sqrt(hakt2-z2*z2))
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jwind=jw1+jwind2
                  wj=lw(jwind)
                  z1=jw1-clw
C  get rest of directional differences
                  az1=z1*z1+z2*z2
C  get rest of directional differences
                  IF (aws.and.az1.ge.hhomi) THEN
                     IF(dp1.gt.1) THEN
                        zz(2)=z1
                        zz(4)=z1*z1
                        zz(5)=z1*z2
                     END IF
                     DO k=1,dp1
                        thij(k)=theta(jind+(k-1)*n)
                     END DO
                     IF(dp1.gt.1) THEN
                        thij(1)=thij(1)-thij(2)*z1-thij(3)*z2
                        IF (dp1.gt.3) THEN
                           thij(1)=thij(1)+thij(4)*zz(4)+thij(5)*
     1                             zz(5)+thij(6)*zz(6)
                           thij(2)=thij(2)-thij(5)*z2-2.d0*thij(4)*z1
                           thij(3)=thij(3)-thij(5)*z1-2.d0*thij(6)*z2
                        END IF
                     END IF
C
C           get difference of thetas
C
                     DO k=1,dp1
                        thij(k)=thi(k)-thij(k)
                     END DO
                     sij=kldistp(dp1,thij,bii,ind)
                     IF (sij.le.1.d0) THEN
                        hfixmax=max(hfixmax,az1)
                        IF (sij.gt.spmin) THEN
                           w(jwind+trl)=wj*(1.d0-spf*(sij-spmin))
                           hhommax=min(hhommax,az1)
                        ELSE
                           w(jwind+trl)=wj
                        END IF
                     ELSE
                        w(jwind+trl)=0.d0
                        hhommax=min(hhommax,az1)
                     END IF
                  ELSE
                  w(jwind+trl)=wj
                  END IF
               END DO
            END DO
C
C      Smooth the weights
C
            call testwgh2(w(1+trl),dlw,dp1,hw,z)
            z=max(.1d0,min(z,hw))
            cc=min(z-1.d0,1.d0/hakt2)
            call smwghts2(w(1+trl),hakt,z,sw(1+trs),dlw,dsw,cc)
            DO k=1,dp2
               swj(k)=0.d0
               swj2(k)=0.d0
               swj0(k)=0.d0
            END DO
            DO k=1,dp1
               swjy(k)=0.d0
            END DO
            DO jw2=1,dsw
               j2=jw2-csw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dsw
               z2=jw2-csw
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               IF(dp1.gt.3) THEN
                  zz(10)=z2*zz(6)
                  zz(15)=z2*zz(10)
               END IF
               ih1=FLOOR(sqrt(hs2-z2*z2))
               DO jw1=csw-ih1,csw+ih1
                  j1=jw1-csw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  z1=jw1-csw
                  jwind=jw1+jwind2
                  jind=j1+jind2
                  lwj=slw(jwind)*si(jind)
                  wj=sw(jwind+trs)*si(jind)
                  if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
                  IF(dp1.gt.1) THEN
                     zz(2)=z1
                     zz(4)=z1*z1
                     zz(5)=z1*z2
                  END IF
                  IF(dp1.gt.3) THEN
                     zz(7)=z1*zz(4)
                     zz(8)=z1*zz(5)
                     zz(9)=z1*zz(6)
                     zz(11)=z1*zz(7)
                     zz(12)=z1*zz(8)
                     zz(13)=z1*zz(9)
                     zz(14)=z1*zz(10)
                  END IF
                  DO k=1,dp2
                     swj0(k)=swj0(k)+lwj*zz(k)
                  END DO
                  if(wj.le.0.d0) CYCLE
                  DO k=1,dp2
                     swj(k)=swj(k)+wj*zz(k)
                     swj2(k)=swj2(k)+wj*wj*zz(k)
                  END DO
                  wjy=wj*y(jind)
                  DO l=1,dp1
                     swjy(l)=swjy(l)+wjy*zz(l)
                  END DO
               END DO
            END DO
            DO k=1,dp1
               ai(iind+(k-1)*n)=swjy(k)
            END DO
            DO k=1,dp2
               bi(iind+(k-1)*n)=swj(k)
               bi2(iind+(k-1)*n)=swj2(k)
               bi0(iind+(k-1)*n)=swj0(k)
            END DO
            hhom(iind)=sqrt(hhommax)
            IF(lfix.and.hakt-sqrt(hfixmax).ge.hnfix) THEN
               fix(iind)=1
            END IF
C         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom,fix)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function kldistp(dp1,thij,bii,ind)
C
C  distance between local polynomial models
C
C     dp1  polynomial degree +1
C     thij parameter estimate in j dim(dp1) for basis in i
C     bii  XTX dim(dp2)
C     ind   index matrix to access the correct elements in bii
C
      implicit none
      integer dp1,ind(dp1,dp1)
      double precision thij(*),bii(*),thijl
      integer l,k
      double precision d
      d=0.d0
      DO l=1,dp1
         thijl=thij(l)
         d=d+bii(ind(l,l))*thijl*thijl
         IF(l.eq.dp1) CYCLE
         DO k=l+1,dp1
            d=d+2.d0*bii(ind(k,l))*thijl*thij(k)
         END DO
      END DO
      kldistp=d
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghts1(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit none
      integer dw,dsw,cw,csw,cdiff
      double precision w(dw),sw(dsw),hw,hakt,cc
      integer i1,ja1,je1,j1,i10
      double precision z,z0,z1,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      zmax=0.d0
      DO i1=1,dsw
         sw(i1)=0.d0
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
            sw(j1+cdiff)=w(j1)
         END DO
      ELSE
         DO i1=1,dsw
            z1=i1-csw
            i10=i1-cdiff
            ja1=max(i1-2*cdiff,1)
            je1=min(i1,dw)
            z=0.d0
            z0=0.d0
            DO j1=ja1,je1
               z1=(i10-j1)
               z1=z1*z1
               if(hw2-z1.lt.0.d0) CYCLE
               ww=(1.d0-z1/hw2)
               if(ww.lt.1.d0) ww=cc*ww
               z=z+ww*w(j1)
            END DO
            sw(i1)=z
            zmax=max(zmax,z)
         END DO
         DO i1=1,dsw
            sw(i1)=sw(i1)/zmax
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghts2(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit none
      integer dw,dsw,cw,csw,cdiff
      double precision w(dw,dw),sw(dsw,dsw),hw,hakt,cc
      integer i1,i2,id,jd,ja1,je1,ja2,je2,j1,j2,i10,i20
      double precision z,z0,z1,z2,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      zmax=0.d0
      DO i1=1,dsw
         DO i2=1,dsw
            sw(i1,i2)=0.d0
         END DO
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
            DO j2=1,dw
               sw(j1+cdiff,j2+cdiff)=w(j1,j2)
            END DO
         END DO
      ELSE
         DO i1=1,dsw
            z1=i1-csw
            i10=i1-cdiff
            ja1=max(i1-2*cdiff,1)
            je1=min(i1,dw)
            id=FLOOR(sqrt(hsw2-z1*z1))
            if(csw-id.lt.1) CYCLE
            DO i2=csw-id,csw+id
               i20=i2-cdiff
               z=0.d0
               z0=0.d0
               DO j1=ja1,je1
                  z1=(i10-j1)
                  z1=z1*z1
                  if(hw2-z1.lt.0.d0) CYCLE
                  jd=FLOOR(sqrt(hw2-z1))
                  ja2=max(i20-jd,1)
                  je2=min(i20+jd,dw)
                  DO j2=ja2,je2
                     z2=(i20-j2)
                     ww=(1.d0-(z1+z2*z2)/hw2)
                     if(ww.lt.1.d0) ww=cc*ww
                     z=z+ww*w(j1,j2)
                  END DO
               END DO
               sw(i1,i2)=z
               zmax=max(zmax,z)
            END DO
         END DO
         DO i1=1,dsw
            DO i2=1,dsw
               sw(i1,i2)=sw(i1,i2)/zmax
            END DO
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      test regulatity of weight matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine testwgh2(w,dlw,dp1,hw,z)
      integer dlw,dp1
      double precision w(dlw,dlw),hw,z
      integer clw,i,ip,im,cp1,cp2,cm1,cm2
      double precision zh,zv
      clw=(dlw+1)/2
      cp1=clw+1
      cm1=clw-1
      z=hw
      IF(clw.gt.2.and.dp1.eq.3) THEN
         cp2=clw+2
         cm2=clw-2
         zh=w(clw,cp1)*w(clw,cp2)+w(clw,cm1)*w(clw,cm2)
         zv=w(cp1,clw)*w(cp2,clw)+w(cm1,clw)*w(cm2,clw)
         IF(zh*zv.gt.0.125d0) THEN
            z=hw-2.d0
         ELSE
            DO i=1,clw-1
               ip=clw+i
               im=clw-i
               zh=zh+w(ip,cp1)*w(ip,cp2)+w(ip,cm1)*w(ip,cm2)+
     1               w(im,cp1)*w(im,cp2)+w(im,cm1)*w(im,cm2)
               zv=zv+w(cp1,ip)*w(cp2,ip)+w(cm1,ip)*w(cm2,ip)+
     1               w(cp1,im)*w(cp2,im)+w(cm1,im)*w(cm2,im)
               IF(zh*zv.gt.0.125d0) THEN
                  z=hw-2.d0
                  CYCLE
               END IF
            END DO
            IF(zh*zv.le.0.125d0) THEN
               zh=w(clw,cp1)+w(clw,cm1)
               zv=w(cp1,clw)+w(cm1,clw)
               DO i=1,clw-1
                  ip=clw+i
                  im=clw-i
                  zh=zh+w(ip,cp1)+w(ip,cm1)+w(im,cp1)+w(im,cm1)
                  zv=zv+w(cp1,ip)+w(cm1,ip)+w(im,cp1)+w(im,cm1)
                  IF(zh*zv.gt.0.5d0) THEN
                     z=hw-1.d0
                     CYCLE
                  END IF
               END DO
            END IF
         END IF
      ELSE IF(clw.gt.1.and.dp1.eq.2) THEN
         zh=w(clw,cp1)+w(clw,cm1)
         zv=w(cp1,clw)+w(cm1,clw)
         IF(zh*zv.gt.0.5d0) THEN
            z=hw-1.d0
         ELSE
            DO i=1,clw-1
               ip=clw+i
               im=clw-i
               zh=zh+w(ip,cp1)+w(ip,cm1)+w(im,cp1)+w(im,cm1)
               zv=zv+w(cp1,ip)+w(cm1,ip)+w(im,cp1)+w(im,cm1)
               IF(zh*zv.gt.0.5d0) THEN
                  z=hw-1.d0
                  CYCLE
               END IF
            END DO
         END IF
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpaws2(n,dp1,dp2,ai,bi,theta,dmat,ind)
C
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y
C     bi         \sum \Psi^T Wi^k \Psi
C     theta      new parameter estimate
C     dmat       working arrays
C     restricted to dp2<=20
      implicit none
      integer n,dp1,dp2
      double precision ai(n,dp1),bi(n,dp2),theta(n,dp1),dmat(dp1,dp1)
      integer i,j,k,info,ind(dp1,dp1)
      double precision aa(20)
      DO i=1,n
         DO k=1,dp1
            DO j=k,dp1
               dmat(k,j)=bi(i,ind(k,j))
            END DO
            aa(k)=ai(i,k)
         END DO
C     now calculate theta as B_i^{-1} A_i
         call dposv("U",dp1,1,dmat,dp1,aa,dp1,info)
C    if info>0 just keep the old estimate
         IF (info.gt.0) CYCLE
         DO j=1,dp1
            theta(i,j)=aa(j)
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpaws1(n,dp1,dp2,ai,bi,theta,dmat,ind)
C
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y
C     bi         \sum \Psi^T Wi^k \Psi
C     theta      new parameter estimate
C     dmat       working arrays
C     restricted to dp2<=20
      implicit none
      integer n,dp1,dp2
      double precision ai(n,dp1),bi(n,dp2),theta(n,dp1),dmat(dp1,dp1)
      integer i,j,k,info,ind(dp1,dp1),iii
      double precision aa(3),h,cii(5)
      DO i=1,n
         cii(1)=1.d0
         h=bi(i,1)
         if(dp1.gt.1) THEN
            DO k=2,dp2
               cii(k)=cii(k-1)*h
            END DO
         END IF
         DO k=1,dp1
            DO j=k,dp1
               iii=ind(k,j)
               dmat(k,j)=bi(i,iii)/cii(iii)
            END DO
            aa(k)=ai(i,k)/cii(k)
         END DO
C     now calculate theta as B_i^{-1} A_i
         call dposv("U",dp1,1,dmat,dp1,aa,dp1,info)
C    if info>0 just keep the old estimate
         IF (info.gt.0) CYCLE
         DO j=1,dp1
            theta(i,j)=aa(j)/cii(j)
         END DO
      END DO
      RETURN
      END
      subroutine vpaws(n,dp2,bi,bi2,var)
C   calculate variance (reduction) of estimates
      integer n,dp2
      double precision bi(n,dp2),bi2(n,dp2),var(n)
      integer i
      double precision a1,a2,a3,d
      if(dp2.eq.3) THEN
         DO i=1,n
            d = bi(i,1)*bi(i,3)-bi(i,2)*bi(i,2)
            IF(d.lt.1d-8) THEN
               var(i)=1.d0/bi(i,1)
C  bi(i,1) should contain either 1 or 1/sigma^2 in this case
               CYCLE
            END IF
            a1 = bi(i,3)/d
            a2 = -bi(i,2)/d
            var(i)=a1*a1*bi2(i,1)+a2*a2*bi2(i,3)+2.d0*a1*a2*bi2(i,2)
         END DO
      ELSE
         DO i=1,n
            d=bi(i,1)*bi(i,4)*bi(i,6)+2.d0*bi(i,2)*bi(i,3)*bi(i,5)-
     1        bi(i,3)*bi(i,3)*bi(i,4)-bi(i,2)*bi(i,2)*bi(i,6)-
     1        bi(i,5)*bi(i,5)*bi(i,1)
            IF(d.lt.1d-8) THEN
               var(i)=1.d0/bi(i,1)
C  bi(i,1) should contain either 1 or 1/sigma^2 in this case
            CYCLE
            END IF
            a1=(bi(i,6)*bi(i,4)-bi(i,5)*bi(i,5))/d
            a2=(bi(i,3)*bi(i,5)-bi(i,2)*bi(i,6))/d
            a3=(bi(i,2)*bi(i,5)-bi(i,3)*bi(i,4))/d
            var(i)=a1*a1*bi2(i,1)+a2*a2*bi2(i,4)+a3*a3*bi2(i,6)+
     1             2.d0*a1*a2*bi2(i,2)+2.d0*a1*a3*bi2(i,3)+
     1             2.d0*a2*a3*bi2(i,5)
         END DO
      END IF
      RETURN
      END
