      subroutine swap(x,k,l)
      implicit none
      integer k,l
      double precision x(*),t
      t = x(k)
      x(k) = x(l)
      x(l) = t
      return
      end

      subroutine qselect(x,n,k)
C
C  partial sorting using a select algorithm
C  output x(k) contains the kth element of sort(x)
C
      implicit none
      integer k,n
      double precision x(n)
      integer lft,rght,cur,i
      double precision guess
      lft=1
      rght=n
      DO WHILE (lft.lt.rght)
         guess = x(k)
         call swap(x,k,rght)
         cur = lft
         DO i=cur,rght-1
            IF(x(i).lt.guess) THEN
               call swap(x,i,cur)
               cur=cur+1
            END IF
         END DO
         call swap(x,rght,cur)
         if(cur.eq.k) EXIT
         if(cur.lt.k) THEN
            lft=cur+1
         ELSE
            rght=cur-1
         END IF
      END DO
      RETURN
      END

      double precision function fmedian(x,n)
C
C  compute the median using a select algorithm instead of sorting
C  used in mediansm
C
      implicit none
      integer n
      double precision x(n)
      integer m
      m = n/2+1
      call qselect(x,n,m)
      fmedian = x(m)
      if(mod(n,2).eq.0) THEN
         m = n-m+1
         call qselect(x,n,m)
         fmedian = (fmedian+x(m))/2.d0
      END IF
      return
      end
      subroutine mediansm(y,mask,n1,n2,n3,ind,nind,work,ncores,yout)
C
C
C   3D median smoother of y with neighborhood defined by ind
C   results in yout
C   size of work needs to be 2*nind
C
      use omp_lib
      implicit none
      integer n1,n2,n3,nind,ind(3,nind),ncores
      integer mask(n1,n2,n3)
      double precision y(n1,n2,n3),yout(n1,n2,n3),work(nind,ncores)
      integer i1,i2,i3,j1,j2,j3,j,k,thrednr
      double precision fmedian
      external fmedian
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      thrednr = 1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i1,i2,i3,j1,j2,j3,j,k,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO i1=1,n1
!$         thrednr = omp_get_thread_num()+1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3).eq.0) THEN
                  yout(i1,i2,i3) = y(i1,i2,i3)
                  CYCLE
               ENDIF
               k=0
               DO j=1,nind
                  j1=i1+ind(1,j)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  j2=i2+ind(2,j)
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  j3=i3+ind(3,j)
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  if(mask(j1,j2,j3).eq.0) CYCLE
                  k=k+1
                  work(k,thrednr)=y(j1,j2,j3)
               END DO
               IF(k.gt.1) THEN
                  yout(i1,i2,i3) = fmedian(work(1,thrednr),k)
               ELSE
                  yout(i1,i2,i3) = y(i1,i2,i3)
               END IF
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(yout)
      return
      end

      subroutine awslchi2(s,ksi,ni,sigma,vpar,L,mask,n1,n2,n3,ind,
     1      w,nw,minni,wad,sad,lambda,nthreds,iL,work,thn,sigman,ksin,
     2      flb,nfb)
C
C  local variance estimation using (adaptive) weighted likelihood
C
C   Takes observed intensities in s and
C     initial estimates of \sigma in sigma
C   perform adaptive smoothing on R^3
C   th containes previous estimates of E S
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value
C
C   lambda   - kritical value for pairwise tests
C   thn      - new estimate sum_j w_a(j) S_j
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:3] are j1-i1,j2-i2 and j3-i3 respectively
C   wad, sad - array for weights>0 and corresponding observed s
C
      use omp_lib
      implicit none
      integer n1,n2,n3,nw,ind(3,nw),nthreds,iL,nfb
      integer mask(n1,n2,n3)
      double precision s(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1  ksi(n1,n2,n3),sigman(n1*n2*n3),lambda,w(nw),sigma(n1,n2,n3),
     2  wad(nw,nthreds),sad(nw,nthreds),L,minni,work(iL,nthreds),
     3  ksin(n1,n2,n3),vpar(6),flb(nfb)
      integer i1,i2,i3,j1,j2,j3,i,j,jj,n,maxit,thrednr
      double precision z,sw,sws,sws2,sj,thi,wj,kval,fnsi,sgi,tol,low,up,
     1       fmin,sgi2,vz,thi2,thj2,fnsj,thj,nii
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
      tol=1d-5
      maxit=100
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& FIRSTPRIVATE(iL,L,minni,n1,n2,n3,maxit)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,sws,sws2,thi,kval,thi2,thj2,
C$OMP& wj,sj,thrednr,fnsi,low,up,tol,sgi,jj,fmin,sgi2,vz,thj,fnsj,nii)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         if(mask(i1,i2,i3).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
         sw=0.d0
         sws=0.d0
         sws2=0.d0
         sgi=sigma(i1,i2,i3)
         sgi2=sgi*sgi
         thi = sqrt(max(1d-16,ksi(i1,i2,i3)/sgi2-2.d0*L))
         thn(i) = thi
         if(thi.gt.vpar(1)) THEN
            thi2 = thi*thi
            vz = vpar(3)*thi+vpar(4)*thi2+vpar(5)*thi*thi2+vpar(6)
            fnsi = max(vpar(2),vz/(vz+1))
         ELSE
            fnsi = vpar(2)
         END IF
         nii = ni(i)
C   thats the estimated standard deviation of s(i1,i2,i3)
         kval = lambda/nii*(sgi/nii+thi)/(0.1d0/nii*sgi+thi)
C allow for increased kval for low SNR and small ni
C correction vanishes for nii -> \infty and thi>0
C without this the propagation condition is violated for very small SNR
         jj = 0
         DO j=1,nw
            wad(j,thrednr)=0.d0
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(mask(j1,j2,j3).eq.0) CYCLE
            wj=w(j)
            thj = sqrt(max(1d-16,ksi(j1,j2,j3)/sgi2-2.d0*L))
            if(thj.gt.vpar(1)) THEN
               thj2 = thj*thj
               vz = vpar(3)*thj+vpar(4)*thj2+vpar(5)*thj*thj2+vpar(6)
               fnsj = max(vpar(2),vz/(vz+1))
            ELSE
               fnsj = vpar(2)
            END IF
            z=thi-thj
            z=z*z/(fnsi+fnsj)
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            sw=sw+wj
            sj=s(j1,j2,j3)
            sws=sws+wj*sj
            sws2=sws2+wj*sj*sj
            jj=jj+1
            wad(jj,thrednr)=wj
            sad(jj,thrednr)=sj
         END DO
         ni(i) = sw
         if(sw.gt.minni) THEN
            ksin(i1,i2,i3) = sws2/sw
C needed for the next iteration
C            low = max(sqrt(ksin(i1,i2,i3)/2.d0/L),sgi/1d1)
             low = sgi/2d0
C  sqrt(ksin(i1,i2,i3)/2.d0/L) is the solution in the central case !
C  old code was still correct but inefficient
C            up = sgi*1d1
            up = min(sqrt(ksin(i1,i2,i3)/2.d0/L),sgi*2d0)
            if(up.le.low) THEN
               sgi = up
            ELSE
               call localmin(low,up,wad(1,thrednr),sad(1,thrednr),L,jj,
     1                    tol,maxit,work(1,thrednr),sgi,fmin,flb,nfb)
            END IF
         END IF
         sigman(i)=sgi
         thn(i) = sqrt(max(1.d-16,ksin(i1,i2,i3)-2.d0*sgi*sgi*L))
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sigman)
      RETURN
      END
C
C   same for Gaussian distribution
C
      subroutine awslgaus(s,th,ni,sigma,mask,n1,n2,n3,ind,
     1      w,nw,minni,lambda,thn,sigman)
C
C  local variance estimation for Gaussian data
C  using (adaptive) weighted likelihood
C
C   Takes observed intensities in s and
C     initial estimates of \sigma in sigma
C   perform adaptive smoothing on R^3
C   th containes previous estimates of E S
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value
C
C   lambda   - kritical value for pairwise tests
C   thn      - new estimate sum_j w_a(j) S_j
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:3] are j1-i1,j2-i2 and j3-i3 respectively
C   wad, sad - array for weights>0 and corresponding observed s
C
      use omp_lib
      implicit none
      integer n1,n2,n3,nw,ind(3,nw)
      integer mask(n1,n2,n3)
      double precision s(n1,n2,n3),ni(n1*n2*n3),thn(n1*n2*n3),
     1 th(n1,n2,n3),sigman(n1*n2*n3),lambda,w(nw),sigma(n1,n2,n3),minni
      integer i1,i2,i3,j1,j2,j3,i,j,n,thrednr
      double precision z,sw,sws,sws2,sj,thi,wj,kval,sgi
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,sws,sws2,thi,kval,
C$OMP& wj,sj,thrednr,sgi)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         if(mask(i1,i2,i3).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
         sw=0.d0
         sws=0.d0
         sws2=0.d0
         sgi=sigma(i1,i2,i3)
         thi = th(i1,i2,i3)
C   thats the estimated standard deviation of s(i1,i2,i3)
         kval = 2.d0*lambda/ni(i)
         DO j=1,nw
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            if(mask(j1,j2,j3).eq.0) CYCLE
            wj=w(j)
            z=(thi-th(j1,j2,j3))/sgi
            z=z*z
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            sw=sw+wj
            sj=s(j1,j2,j3)
            sws=sws+wj*sj
            sws2=sws2+wj*sj*sj
         END DO
         ni(i) = sw
         if(sw.gt.minni) THEN
          sigman(i)= sqrt(max(0.d0,(sws2-sws/sw*sws)/(sw-1.d0)))
C  using sigma not sigma^2
         ELSE
            sigman(i)=sgi
         END IF
         thn(i) = sws/sw
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sigman)
      RETURN
      END

      subroutine awsvchi(y,th,ni,fns,mask,n1,n2,n3,ind,w,nw,lambda,
     1                    sigma,thn,sy)
C   Takes noncentral Chi values in y
C   perform adaptive smoothing on R^3
C   th containes previous estimates
C   ni containes previous sum of weights divided by variance of Chi(2,th/sigma)
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value
C   lambda   - kritical value for pairwise tests
C   sigma    - actual estimate of sigma
C   thn      - new estimate sum_j w_a(j) Y_j
C   th2      - sum_j w_a(j) Y_j^2
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
      implicit none
      integer n1,n2,n3,nw,ind(3,nw)
      integer mask(n1,n2,n3)
      double precision y(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),
     1  thn(n1*n2*n3),sy(n1*n2*n3),lambda,w(nw),sigma,fns(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,i,j,n
      double precision z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw,fnsi
      n = n1*n2*n3
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw,fnsi)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         if(mask(i1,i2,i3).eq.0) CYCLE
         sw=0.d0
         swy=0.d0
         sw2=0.d0
         swy2=0.d0
         thi = th(i1,i2,i3)
         fnsi = fns(i1,i2,i3)
C   thats the estimated standard deviation of y(i1,i2,i3)
         kval = lambda/ni(i)*sigma*sigma
         Do j=1,nw
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            wj=w(j)
            z=thi-th(j1,j2,j3)
            z=z*z/(fnsi+fns(j1,j2,j3))
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            sw=sw+wj
            sw2=sw2+wj*wj
            yj=y(j1,j2,j3)
            swy=swy+wj*yj
            swy2=swy2+wj*yj*yj
         END DO
         thi = swy/sw
         z = swy2/sw
C  z2-thi^2  is an estimate of the variance of y(i)
         cw = 1.d0-sw2/sw/sw
         IF(cw.gt.0.d0) THEN
            sy(i) = sqrt((z-thi*thi)/cw)
C  sy(i)  is an estimate of sigma corrected for
C       simultaneously estimating the mean and for non-central chi-bias
         ELSE
            sy(i) = 0.d0
C  case ni(i) = 1
         END IF
         thn(i) = thi
         ni(i) = sw
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sy)
      RETURN
      END
      subroutine awsadchi(y,th,ni,fns,mask,n1,n2,n3,ind,w,nw,lambda,
     1                    sigma,wad,nthreds,thn,sy)
C   Takes noncentral Chi values in y
C   perform adaptive smoothing on R^3
C   th containes previous estimates
C   ni containes previous sum of weights
C   mask - logical mask (use if mask==TRUE)
C   n1,n2,n3 - dimensions
C   ind  - integer array dim (3,n) containing relative indices in xyz
C   w    - vector of corresponding location weights
C   nw   - number of positive weights (initial value
C   lambda   - kritical value for pairwise tests
C   sigma    - actual estimate of sigma
C   thn      - new estimate sum_j w_a(j) Y_j
C   th2      - sum_j w_a(j) Y_j^2
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
      use omp_lib
      implicit none
      integer n1,n2,n3,nw,ind(3,nw),nthreds
      integer mask(n1,n2,n3)
      double precision y(n1,n2,n3),th(n1,n2,n3),ni(n1*n2*n3),
     1       thn(n1*n2*n3),fns(n1,n2,n3),sy(n1*n2*n3),lambda,w(nw),
     2       sigma,wad(nw,nthreds)
      integer i1,i2,i3,j1,j2,j3,i,j,n,thrednr
      double precision z,sw,sw2,swy,swy2,yj,thi,wj,kval,cw,fnsi
C!$      integer omp_get_thread_num
C!$      external omp_get_thread_num
      n = n1*n2*n3
      thrednr = 1
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,i1,i2,i3,j1,j2,j3,z,sw,swy,sw2,swy2,thi,kval,
C$OMP& wj,yj,cw,thrednr,fnsi)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         if(mask(i1,i2,i3).eq.0) CYCLE
!$         thrednr = omp_get_thread_num()+1
         sw=0.d0
         swy=0.d0
         sw2=0.d0
         swy2=0.d0
         thi = th(i1,i2,i3)
         fnsi = fns(i1,i2,i3)
C   thats the estimated standard deviation of y(i1,i2,i3)
         kval = lambda/ni(i)*sigma*sigma
         DO j=1,nw
            wad(j,thrednr)=0.d0
            j1=i1+ind(1,j)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,j)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,j)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            wj=w(j)
            z=thi-th(j1,j2,j3)
            z=z*z/(fnsi+fns(j1,j2,j3))
            if(z.ge.kval) CYCLE
            wj=wj*min(1.d0,2.d0-2.d0*z/kval)
            wad(j,thrednr)=wj
            sw=sw+wj
            sw2=sw2+wj*wj
            yj=y(j1,j2,j3)
            swy=swy+wj*yj
         END DO
         thi = swy/sw
         DO j=1,nw
            wj=wad(j,thrednr)
            if(wj.le.1d-8) CYCLE
            j1=i1+ind(1,j)
            j2=i2+ind(2,j)
            j3=i3+ind(3,j)
C no need to test for grid coordinates since wj>0
            swy2=swy2+wj*abs(thi-y(j1,j2,j3))
         END DO
         z = swy2/sw/.8d0
C  z  is an estimate of the standard deviation of y(i) by mean absolute deviation
         cw = 1.d0-sw2/sw/sw
         IF(cw.gt.0.d0) THEN
            sy(i) = z/sqrt(cw)
C  sy(i)  is an estimate of sigma corrected for
C       simultaneously estimating the mean and for non-central chi-bias
         ELSE
            sy(i) = 0.d0
C  case ni(i) = 1
         END IF
         thn(i) = thi
         ni(i) = sw
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni,sy)
      RETURN
      END
      subroutine afmodevn(y,n1,n2,n3,mask,h,vext,sigma)
C
C   Aja-Fernandez Mode Vn (6)
C
      implicit none
      integer n1,n2,n3
      double precision y(n1,n2,n3),sigma(n1,n2,n3),h,vext(2)
      integer mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      double precision m1,m2,z
      ih1=int(h)
      ih2=int(h*vext(1))
      ih3=int(h*vext(2))
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               if(mask(i1,i2,i3).ne.0) THEN
                  ni=0
                  m1=0.d0
                  m2=0.d0
                  DO j1=i1-ih1,i1+ih1
                     if(j1.le.0.or.j1.gt.n1) CYCLE
                     DO j2=i2-ih2,i2+ih2
                        if(j2.le.0.or.j2.gt.n2) CYCLE
                        DO j3=i3-ih3,i3+ih3
                           if(j3.le.0.or.j3.gt.n3) CYCLE
                           if(mask(j1,j2,j3).ne.0) THEN
                              z=y(j1,j2,j3)
                              m1=m1+z
                              m2=m2+z*z
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  m1=m1/ni
                  m2=m2/ni
                  z=m2-m1*m1
                  if(ni.gt.1) THEN
                     sigma(i1,i2,i3)=ni*z/(ni-1)
                  ELSE
                     sigma(i1,i2,i3)=0.d0
                  ENDIF
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine afmodem1(y,n1,n2,n3,mask,h,vext,sigma)
C
C   Aja-Fernandez Mode Vn (6)
C
      implicit none
      integer n1,n2,n3
      double precision y(n1,n2,n3),sigma(n1,n2,n3),h,vext(2)
      integer mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      double precision m1
      ih1=int(h)
      ih2=int(h*vext(1))
      ih3=int(h*vext(2))
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               if(mask(i1,i2,i3).ne.0) THEN
                  ni=0
                  m1=0.d0
                  DO j1=i1-ih1,i1+ih1
                     if(j1.le.0.or.j1.gt.n1) CYCLE
                     DO j2=i2-ih2,i2+ih2
                        if(j2.le.0.or.j2.gt.n2) CYCLE
                        DO j3=i3-ih3,i3+ih3
                           if(j3.le.0.or.j3.gt.n3) CYCLE
                           if(mask(j1,j2,j3).ne.0) THEN
                              m1=m1+y(j1,j2,j3)
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  sigma(i1,i2,i3)=m1/ni
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine afmodem2(y,n1,n2,n3,mask,h,vext,sm)
C
C   Aja-Fernandez Mode Vn (6)
C
      implicit none
      integer n1,n2,n3
      double precision y(n1,n2,n3),sm(n1,n2,n3),h,vext(2)
      integer mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      double precision m2
      ih1=int(h)
      ih2=int(h*vext(1))
      ih3=int(h*vext(2))
      Do i1=1,n1
         Do i2=1,n2
            Do i3=1,n3
               if(mask(i1,i2,i3).ne.0) THEN
                  ni=0
                  m2=0.d0
                  DO j1=i1-ih1,i1+ih1
                     if(j1.le.0.or.j1.gt.n1) CYCLE
                     DO j2=i2-ih2,i2+ih2
                        if(j2.le.0.or.j2.gt.n2) CYCLE
                        DO j3=i3-ih3,i3+ih3
                           if(j3.le.0.or.j3.gt.n3) CYCLE
                           if(mask(j1,j2,j3).ne.0) THEN
                              m2=m2+y(j1,j2,j3)*y(j1,j2,j3)
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  sm(i1,i2,i3)=m2/ni
               ELSE
                  sm(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END

      subroutine paramw3(h,wght,ind,w,n)
C  compute a description of local weights
C  h    - bandwidth
C  wght - vector (length 2) of relative voxel extensions
C  ind  - integer array dim (3,n) containing relative indices in xyz
C  w    - vector of corresponding weights
C  n    - number of positive weights (initial value
C         (2 int(h)+1)*(2 int(h/vext(1))+1)*(2 int(h/vext(2))+1)
      integer n,ind(3,n)
      double precision h,wght(2),w(n)
      integer i,i1,i2,i3,ih1,ih2,ih3
      double precision hsq,z1,z2,z3
      hsq=h*h
      ih1 = int(h)
      ih2 = int(h/wght(1))
      ih3 = int(h/wght(2))
      i=1
      DO i1=-ih1,ih1
         z1=i1*i1
         DO i2=-ih2,ih2
            z2=i2*wght(1)
            z2=z1+z2*z2
            IF(z2.ge.hsq) CYCLE
            DO i3=-ih3,ih3
               z3=i3*wght(2)
               z3=z2+z3*z3
               IF(z3.ge.hsq) CYCLE
               ind(1,i)=i1
               ind(2,i)=i2
               ind(3,i)=i3
               w(i)=1.d0-z3/hsq
               i=i+1
            END DO
         END DO
      END DO
      n=i-1
      RETURN
      END
