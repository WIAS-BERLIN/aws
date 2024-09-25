      subroutine ghfse3i(i4,kstar,k456,ng,
     1                    kappa,vext,h,varred,n,dist)
C
C   compute bandwidth sequence for given kappa and gradients
C   lkfse3i0 computes weighting schemes and corresponding variance reduction
C   k456,ng (input) contain auxiliary statistics(gradients)
C
      implicit none
      integer ng,i4,kstar,n,dist
      double precision k456(3,ng,ng),vext(2),
     1       kappa,h(kstar),varred(kstar)
      integer k,n0,maxn
      double precision hakt,hakt0,vr,ch,chk,vred,v0r
      ch=1.25d0
      hakt=1.d0
C   initialize kappa
C   loop over steps
      call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
      chk=ch*vr
      maxn = 1
      DO k=1,kstar
         call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
C  search for new hakt
         vred=vr/chk
         DO WHILE (vred.lt.1)
            hakt=hakt*1.05
            call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
            vred=vr/chk
         END DO
         DO WHILE (vred.gt.1.01)
            hakt0=hakt
            v0r=vr
            n0=n
            hakt=hakt/1.005
            call lkfse3i0(hakt,kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
            vred=vr/chk
            If (vred.lt.1) THEN
               hakt=hakt0
               vr=v0r
               n=n0
            END IF
         END DO
         h(k) = hakt
         varred(k) = vr
         chk=chk*ch
         maxn=max(maxn,n)
         if(k.eq.kstar) THEN
            call lkfse3i0(h(k),kappa/hakt,i4,k456,ng,
     1                   vext,vr,n,dist)
         END IF
C  number of positive weights for last step in n
      END DO
      n=maxn
      RETURN
      END
      subroutine exppm6(p,ex)
      implicit none
      double precision p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=dsin(p)
      ex(1,3)=0.d0
      ex(2,1)=-dsin(p)
      ex(2,2)=dcos(p)
      ex(2,3)=0.d0
      ex(3,1)=0.d0
      ex(3,2)=0.d0
      ex(3,3)=1.d0
      RETURN
      END
      subroutine exppm5(p,ex)
      implicit none
      double precision p,ex(3,3)
      ex(1,1)=dcos(p)
      ex(1,2)=0.d0
      ex(1,3)=-dsin(p)
      ex(2,1)=0.d0
      ex(2,2)=1.d0
      ex(2,3)=0.d0
      ex(3,1)=dsin(p)
      ex(3,2)=0.d0
      ex(3,3)=dcos(p)
      RETURN
      END
      subroutine exppm4(p,b,ex)
      implicit none
      double precision b,p,ex(3,3)
      double precision D,sqr2,sb,pDs,D2,cpds,spds,spdsh
      sqr2 = dsqrt(2.d0)
      sb = dsin(b)
      D = dsqrt(3.d0-dcos(2.d0*b))
      pDs = p*D/sqr2
      cpds = dcos(pDs)
      spds = dsin(pDs)
      spdsh = dsin(pDs/2.d0)
      D2 = D*D
      ex(1,1)=2.d0*(1.d0+cpDs*sb*sb)/D2
      ex(1,2)=-sqr2/D*sb*spDs
      ex(1,3)=-4.d0/D2*sb*spdsh*spdsh
      ex(2,1)=-ex(1,2)
      ex(2,2)=cpDs
      ex(2,3)=sqr2/D*spDs
      ex(3,1)=ex(1,3)
      ex(3,2)=-ex(2,3)
      ex(3,3)=1.d0-2.d0/D2*(1.d0-cpDs)
      RETURN
      END
      subroutine k456krb(par,b,matm,erg)
C
C  Solve exponential equation for dicrepance parameters
C  compute ||\prod_{i=4}^6 exp(par[i] m_i) - matm||^2
C
      implicit none
      double precision par(3),b,matm(3,3),erg
      integer i1,i2
      double precision s,z,em4(3,3),em5(3,3),em6(3,3),am4(3,3),am5(3,3)
      call exppm4(par(1),b,em4)
      call exppm5(par(2),em5)
      call exppm6(par(3),em6)
      DO i1=1,3
         DO i2=1,3
            am5(i1,i2)=em5(i1,1)*em6(1,i2)+em5(i1,2)*em6(2,i2)+
     1                                     em5(i1,3)*em6(3,i2)
            END DO
         END DO
      DO i1=1,3
         DO i2=1,3
            am4(i1,i2)=em4(i1,1)*am5(1,i2)+em4(i1,2)*am5(2,i2)+
     1                                     em4(i1,3)*am5(3,i2)
         END DO
      END DO
      s=0.d0
      DO i1=1,3
         DO i2=1,3
            z=matm(i1,i2)-am4(i1,i2)
            s=s+z*z
         END DO
      END DO
      erg=s
      RETURN
      END
      subroutine abofg(g,n,bg)
C
C  compute spherical coordinates for gradient vectors
C
      implicit none
      integer n
      double precision g(3,n),bg(2,n)
      integer i
      double precision z,g1,g2,g3,beta,gamma,onemeps
      onemeps=1.d0-1d-10
      DO i=1,n
         g1=g(1,i)
         g2=g(2,i)
         g3=g(3,i)
C standardize to length 1
         z=g1*g1+g2*g2+g3*g3
         z=sqrt(z)
         g1=g1/z
         g2=g2/z
         g3=g3/z
         beta=asin(g1)
         if(abs(g1).lt.onemeps) THEN
            z=g3/cos(beta)
            if(abs(z).lt.onemeps) THEN
               gamma=acos(z)
            ELSE
               gamma=1.570796327d0-sign(1.570796327d0,z)
            END IF
         ELSE
            gamma=0.d0
         END IF
         if(g2.gt.0.d0) gamma = -gamma
         bg(1,i)=beta
         bg(2,i)=gamma
      END DO
      RETURN
      END
      subroutine bgstats(g,n,bg,bghat)
      implicit none
      integer n
      double precision g(3,n),bg(2,n),bghat(2,n,n)
      integer i1,i2
      double precision dgamma,cb1,sb1,cb2,sb2,betah,cbh,z,gammah,cdg
C   first get sperical coordinates in bg
      call abofg(g,n,bg)
      DO i1=1,n
         sb1=sin(bg(1,i1))
         cb1=cos(bg(1,i1))
C    die normalen-vektoren der Gradienten
         DO i2=1,n
            dgamma=bg(2,i1)-bg(2,i2)
            cdg=cos(dgamma)
            if(abs(cdg).gt.1-1d-8) THEN
               bghat(1,i1,i2)=asin(sin(bg(1,i1)-
     1                             sign(1.d0,cdg)*bg(1,i2)))
               bghat(2,i1,i2)=0.d0
               CYCLE
            END IF
            sb2=sin(bg(1,i2))
            cb2=cos(bg(1,i2))
            z=sb1*cb2-cb1*sb2*cdg
            betah=asin(z)
            cbh=cos(betah)
            IF(abs(cbh).gt.1d-8) THEN
               z=cb1*sin(dgamma)/cbh
               z=sign(min(1.d0,abs(z)),z)
               gammah=asin(z)
            ELSE
               if(abs(cb1).gt.1d-6) THEN
C   this should not happen
                  call dblepr("cb1",3,cb1,1)
                  call dblepr("cbh",3,cbh,1)
               END IF
               gammah = dgamma*sign(1d0,cb1*cbh)
            END IF
            bghat(1,i1,i2)=betah
            bghat(2,i1,i2)=gammah
C    die spherischen Koordinaten der Gradientenpaare (Parameter der Rotationsmatix)
         END DO
      END DO
      RETURN
      END
      subroutine lkfse3i(h,kappa,i4,k456,ng,
     1                   vext,ind,wght,n,dist)
      implicit none
      integer ng,n,ind(5,n),i4,dist
      double precision h,kappa,k456(3,ng,ng),vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3,j4
      double precision h2,kap2,x1,x2,x3,k4,k5,k6,z,z1,vd2,vd3
      ih1 = int(max(1.d0,h))
      ih2 = int(max(1.d0,h/vext(1)))
      ih3 = int(max(1.d0,h/vext(2)))
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)*vext(1)
      vd3 = vext(2)*vext(2)
      i = 1
      z = 0.d0
      k5 = 0.d0
      k6 = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         IF(dist.lt.3) THEN
            k5 = k456(2,i4,j4)
            k6 = k456(3,i4,j4)
         END IF
         SELECT CASE (dist)
            CASE (1)
               z = (k4*k4+k5*k5+abs(k6))/kap2
            CASE (2)
               z = (k4*k4+k5*k5+k6*k6)/kap2
            CASE (3)
               z = k4*k4/kap2
            CASE (4)
               z = abs(k4)/kappa
            CASE DEFAULT
               z = abs(k4)/kappa
         END SELECT
         if(dist.le.3) THEN
            if(z.gt.h2) CYCLE
C   last three komponents already to large
            DO j1 = 0,ih1
               x1 = j1
               x1 = z + x1*x1
               if(x1.gt.h2) CYCLE
               DO j2 = -ih2,ih2
                  x2 = j2
                  x2 = x1 + vd2*x2*x2
                  if(x2.gt.h2) CYCLE
                  DO j3 = -ih3,ih3
                     x3 = j3
                     z1 = x2+vd3*x3*x3
                     if(z1.ge.h2) CYCLE
                     if(i.gt.n) THEN
                        call intpr("Exceeded max i",14,i,1)
                        call intpr("for i4",6,i4,1)
                        n = i-1
                        return
                     END IF
                     wght(i)= (1.d0-z1/h2)
                     ind(1,i) = j1
                     ind(2,i) = j2
                     ind(3,i) = j3
                     ind(4,i) = i4
                     ind(5,i) = j4
                     i = i+1
                  END DO
                  call rchkusr()
               END DO
            END DO
         ELSE
C dist=4
            if(z.gt.h) CYCLE
            DO j1 = 0,ih1
               x1 = j1
               x1 = x1*x1
               DO j2 = -ih2,ih2
                  x2 = j2
                  x2 = x1 + vd2*x2*x2
                  DO j3 = -ih3,ih3
                     x3 = j3
                     z1 = x2+vd3*x3*x3
                     z1=z+sqrt(z1)
                     if(z1.ge.h) CYCLE
                     if(i.gt.n) THEN
                        call intpr("Exceeded max i",14,i,1)
                        call intpr("for i4",6,i4,1)
                        n = i-1
                        return
                     END IF
                     wght(i)= (1.d0-z1*z1/h2)
                     ind(1,i) = j1
                     ind(2,i) = j2
                     ind(3,i) = j3
                     ind(4,i) = i4
                     ind(5,i) = j4
                     i = i+1
                  END DO
                  call rchkusr()
               END DO
            END DO
         ENDIF
      END DO
      n = i-1
      RETURN
      END
      subroutine lkfse3i0(h,kappa,i4,k456,ng,vext,vred,n,dist)
      implicit none
      integer ng,n,i4,dist
      double precision h,kappa,k456(3,ng,ng),vext(2),vred
      integer ih1,ih2,ih3,j1,j2,j3,j4
      double precision x1,x2,x3,k4,k5,k6,z,z1,
     1       sw,sw2,wght,anz,h2,kap2,vd2,vd3
      ih1 = int(max(1.d0,h))
      ih2 = int(max(1.d0,h/vext(1)))
      ih3 = int(max(1.d0,h/vext(2)))
      sw=0.d0
      sw2=0.d0
      h2 = h*h
      kap2 = kappa*kappa
      vd2 = vext(1)*vext(1)
      vd3 = vext(2)*vext(2)
      n = 0
      z = 0.d0
      k5 = 0.d0
      k6 = 0.d0
C  just to prevent compiler warnings
      DO j4 = 1,ng
         k4 = k456(1,i4,j4)
         IF(dist.lt.3) THEN
            k5 = k456(2,i4,j4)
            k6 = k456(3,i4,j4)
         END IF
         SELECT CASE (dist)
            CASE (1)
               z = (k4*k4+k5*k5+abs(k6))/kap2
            CASE (2)
               z = (k4*k4+k5*k5+k6*k6)/kap2
            CASE (3)
               z = k4*k4/kap2
            CASE (4)
               z = abs(k4)/kappa
            CASE DEFAULT
               z = abs(k4)/kappa
         END SELECT
         if(dist.le.3) THEN
            if(z.gt.h2) CYCLE
C   last three komponents already to large
            DO j1 = 0,ih1
C   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice
               if(j1.eq.0) THEN
                  anz=1.d0
               ELSE
                  anz=2.d0
               ENDIF
               x1 = j1
               x1 = z + x1*x1
               if(x1.gt.h2) CYCLE
               DO j2 = -ih2,ih2
                  x2 = j2
                  x2 = x1 + vd2*x2*x2
                  if(x2.gt.h2) CYCLE
                  DO j3 = -ih3,ih3
                     x3 = j3
                     z1 = x2+vd3*x3*x3
C   corrected from vd2 to vd3 J.P. 29.8.2013
                     if(z1.gt.h2) CYCLE
                     wght= (1.d0-z1/h2)
                     sw=sw+anz*wght
                     wght=wght*wght
                     sw2=sw2+anz*wght
                     n=n+1
                  END DO
                  call rchkusr()
               END DO
            END DO
         ELSE
            if(z.gt.h) CYCLE
C   last three komponents already to large
            DO j1 = 0,ih1
C   if j1>0  (-j1,-j2,-j3) gets the same weight, so count it twice
               if(j1.eq.0) THEN
                  anz=1.d0
               ELSE
                  anz=2.d0
               ENDIF
               x1 = j1
               x1 = x1*x1
               DO j2 = -ih2,ih2
                  x2 = j2
                  x2 = x1 + vd2*x2*x2
                  if(x2.gt.h2) CYCLE
                  DO j3 = -ih3,ih3
                     x3 = j3
                     z1 = x2+vd3*x3*x3
                     z1=z+sqrt(z1)
                     if(z1.gt.h) CYCLE
                     wght= (1.d0-z1*z1/h2)
                     sw=sw+anz*wght
                     wght=wght*wght
                     sw2=sw2+anz*wght
                     n=n+1
                  END DO
                  call rchkusr()
               END DO
            END DO
         END IF
      END DO
      vred = sw*sw/sw2
      RETURN
      END
      subroutine lkfulse3(h,kappa,k456,ng,vext,ind,wght,n,dist)
      implicit none
      integer ng,n,ind(5,n),dist
      double precision h(ng),kappa(ng),k456(3,ng,ng),vext(2),wght(n)
      integer ns,ni,i
      ns = 0
      DO i = 1,ng
         ni = n-ns
         call lkfse3i(h(i),kappa(i),i,k456,ng,
     1                vext,ind(1,ns+1),wght(ns+1),ni,dist)
         ns = ns+ni
      END DO
      n = ns
      RETURN
      END
      subroutine adsmse3p(y,th,ni,pos,nv,n1,n2,n3,ngrad,lambda,ncoils,
     1                    ncores,ind,w,n,thn,ldf,sw,swy,model)
C   model=1 takes noncentral Chi-sq values in y
C   model=0 takes noncentral Chi values in y
C
C   perform adaptive smoothing on SE(3)
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
c   model=0  Chi^2-based KL-distance, y ~ Chi, th on same scale, smooth y
c   model=1  Chi^2-based KL-distance, y ~ Chi^2, th on same scale, smooth y
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2
C!$      use omp_lib
      implicit none
      integer n1,n2,n3,ngrad,n,ind(5,n),ncoils,model,ncores
      integer pos(*),nv
      double precision y(nv,ngrad),thn(nv,ngrad),
     1     ni(nv,ngrad),ldf(nv,ngrad),th(nv,ngrad)
      double precision lambda,w(n),sw(ngrad,ncores),swy(ngrad,ncores),
     2       lgfi,dgfi,fici,df
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,
     1       jind,iindp,jindp,n12
      double precision z,thi,nii,thj,ldfi,ldfj,yj
      double precision kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      nii=1.d0
      thrednr = 1
      n12=n1*n2
C just to prevent a compiler warning
C  precompute values of lgamma(corrected df/2) in each voxel
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ngrad,ncores,pos,nv,y,thn,ni,ldf,th,w,sw,swy,
C$OMP& ind,ncoils,model,df,n,lambda,n12)
C$OMP& FIRSTPRIVATE(nii)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,z,thi,thj,
C$OMP& ldfi,ldfj,dgfi,fici,lgfi,yj,jind,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp = pos(iind)
         if(iindp.eq.0) CYCLE
         if(model.lt.2) THEN
            i1=mod(iind,n1)
            if(i1.eq.0) i1=n1
            i2=mod((iind-i1)/n1+1,n2)
            if(i2.eq.0) i2=n2
            i3=(iind-i1-(i2-1)*n1)/n1/n2+1
C  not needed for Gauss approximation
            DO i4=1,ngrad
               thi=th(iindp,i4)
               call lgstats(thi,df,model,ldfi)
               ldf(iindp,i4)=ldfi
            END DO
         END IF
      END DO
C$OMP END DO
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
C!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         iindp = pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               thi = th(iindp,i4)
               IF(model.lt.2) THEN
                  ldfi=ldf(iindp,i4)
                  call ncstats0(thi,ldfi,df,model,lgfi,dgfi,fici)
               END IF
               nii = ni(iindp,i4)/lambda
           END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
C adaptation
            if(lambda.lt.1d10) THEN
               thj=th(jindp,j4)
               if(model.ge.2) THEN
                  z=(thi-thj)
                  z=nii*z*z
               ELSE
                  ldfj=ldf(jindp,j4)
                  z=nii*kldisnc1(lgfi,dgfi,fici,thj,ldfj,df,model)
               END IF
C  do not adapt on the sphere !!!
            ELSE
            z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            yj=y(jindp,j4)
            if(model.eq.2) yj=yj*yj
            swy(i4,thrednr)=swy(i4,thrednr)+z*yj
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               thi = th(iindp,i4)
               IF(model.lt.2) THEN
                  ldfi=ldf(iindp,i4)
                  call ncstats0(thi,ldfi,df,model,lgfi,dgfi,fici)
               END IF
               nii = ni(iindp,i4)/lambda
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               thj=th(jindp,j4)
               if(model.ge.2) THEN
                  z=(thi-thj)
                  z=nii*z*z
               ELSE
                  ldfj=ldf(jindp,j4)
                  z=nii*kldisnc1(lgfi,dgfi,fici,thj,ldfj,df,model)
               END IF
C  do not adapt on the sphere !!!
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            yj=y(jindp,j4)
            if(model.eq.2) yj=yj*yj
            swy(i4,thrednr)=swy(i4,thrednr)+z*yj
         END DO
         DO i4=1,ngrad
            z=swy(i4,thrednr)/sw(i4,thrednr)
            if(model.eq.2) z=sqrt(z)
C            thn(i1,i2,i3,i4) = max(z,minlev)
            thn(iindp,i4) = z
            ni(iindp,i4) = sw(i4,thrednr)
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
      subroutine adsmse3m(y,th,ni,sthi,pos,nv,ns,n1,n2,n3,ngrad,lambda,
     1             ncores,ind,w,n,thn,sw,swy,si,thi)
C
C  Multi-shell version (differs in dimension of th
C  KL-distance based on all spheres and Gauss-approximation only
C
C   perform adaptive smoothing on SE(3)
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2
C!$      use omp_lib
      implicit none
      integer ns,n1,n2,n3,ngrad,n,ind(5,n),ncores
      integer pos(*),nv
      double precision y(nv,ngrad),thn(nv,ngrad),
     1 ni(nv,ngrad),th(ns,nv,ngrad),sthi(ns,nv,ngrad)
      double precision lambda,w(n),sw(ngrad,ncores),swy(ngrad,ncores)
      integer iind,jind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,n12,
     1        iindp,jindp
      double precision sz,z,nii,si(ns,ncores),thi(ns,ncores)
      double precision kldisnc1
      external kldisnc1
      nii=1.d0
      thrednr = 1
      n12=n1*n2
C just to prevent a compiler warning
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,ngrad,ncores,pos,nv,y,thn,ni,th,sthi,w,sw,swy,
C$OMP& ind,n,lambda,ns,thi,si,n12)
C$OMP& FIRSTPRIVATE(nii)
C$OMP& PRIVATE(iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,z,sz,
C$OMP& jind,iindp,jindp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp = pos(iind)
         if(iindp.eq.0) CYCLE
C!$       thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,iindp,i4)
C   thast the approximated standard deviation
                  si(k,thrednr) = sthi(k,iindp,i4)
               END DO
               nii = ni(iindp,i4)/lambda
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
C adaptation
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,jindp,j4))/si(k,thrednr)
                  sz=sz+z*z
               END DO
               z=nii*sz
C  do not adapt on the sphere !!!
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(jindp,j4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               DO k=1,ns
                  thi(k,thrednr) = th(k,iindp,i4)
                  si(k,thrednr) = sthi(k,iindp,i4)
               END DO
               nii = ni(iindp,i4)/lambda
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,jindp,j4))/si(k,thrednr)
                  sz=sz+z*z
               END DO
               z=nii*sz
C  do not adapt on the sphere !!!
            ELSE
               z=0.d0
            END IF
            if(z.ge.1.d0) CYCLE
            z=w(i)*min(1.d0,2.d0-2.d0*z)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(jindp,j4)
         END DO
         DO i4=1,ngrad
            thn(iindp,i4) = swy(i4,thrednr)/sw(i4,thrednr)
            ni(iindp,i4) = sw(i4,thrednr)
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
      subroutine adsmse3s(y,y0,th,ni,th0,ni0,fsi2,fsi02,pos,nv,ns,n1,
     1                n2,n3,ngrad,lambda,ws0,ind,w,n,ind0,w0,
     2                n0,thn,nin,th0n,ni0n,sw,swy,thi,nii,fsi2i)
C
C  Multi-shell version (differs in dimension of th
C  KL-distance based on all spheres and Gauss-approximation only
C  see ~polzehl/latex/1211_simplemetric/Rcode/Figuregaussapprox.r
C  for approximative KL-distance between non-central chi distributions
C
C   perform adaptive smoothing on SE(3) multishell including s0
C   y  -  si images
C   y0 -  mean s0 image
C   th -  estimated/interpolated \E si on all shells (including 0 shell)
C   ni -  corresponding sum of weights
C   th0 -  estimated/interpolated \E s0 and mean_g(\E si) on all other shells
C   ni0 -  corresponding sum of weights
C   pos - index of voxel in head mask
C   ns   - number of shells (including 0 shell)
C   n1,n2,n3,ngrad - dimensions, number of gradients (bv!=0)
C   lambda - skale parameter
C   ws0  - relative weight for information from s0 images (should be in [0,1])
C   ncoils - df/2 of \chi distributions
C   minlev - expectation of central chi distr. (needed in variance estimates)
C   ind    - index vectors for si weighting schemes
C   w    - corresponding weights
C   n    - number of si weights
C   ind0    - index vectors for s0 weighting schemes
C   w0    - corresponding weights
C   n0    - number of s0 weights
C   thn   - new estimates of \E si
C   thn0   - new estimates of \E s0
C   nin   - new sum of weight for si
C   ni0n   - new sum of weight for s0
C...sw,swy,si,thi,nii - working areas
C   ind(.,i) contains coordinate information corresponding to positive
C   location weights in w(i) for si images
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
C!$      use omp_lib
      implicit none
      integer nv,ns,n1,n2,n3,ngrad,n,n0,ind(5,n),ind0(3,n0),pos(*)
      double precision y(*),y0(nv),th(ns,*),ni(ns,*),th0(ns,nv),
     1  ni0(ns,nv),fsi2(ns,*),fsi02(ns,nv),thn(*),th0n(nv),nin(*),
     2  ni0n(nv)
C  * refers to n1*n2*n3*ngrad for y,th,ni,thn,fsi2,nin and to
C              n1*n2*n3 for y0,th0,ni0,th0n,ni0n,fsi02,mask
      double precision w(n),w0(n0),lambda,thi(ns,*),ws0,fsi2i(ns,*),
     1       sw(ngrad,*),swy(ngrad,*),nii(ns,*)
C  * refers to ns*ncores in thi, fsi2i, nii and to
C              ngrad*ncores in sw and swy
      integer iind,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,jind,iind4,
     1        jind4,n12,iindp,jindp
      double precision sz,z,sw0,swy0
      thrednr = 1
      n12 = n1*n2
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ns,n1,n2,n3,ngrad,n,n0,ind,ind0,y,y0,nv,
C$OMP&       th,ni,th0,ni0,w,w0,thn,th0n,nin,ni0n,thi,sw,swy,nii,
C$OMP&       lambda,pos,ws0,fsi2,fsi02,fsi2i,n12)
C$OMP& PRIVATE(iind,iindp,i,i1,i2,i3,i4,j1,j2,j3,j4,thrednr,k,sz,z,
C$OMP&       sw0,swy0,jind,jindp,iind4,jind4)
C$OMP DO SCHEDULE(GUIDED)
C  First si - images
      DO iind=1,n1*n2*n3
         iindp = pos(iind)
         if(iindp.eq.0) CYCLE
C!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         DO i4=1,ngrad
            sw(i4,thrednr)=0.d0
            swy(i4,thrednr)=0.d0
         END DO
         i4=0
         DO i=1,n
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               iind4 = iindp+(i4-1)*nv
               DO k=1,ns
                  fsi2i(k,thrednr)=fsi2(k,iind4)
                  thi(k,thrednr) = th(k,iind4)
                  nii(k,thrednr) = ni(k,iind4)/lambda
               END DO
               nii(1,thrednr)=ws0*nii(1,thrednr)
            END IF
            j1=i1+ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
            jind4=jindp+(j4-1)*nv
C adaptation
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,jind4))
                  sz=sz+nii(k,thrednr)*z*z/
     1                        (fsi2(k,jind4)+fsi2i(k,thrednr))
                  if(sz.ge.1.d0) EXIT
               END DO
C  do not adapt on the sphere !!!
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(jind4)
         END DO
C  now opposite directions
         DO i=1,n
            if(ind(1,i).eq.0) CYCLE
            if(ind(4,i).ne.i4) THEN
C   by construction ind(4,.) should have same values consequtively
               i4 = ind(4,i)
               iind4 = iindp+(i4-1)*nv
               DO k=1,ns
                  fsi2i(k,thrednr)=fsi2(k,iind4)
                  thi(k,thrednr) = th(k,iind4)
                  nii(k,thrednr) = ni(k,iind4)/lambda
               END DO
               nii(1,thrednr)=ws0*nii(1,thrednr)
C  first component corresponds to S0 image, ws0 is used to downweight its influence
C  when smoothing diffusion weighted data
            END IF
C
C   handle case j1-i1 < 0 which is not contained in ind
C   using axial symmetry
C
            j1=i1-ind(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            j4=ind(5,i)
            jind4=jindp+(j4-1)*nv
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th(k,jind4))
                  sz=sz+nii(k,thrednr)*z*z/
     1                        (fsi2(k,jind4)+fsi2i(k,thrednr))
                  if(sz.ge.1.d0) EXIT
               END DO
C  do not adapt on the sphere !!!
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            sw(i4,thrednr)=sw(i4,thrednr)+z
            swy(i4,thrednr)=swy(i4,thrednr)+z*y(jind4)
         END DO
         DO i4=1,ngrad
            iind4 = iindp+(i4-1)*nv
            thn(iind4) = swy(i4,thrednr)/sw(i4,thrednr)
            nin(iind4) = sw(i4,thrednr)
         END DO
C    now the s0 image in iind
         sw0=0.d0
         swy0=0.d0
         DO k=1,ns
            thi(k,thrednr) = th0(k,iindp)
            nii(k,thrednr) = ni0(k,iindp)/lambda
            fsi2i(k,thrednr)=fsi02(k,iindp)
         END DO
         DO i=1,n0
            j1=i1+ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2+ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3+ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
C adaptation
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th0(k,jindp))
                  sz=sz+nii(k,thrednr)*z*z/
     1                        (fsi02(k,jindp)+fsi2i(k,thrednr))
                  if(sz.ge.1.d0) EXIT
               END DO
C  do not adapt on the sphere !!!
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(jindp)
         END DO
C  now opposite directions
         DO i=1,n0
            if(ind0(1,i).eq.0) CYCLE
C
C   handle case j1-i1 < 0 which is not contained in ind
C   using axial symmetry
C
            j1=i1-ind0(1,i)
            if(j1.le.0.or.j1.gt.n1) CYCLE
            j2=i2-ind0(2,i)
            if(j2.le.0.or.j2.gt.n2) CYCLE
            j3=i3-ind0(3,i)
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            jindp=pos(jind)
            if(lambda.lt.1d10) THEN
               sz=0.d0
               DO k=1,ns
                  z=(thi(k,thrednr)-th0(k,jindp))
                  sz=sz+nii(k,thrednr)*z*z/
     1                        (fsi02(k,jindp)+fsi2i(k,thrednr))
                  if(sz.ge.1.d0) EXIT
               END DO
C  do not adapt on the sphere !!!
            ELSE
               sz=0.d0
            END IF
            if(sz.ge.1.d0) CYCLE
            z=w0(i)
            if(sz.gt.0.5d0) z=z*(2.d0-2.d0*sz)
            sw0=sw0+z
            swy0=swy0+z*y0(jindp)
         END DO
         th0n(iindp) = swy0/sw0
         ni0n(iindp) = sw0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,nin,th0n,ni0n)
      RETURN
      END
      subroutine asmse30p(y,th,ni,pos,nv,n1,n2,n3,lambda,ncoils,ind,w,
     1                    n,starts,nstarts,thn,ldf,swi,model)
C   model=1 takes noncentral Chi-sq values in y0
C   model=0 takes noncentral Chi values in y0
C   perform adaptive smoothing on R^3
C   ind(.,i) contains coordinate indormation corresponding to positive
C   location weights in w(i)
C   ind(.,i)[1:5] are j1-i1,j2-i2,j3-i3, i4 and j4 respectively
C
c   model=0  Chi^2-based KL-distance, y ~ Chi, th on same scale, smooth y
c   model=1  Chi^2-based KL-distance, y ~ Chi^2, th on same scale, smooth y
c   model=2  Gauss-based KL-distance, y ~ Chi, th on same scale, smooth y^2

      implicit none
      integer n1,n2,n3,n,ind(5,n),starts(*),nstarts,ncoils,model
      integer pos(*),nv
      double precision y(nv),th(nv),ni(nv),thn(nv),
     1       lambda,w(n),sw0,swy0,swi(nstarts),ldf(nv)
      integer i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,nn,iindp,jindp,jind,n12
      double precision z,ldfi,lgfi,dgfi,fici,sc,df,yj,maxswi
      double precision kldisnc1
      external kldisnc1
      df=2.d0*ncoils
      nn=n1*n2*n3
      n12=n1*n2
C
C  compute location weights in swi(i0) as sum of location weights used for same
C  relative voxel in R3 and coinciding gradients (sum over gradients)
C
      maxswi = 0.d0
      DO i0=1,nstarts
         swi(i0)=0.d0
         DO i=starts(i0)+1,starts(i0+1)
            swi(i0)=swi(i0)+w(i)
         END DO
         swi(i0)=swi(i0)/(starts(i0+1)-starts(i0))
      END DO
      maxswi=1.d0
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,nn,ncoils,pos,nv,y,th,ni,thn,ind,ldf,n,nstarts,
C$OMP& starts,w,swi,df,n12)
C$OMP& FIRSTPRIVATE(model,lambda,maxswi)
C$OMP& PRIVATE(sw0,swy0,i,i0,i1,i2,i3,j1,j2,j3,l1,l2,l3,z,lgfi,dgfi,
C$OMP& ldfi,fici,sc,yj,iindp,jindp,jind)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
         iindp = pos(i)
         if(iindp.eq.0) CYCLE
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         thn(iindp) = th(iindp)
         if(model.lt.2) THEN
            call lgstats(th(iindp),df,model,ldfi)
            ldf(i) = ldfi
         END IF
      END DO
C$OMP END DO
C$OMP BARRIER
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nn
         iindp = pos(i)
         if(iindp.eq.0) CYCLE
         i1=mod(i,n1)
         if(i1.eq.0) i1=n1
         i2=mod((i-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(i-i1-(i2-1)*n1)/n1/n2+1
         sw0=0.d0
         swy0=0.d0
         sc = lambda/ni(i)
         IF(model.lt.2) THEN
            call ncstats0(th(iindp),ldf(iindp),df,model,lgfi,dgfi,fici)
         END IF
         DO i0=1,nstarts
            l1=ind(1,1+starts(i0))
            j1=i1+l1
            if(j1.le.0.or.j1.gt.n1) CYCLE
            l2=ind(2,1+starts(i0))
            j2=i2+l2
            if(j2.le.0.or.j2.gt.n2) CYCLE
            l3=ind(3,1+starts(i0))
            j3=i3+l3
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            if(model.ge.2) THEN
               z=th(iindp)-th(jindp)
               z=z*z
            ELSE
               z=kldisnc1(lgfi,dgfi,fici,th(jindp),
     1                        ldf(jindp),df,model)
            END IF
            if(z.ge.sc) CYCLE
            if(z.lt.sc) THEN
               z=swi(i0)*min(1.d0,2.d0-2.d0*z/sc)
               sw0=sw0+z
               yj=y(jindp)
               if(model.eq.2) yj=yj*yj
               swy0=swy0+z*yj
            END IF
         END DO
         DO i0=1,nstarts
            l1=ind(1,1+starts(i0))
            if(l1.eq.0) CYCLE
C  this part for negative l1 only (opposite directions)
            j1=i1-l1
            if(j1.le.0.or.j1.gt.n1) CYCLE
            l2=ind(2,1+starts(i0))
            j2=i2-l2
            if(j2.le.0.or.j2.gt.n2) CYCLE
            l3=ind(3,1+starts(i0))
            j3=i3-l3
            if(j3.le.0.or.j3.gt.n3) CYCLE
            jind=j1+(j2-1)*n1+(j3-1)*n12
            jindp=pos(jind)
            if(jindp.eq.0) CYCLE
            if(model.ge.2) THEN
               z=th(iindp)-th(jindp)
               z=z*z
            ELSE
               z=kldisnc1(lgfi,dgfi,fici,th(jindp),
     1         ldf(jindp),df,model)
            END IF
            if(z.ge.sc) CYCLE
            if(z.lt.sc) THEN
               z=swi(i0)*min(1.d0,2.d0-2.d0*z/sc)
               sw0=sw0+z
               yj=y(jindp)
               if(model.eq.2) yj=yj*yj
               swy0=swy0+z*yj
            END IF
         END DO
         z = swy0/sw0
         if(model.eq.2) z=sqrt(z)
         thn(iindp) = z
         ni(iindp) = sw0/maxswi
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,ni)
      RETURN
      END
C
C   Kullback-leibler distance for noncentral-Chi-distributions
C   with parameters thi and thj and 2*nc degrees of freedom
C
C
C   variant with precomputed quantities
C
      subroutine lgstats(thi,df,model,lgfi)
      implicit none
      integer model
      double precision thi,lgfi,df
      double precision mu2i,z1,z2,fi
      double precision lgammaf,digammaf
      external lgammaf, digammaf
      mu2i = thi
      if(model.eq.0) mu2i = mu2i*mu2i
      mu2i = max(0.d0,mu2i-df)
      z1 = df + mu2i
      z2 = z1 + mu2i
      fi = z1*z1/z2
      lgfi = lgammaf(fi/2.d0)
      RETURN
      END
      subroutine ncstats0(thi,lgfi0,df,model,lgfi,dgfi,fici)
      implicit none
      integer model
      double precision thi,lgfi0,lgfi,dgfi,fici,df
      double precision mu2i,z1,z2,dlci,fi,ci
      double precision digammaf
      external digammaf
      mu2i = thi
      if(model.eq.0) mu2i = mu2i*mu2i
      mu2i = max(0.d0,mu2i-df)
      z1 = df + mu2i
      z2 = z1 + mu2i
      ci = z2/z1
      fi = z1/ci
      fici = fi*ci
      dlci = dlog(ci)
      dgfi = digammaf(0.5d0*fi)+dlci
      lgfi = lgfi0+0.5d0*(fi*dlci+fi-fi*dgfi)
      RETURN
      END
      double precision function kldisnc1(lgfi,dgfi,fici,thj,lgfj,df,
     1                          model)
C    for smoothing noncentral Chi values
C    thi,thj  current estimates
C    df= 2 * number of coils
C    model = 0   smoothing of chi values
C    model = 1   smoothing of chi^2 values
      implicit none
      integer model
      double precision lgfi,dgfi,fici,thj,df,lgfj
      double precision mu2j,fj,cj,z1,z2
C  use Chi^2 instead of Chi for KL distance
      mu2j = thj
      if(model.eq.0) mu2j = mu2j*mu2j
      mu2j = max(0.d0,mu2j-df)
C  Approximation by Patnaik (1949)
      z1 = df + mu2j
      z2 = z1 + mu2j
      cj = z2/z1
      fj = z1/cj
      kldisnc1 = lgfj-lgfi+0.5d0*(fj*dlog(cj)+
     1                 fici/cj-fj*dgfi)
      RETURN
      END
      subroutine lkfuls0(h,vext,ind,wght,n)
      implicit none
      integer n,ind(3,n)
      double precision h,vext(2),wght(n)
      integer ih1,ih2,ih3,i,j1,j2,j3
      double precision h2,x1,x2,x3,z,z1,vd2,vd3
      vd2 = vext(1)
      vd3 = vext(2)
      ih1 = int(max(1.d0,5.0d0*h))
      ih2 = int(max(1.d0,5.0d0*h/vd2))
      ih3 = int(max(1.d0,5.0d0*h/vd3))
      h2 = h*h
      i = 1
      z = 0.d0
C  just to prevent compiler warnings
      DO j1 = 0,ih1
         x1 = j1
         DO j2 = -ih2,ih2
            x2 = vd2*j2
            DO j3 = -ih3,ih3
               x3 = vd3*j3
               z1 = z+x1*x1+x2*x2+x3*x3
               if(z1.ge.h2) CYCLE
               if(i.gt.n) THEN
                  call intpr("Exceeded max i",14,i,1)
                  n = i-1
                  return
               END IF
               wght(i)= (1.d0-z1/h2)
               ind(1,i) = j1
               ind(2,i) = j2
               ind(3,i) = j3
               i = i+1
            END DO
         END DO
      END DO
      n = i-1
      RETURN
      END
      subroutine ipolsp1(theta,th0,ni,ni0,n,ng,gind,gw,nbv,nbvp1,
     1                    msth,msni)
C   interpolate values of theta on spheres where it was not observed
C   this only gets voxel within mask
      implicit none
      integer n,ng,nbv,nbvp1,gind(3,nbv,ng)
      double precision theta(n,ng),th0(n),ni(n,ng),ni0(n),gw(3,nbv,ng),
     1       msth(nbvp1,n,ng),msni(nbvp1,n,ng)
      integer i,j,k,i1,i2,i3,ip1
      double precision w1,w2,w3
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,k,w1,w2,w3,i1,i2,i3,ip1)
C$OMP DO SCHEDULE(GUIDED)
      DO j=1,ng
         DO k=1,n
            msth(1,k,j) = th0(k)
            msni(1,k,j) = ni0(k)
            DO i=1,nbv
               ip1 = i+1
               i1 = gind(1,i,j)
               IF(i1.eq.j) THEN
C  same shell just copy
                  msth(ip1,k,j) = theta(k,j)
                  msni(ip1,k,j) = ni(k,j)
               ELSE
C  different shell need to interpolate
                  i2 = gind(2,i,j)
                  i3 = gind(3,i,j)
                  w1 = gw(1,i,j)
                  w2 = gw(2,i,j)
                  w3 = gw(3,i,j)
                  msth(ip1,k,j) = theta(k,i1)*w1+theta(k,i2)*w2+
     1                            theta(k,i3)*w3
                  msni(ip1,k,j) = 1.d0/
     1                           (w1/ni(k,i1)+w2/ni(k,i2)+w3/ni(k,i3))
               END IF
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(msth,msni)
      RETURN
      END
      subroutine getmsth0(theta,n,lindi,msth)
      implicit none
      integer n,lindi
      double precision theta(n,lindi),msth(n)
      double precision s
      integer i,j
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,s,j)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         s=0.d0
         DO j=1,lindi
            s=s+theta(i,j)
         END DO
         msth(i)=s/lindi
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
      return
      end
      subroutine getmsni0(ni,n,lindi,msni)
      implicit none
      integer n,lindi
      double precision ni(n,lindi),msni(n)
      double precision s
      integer i,j
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,s,j)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         s=0.d0
         DO j=1,lindi
            s=s+1.d0/ni(i,j)
         END DO
         msni(i)=lindi/s
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
      return
      end
