c     --------------------------------------------------------------
c     SCATTERING SUBROUTINES FOR INTERPOED VERSION
c     --------------------------------------------------------------
c     
      
      subroutine findmz (emax,vmin,nsf,zmin,zmax,mz)
      implicit double precision (a-h,o-z)
      
c     ----------------------------------------------------------------- 
c     This subroutine determines the number of translational grid
c     points mz that are to be used in the scattering calculation.
c     
c     The number of the Lobatto functions per de Broglie wave length
c     is depend on the significant figures of accuracy.  Here the mz 
c     was determined by the critior from Fig.2 of DEM,  "Lobatto shape 
c     functions", in "Numerical Grid Methods and Their Application to 
c     Schrodinger's Equation" (C. Cerjan, ed., NATO/ASI, Kluwer, 
c     Dordrecht, 1993)).
c     
c     emax     maximum energy for z component
c     vmin     estimated surface potential depth (latteral averaged 
c     potential)
c     mz       the total number of Labotto functions
c     
c     ----------------------------------------------------------------- 
      
      common /const/ rmlmda
      DATA   Pi /3.141592653589793d0/


      em = emax+dabs(vmin)
      pz = sqrt(rmlmda*em)
      wz = 2.0d0*Pi/pz
      mz = (2.75d0+0.125d0*nsf)*(zmax-zmin)/wz+1
      print *, 'mz = ', mz, ', setting mz = 550'
      mz = 550
      return
      end subroutine findmz
      

      subroutine basis(d,ix,iy,n,n00,dmax,imax)
      implicit double precision (a-h,o-z)
c     
c     calculate reciprocal lattice
c     calculate d (z-component of energy of outgoing wave) for 
c     each channel
c     
c     gax, gbx  = the unit vector of reciprocal lattice along symmetry 
c     direction
c     
c     gay, gby  = y components of unit vector of reciprocal lattice
c     along symmetry direction 2
c     
c     a1     length of real space lattice vector along axis
c     a2     x coordinate of other real space lattice vector
c     b2     y coordinate of other real space lattice vector
c     2
c                                    y|   _____a1_____ 
c     note that the a1,               |  /|          /
c     and b2 are the                  | / |b2       /
c     paramenters of the              |/__|________/____ x
c     unit cell of substrate           a2    
c     
c
c     For each scattered channel:
c     ered = square of incident wavevector for the channel
c     eint = squre of surface component of the wavevector for the 
c     channel
c     
c     -d(i) = square of the z component of the wavevector for the channel
c     energy conservation: ered = eint + (-d(i))
c     if d(i) < 0, channel open, possible diffraction spot
c     if d(i) > 0, channel closed, no spot
c     
      include 'multiscat.inc'
      dimension d(nmax), ix(nmax), iy(nmax)
      
      common /cells/ a1,a2,b2,ei,theta,phi,a0,gax,gay,gbx,gby
      common /const/ rmlmda ! = 2m/h^2
      DATA   Pi /3.141592653589793d0/
           
      ax=a1
      ay=0
      bx=a2
      by=b2
      
      Auc=dabs(ax*by)
      RecUnit=2*Pi/Auc
      gax =  by*RecUnit
      gay = -bx*RecUnit
      gbx = -ay*RecUnit
      gby =  ax*RecUnit

      ered   = rmlmda*ei ! ered is just k_i^2
      thetad = theta*pi/180.0d0
      phid   = phi*pi/180.0d0 

      pkx = sqrt(ered)*sin(thetad)*cos(phid)
      pky = sqrt(ered)*sin(thetad)*sin(phid)
      
      n=0
      do i1 = -imax,imax
         do i2 = -imax,imax
            gx = gax*i1 + gbx*i2
            gy = gay*i1 + gby*i2
            eint = (pkx+gx)**2 + (pky+gy)**2
            di = eint-ered
            if (di.lt.dmax) then 
               n=n+1
               if (n.le.nmax) then
                  ix(n)=i1
                  iy(n)=i2
                  d(n)=di
                  if ((i1.eq.0) .and. (i2.eq.0)) n00=n
               else
                  stop 'ERROR: n too big! (basis)'
            end if
            end if
         end do
      end do
      
      
      return
      end subroutine basis
         
      
      subroutine tshape (a,b,m,w,x,t)
      implicit double precision (a-h,o-z)
      include 'multiscat.inc'
c     
c     ----------------------------------------------------------------- 
c     This subroutine calculates the kinetic energy matrix, T, 
c     in a normalised Lobatto shape function basis.
c
c     Formula for this are taken from: 
c     "QUANTUM SCATTERING VIA THE LOG DERIVATIVE VERSION
c     OF THE KOHN VARIATIONAL PRINCIPLE" 
c     D. E. Manolopoulos and R. E. Wyatt, 
c     Chem. Phys. Lett., 1988, 152,23
c
c     In that paper Lobatto shape functions (Lsf) are defined
c     -----------------------------------------------------------------  
c     
      dimension w(m),x(m),t(m,m)

      dimension ww(mmax+1),xx(mmax+1),tt(mmax+1,mmax+1)
      if (m .gt. mmax) stop 'tshape 1'

c     I think, that this is needed for the sum defined in Lsf to work
      n = m+1
c     Get points and weights for n point Lobatto quadrature in (a,b)
      call lobatto (a,b,n,ww,xx)
      
c     No idea why it's done
      do 1 i = 1,n
         ww(i) = sqrt(ww(i))
 1    continue
 
      do 4 i = 1,n
         ff = 0.0d0
         do 3 j = 1,n
c			gg of i = j is trivially = 0, so no need for loops       
            if (j .eq. i) go to 3

c           gg will be value of derivative of j-th Lsf at
c           i-th root, which is: i-th Lsf evaluated at j-th
c           root divided by ( i-th root minus j-th root )
            gg = 1.0d0/(xx(i)-xx(j))
            ff = ff+gg
            
            do 2 k = 1,n
            
c              This loop multiplies gg defined above by j-th Lsf
c              evaluated at i-th root, which is itself a Lagrangian interpolation
               if (k.eq.j .or. k.eq.i) go to 2
               gg = gg*(xx(j)-xx(k))/(xx(i)-xx(k))
               
 2          continue

c           Write into tt value of derivative of j-th Lsf
c           evaluated at i-th root. This relation is described in the paper mentioned

            tt(j,i) = ww(j)*gg/ww(i)   
c           this appears wrong going by the given paper 
            
 3       continue
c        tt of i,i is 0 as the i-th Lsf has a maximum at the i-th root, unless i=1 or i=n
         tt(i,i) = ff
         
 4    continue
      do 7 i = 1,m
c        In this approach 1: roots are in decreasing order
c        2: last root ( 0 ) doesn't get included in calculations
c        but subroutine lobatto returns roots and weights in 
c        increasing order so this has to be done manually
         w(i) = ww(i+1)
         x(i) = xx(i+1)

         do 6 j = 1,i
         
            hh = 0.0d0
            do 5 k = 1,n
c			   Entries in T matrix are defined as a sum over all k from 0 
c           to n+1 of: 

c           [ k-th weight ]*[ derivative 
c			   of i-th Lsf at k-th root ] *[ derivative 
c			   of j-th Lsf at k-th root ]

               hh = hh + tt(k,i+1)*tt(k,j+1)
 5          continue
c           t is symmetric
            t(i,j) = hh
            t(j,i) = hh
            
 6       continue
 
 7    continue
      return
      end
      


      subroutine lobatto (a,b,n,w,x)
      implicit double precision (a-h,o-z)   
c     ----------------------------------------------------------------- 
c     This subroutine calculates an n-point Gauss-Lobatto
c     quadrature rule in the interval a < x < b.
c
c     Function localizes zeros of the derivative of n-1 th Legendre 
c     polynomial and calculates associated weights.
c
c     All recursive formulas are on Wikipedia except for P[n]''
c     but it can be derived using others
c
c     ----------------------------------------------------------------- 
c     
      dimension w(n),x(n)
  
c     shift and scale have to be used to change integral in (a,b)
c     to (-1,1) where lobatto quadrature works.
c     I denote n-th Legendre polynomial as P[n] and its derivative as P[n]'
c
      l = (n+1)/2
c     Isn't pi defined above as a double ?!?!?!
      pi = acos(-1.0d0)
c     ---<See Docs/GaussianQuadrature.pdf, top of page 2>
      shift = 0.5d0*(b+a)
      scale = 0.5d0*(b-a)
      weight = (b-a)/(n*(n-1))

      
c     Specific to Lobatto quadrature, First point is a
      x(1) = a
      w(1) = weight
      
      do 3 k = 2,l
c     As zeros are symmetric, there is only need to find positive ones
c     z is approximated zero of P[n-1] using Francesco Tricomi approximation
c     then accuracy of the zero is improved using Newton-Raphson
c     few times ( arbitrary so far ). The cosine equation finds a point
c     roughly halfway between 2 zeros of P[n].
         z = cos(pi*(4*k-3)/(4*n-2))
c        Calculate value of P[n-1] at z using Bonnets recursive formula
c        p1 is P[j], p2 = P[j-1], p3 = P[j-2]
         do 2 i = 1,7
            p2 = 0.0d0
            p1 = 1.0d0
            do 1 j = 1,n-1
               p3 = p2
               p2 = p1
               p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
 1          continue
c           p2 gets overwritten to be P[n-1]'
            p2 = (n-1)*(p2-z*p1)/(1.0d0-z*z)
c           p3 gets overwritten to be P[n-1]''
            p3 = (2.0d0*z*p2-n*(n-1)*p1)/(1.0d0-z*z)
c           Actual Newton-Raphson step
            z = z-p2/p3
 2       continue
c        Write in shifted and scaled zeros and weights
         x(k) = shift-scale*z
c        Write in zero with other sign
         x(n+1-k) = shift+scale*z
c        Write in weights (they are always positive)
         w(k) = weight/(p1*p1)
         w(n+1-k) = w(k)
 3    continue
c     Specific to Lobatto quadrature, last point is b
      x(n) = b
      w(n) = weight
      return
      end
       
      subroutine waves (w0,a,b,c,zmax)
      implicit double precision (a-h,o-z)
c
c     ------------------------------------------------------------------
c     Construction and storage of the diagonal matrices a,b and c
c     that enter the log derivative Kohn expression for the S-matrix.
c     ------------------------------------------------------------------
c
      complex*16 a, b, c
c
      dk = sqrt (dabs(w0))
      if (w0 .lt. 0.0d0) then
         theta = dk*zmax
         bcc   = cos(2.0d0*theta)
         bcs   = sin(2.0d0*theta)
         cc    = cos(theta)
         cs    = sin(theta) 
         a  = cmplx(bcc,-bcs)
         b  = (dk**0.5d0)*cmplx(cc,-cs)
         c  = cmplx(0.0d0,dk)
      else
         a = (0.0d0,0.0d0)
         b = (0.0d0,0.0d0)
         c = cmplx(-dk,0.0d0)
      endif
      return
      end
 
      subroutine precon (m,n,vfc,nfc,nfc00,d,e,f,t)
      implicit double precision (a-h,o-z)
      include 'multiscat.inc'
c
c     ------------------------------------------------------------------
c     This subroutine constructs the matrix factors that are required
c     for the block lower triangular preconditioner used in GMRES.
c
c	   t is the matrix from tshapes
c     ------------------------------------------------------------------
c
      complex*16 vfc(m,nfc)
      dimension d(n), e(m), f(m,n), t(m,m)
      dimension g(mmax)

c	   if m exceeds maximum size of m program terminates printing out precon 1
      if (m .gt. mmax) stop 'precon 1'
c
      do k = 1,m
c        dble(x) transfroms input into double precision real
c        for complex numbers it will return only real part
c        This overwrites t to be H0 (as named in '90 Kohn paper)
         t(k,k) = t(k,k) + dble(vfc(k,nfc00))
      enddo

c	   get eigenvalues [e] and eigenvectors [ overwrite them on t] of t, 
c	   which is m x m symmetric real matrix. f is temporary storage array
      call rs (m,m,t,e,t,f,ierr)     
c	   if ierr != 0 program terminates printing out precon 2, meaning rs failed
      if (ierr .ne. 0) stop 'precon 2' 
      
      do j = 1,n
      
         do k = 1,m
            g(k) = t(m,k)/(d(j)+e(k))
            f(k,j) = 0.0d0
         enddo
         
         do i = 1,m
            do k = 1,m
               f(k,j) = f(k,j) + t(k,i)*g(i)
            enddo
         enddo
         
      enddo
      return
      end
 
      subroutine gmres (x,xx,y,m,ix,iy,n,n00,vfc,ivx,ivy,nfc,
     +                  a,b,c,d,e,f,p,s,t,eps,ipc,ifail)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Complex Generalised Minimal Residual Algorithm (GMRES)
c     This version written by DEM, 6/12/94
c     ----------------------------------------------------------------- 
c
      complex*16 x(m*n), y(m*n), vfc(m,nfc)
      complex*16 a(n), b(n), c(n), s(n)
      dimension d(n), e(m), f(m,n), p(n), t(m,m)
      dimension ix(n), iy(n), ivx(nfc), ivy(nfc)
c
c     NB:
c     This subroutine implements a preconditioned version of GMRES(l).
c     The rate of convergence can be improved (at the expense of a 
c     greater disk space requirement and more cpu time per iteration) 
c     by increasing the following parameter:
c
c      parameter (l = 100)  increased l to 200, 22/4/1999, APJ
      parameter (l = 2000)
      complex*16 h(l+1,l+1), g(l+1), z(l+1)
      complex*16 co(l+1), si(l+1), temp
c
c     store x matrices in xx rather than write to disk.
c
      complex*16 xx(m*n,l+1)
c
c     Setup for GMRES(l):
c
      mn = m*n
      do i = 1,mn
         x(i) = (0.0d0,0.0d0) 
      enddo

c
c     Initial step:
c
      kount = 0
      xx(1:mn,1)=x(1:mn)
      do i = 1,mn
         y(i) = x(i)
      enddo
      call upper (x,m,ix,iy,n,vfc,ivx,ivy,nfc)
      do i = 1,mn
         x(i) = -x(i)
      enddo
      x(m*n00) = b(n00)+x(m*n00)
      call lower (x,m,ix,iy,n,vfc,ivx,ivy,nfc,c,d,e,f,t)
      do i = 1,mn
         x(i) = x(i)-y(i)
      enddo
      if (ipc .eq. 1) then
         do i = 1,mn
            y(i) = x(i)
         enddo
         call upper (x,m,ix,iy,n,vfc,ivx,ivy,nfc)
         call lower (x,m,ix,iy,n,vfc,ivx,ivy,nfc,c,d,e,f,t)
         do i = 1,mn
            x(i) = y(i)-x(i)
         enddo
      endif
      xnorm = 0.0d0
      do i = 1,mn
         xnorm = xnorm+conjg(x(i))*x(i)
      enddo
      xnorm = sqrt(xnorm)
      g(1) = xnorm
c
c     Generic recursion: (you mean iteration?)
c
      kconv = 0
      do j = 1,n
         p(j) = 0.0d0
      enddo
      do k = 1,l
         kount = kount+1
         do i = 1,mn
            x(i) = x(i)/xnorm
         enddo
         xx(1:mn,k+1)=x(1:mn)
         do i = 1,mn
            y(i) = x(i)
         enddo
         call upper (x,m,ix,iy,n,vfc,ivx,ivy,nfc)
         call lower (x,m,ix,iy,n,vfc,ivx,ivy,nfc,c,d,e,f,t)
         do i = 1,mn
            x(i) = y(i)+x(i)
         enddo
         if (ipc .eq. 1) then
            do i = 1,mn
               y(i) = x(i)
            enddo
            call upper (x,m,ix,iy,n,vfc,ivx,ivy,nfc)
            call lower (x,m,ix,iy,n,vfc,ivx,ivy,nfc,c,d,e,f,t)
            do i = 1,mn
               x(i) = y(i)-x(i)
            enddo
         endif
         y(1:mn)=xx(1:mn,1)
         do i = 1,n
            s(i) = y(m*i)
         enddo
         do j = 1,k
            y(1:mn)=xx(1:mn,j+1)
            h(j,k) = (0.0d0,0.0d0)
            do i = 1,mn
               h(j,k) = h(j,k)+conjg(y(i))*x(i)
            enddo
            do i = 1,mn
               x(i) = x(i)-y(i)*h(j,k)
            enddo
            if (j .lt. k) then
               do i = 1,n
                  s(i) = s(i)+y(m*i)*z(j)
               enddo
            endif
         enddo
         do i = 1,n
            s(i) = (0.0d0,2.0d0)*b(i)*s(i)
         enddo
         s(n00) = a(n00)+s(n00)
         xnorm = 0.0d0
         do i = 1,mn
            xnorm = xnorm+conjg(x(i))*x(i)
         enddo
         xnorm = sqrt(xnorm)
         h(k+1,k) = xnorm
         do j = 1,k-1
            temp = co(j)*h(j,k)+conjg(si(j))*h(j+1,k)
            h(j+1,k) = conjg(co(j))*h(j+1,k)-si(j)*h(j,k)
            h(j,k) = temp 
         enddo
         call zrotg (h(k,k),h(k+1,k),co(k),si(k))
         g(k+1) = -si(k)*g(k)
         g(k) = co(k)*g(k)
         do j = 1,k
            z(j) = g(j)
         enddo
         do j = k,1,-1
            z(j) = z(j)/h(j,j)
            do i = 1,j-1
               z(i) = z(i)-h(i,j)*z(j)
            enddo
         enddo
c
c        Convergence test:
c
         unit = 0.0d0
         diff = 0.0d0
         do j = 1,n
            pj = conjg(s(j))*s(j)
            unit = unit+pj
            diff = max(diff,abs(pj-p(j)))
            p(j) = pj
         enddo
         diff = max(diff,abs(unit-1.0d0))
         if (diff .lt. eps) then
            kconv = kconv+1
         else
            kconv = 0
         endif
         kk = k
         if (kconv.eq.3 .or. xnorm.eq.0.0d0) go to 2
      enddo
   2  continue
c
c     back substitution for x:
c
      x(1:mn)=xx(1:mn,1)
      do j = 1,kk
         y(1:mn)=xx(1:mn,j+1)
         do i = 1,mn
            x(i) = x(i)+y(i)*z(j)
         enddo
      enddo
c
c     all done? 
c
      if (kconv.lt.3 .and. xnorm.gt.0.0d0) then
         ifail=1
      endif
c
c     yes!
c
      return
      end

      subroutine upper (x,m,ix,iy,n,vfc,ivx,ivy,nfc)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------
c     This subroutine performs the block upper triangular 
c     matrix multiplication y = U*x, where A = L+U. 
c     The result y is overwritten on x on return. 
c     ----------------------------------------------------------
c
      complex*16 x(m,n), vfc(m,nfc)
      dimension ix(n), iy(n), ivx(nfc), ivy(nfc)
c
      do j = 1,n
         do k = 1,m
            x(k,j) = (0.0d0,0.0d0)
         enddo
         do i = j+1,n
            do l = 1,nfc
               if (ix(i) + ivx(l) .ne. ix(j)) go to 1
               if (iy(i) + ivy(l) .ne. iy(j)) go to 1
                  do k = 1,m
                     x(k,j) = x(k,j)+vfc(k,l)*x(k,i)
                  enddo
   1           continue
            enddo
         enddo
      enddo
      return
      end

      subroutine lower (x,m,ix,iy,n,vfc,ivx,ivy,nfc,c,d,e,f,t)
      implicit double precision (a-h,o-z)
      include 'multiscat.inc'
c
c     ----------------------------------------------------------
c     This subroutine solves the block lower triangular 
c     linear equation L*y = x, where A = L+U. 
c     The result y is overwritten on x on return.
c     ----------------------------------------------------------
c
      complex*16 x(m,n), vfc(m,nfc), c(n)
      dimension ix(n), iy(n), ivx(nfc), ivy(nfc)
      dimension d(n), e(m), f(m,n), t(m,m)
c
      
      !parameter (mmax = 200)
      complex*16 y(mmax), fac
      if (m .gt. mmax) stop 'lower 1'
c
      do j = 1,n
         do i = 1,j-1
            do l = 1,nfc
               if (ix(i) + ivx(l) .ne. ix(j)) go to 1
               if (iy(i) + ivy(l) .ne. iy(j)) go to 1
                  do k = 1,m
                     x(k,j) = x(k,j)-vfc(k,l)*x(k,i)
                  enddo
   1           continue
            enddo
         enddo
         do k = 1,m
            y(k) = (0.0d0,0.0d0)
            do l = 1,m
               y(k) = y(k)+x(l,j)*t(l,k)
            enddo
            y(k) = y(k)/(d(j)+e(k))
         enddo
         do k = 1,m
            x(k,j) = (0.0d0,0.0d0)
         enddo
         do l = 1,m
            do k = 1,m
               x(k,j) = x(k,j)+t(k,l)*y(l)
            enddo
         enddo
         fac = x(m,j)*c(j)/(1.0d0-f(m,j)*c(j))
         do k = 1,m
            x(k,j) = x(k,j)+fac*f(k,j)
         enddo
      enddo
      return
      end
 
      subroutine zrotg (a,b,c,s)
      implicit complex*16 (a-h,o-z)
      double precision scale,r
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs and performs a complex Givens rotation
c     ----------------------------------------------------------------- 
c
      rho = b
      if (abs(a) .gt. abs(b)) rho = a
      scale = abs(a) + abs(b)
      if (scale .eq. 0.0d0) then
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
      else 
         r = (a/scale)*conjg(a/scale) + (b/scale)*conjg(b/scale)
         r = scale*sqrt(r)
         if (dble(rho) .lt. 0.0d0) r = -r
         c = conjg(a/r)
         s = b/r
      end if
      z = 1.0d0
      if (abs(a) .gt. abs(b)) z = s
      if (abs(b) .ge. abs(a) .and. abs(c) .ne. 0.0d0) z = 1.0d0/c
      a = r
      b = z
      return
      end
 
      subroutine output (ei,theta,phi,ix,iy,n,n00,d,p,itest)

      implicit double precision (a-h,o-z)
c
c     -------------------------------------------------------------------
c     This subroutine outputs specular scattering probabilities 
c     for a given set of beam parameters ei, theta, and phi.
c     
c     The output file must have reference 21
c     -------------------------------------------------------------------
c
      dimension ix(n), iy(n), d(n), p(n)
      sum = 0.0d0
      if (itest.eq.1) then
         write (21,601) ei,theta,phi
      endif 
      do j = 1,n
         jx = ix(j)
         jy = iy(j)
         if (d(j) .lt. 0.0d0) then
            sum = sum + p(j)
            if (itest.eq.1) then
               write (21,602) jx,jy,p(j)
            endif
         endif
      end do
      if (itset.eq.1) write (21,612) sum
         write (*,'(5e14.6)') ei,theta,phi,p(n00),sum
      return
 601  format ('# Diffraction Probabilities at:'/'#',1x,36('-')/1x,
     + '# Beam energy           ei = ',e10.4/1x,
     + '# Polar angle        theta = ',5e14.6/1x,
     + '# Azimuthal angle      phi = ',e10.4/1x,36('-')//1x,
     + '# Diffraction Spot         Intensity'/'#',1x,36('-'))
 602  format ('#',i7,i6,e22.6)
 612  format ('#',37('-')/5x,'Unitarity',e22.6/'#',1x,36('-'))
      end