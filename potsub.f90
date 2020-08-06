!Subroutines for linking multiscat core with periodic stepped potential data
!Andy Jardine, 12 January 1999 -->

!Modified fet Dec 07 and Dec 08

! ********************************************************************************
! loadfixedpot - loads potential data from file produced by Matlab and stores data
!                in memory.  'potent' is used to interpolate this data to the values
!                required by the rest of the program.  Also loads the number of
!                fourier components available and number of z points.
!
! The data is to be stored as two numbers (a single complex value) per line.  
! Fortran representation is used - (a, b) - where a is the real part and b is the 
! imaginary part.  A total of nzfixed*nfcfixed lines should be present, produced 
! by the matlab script 'four.m'
! DOCUMENTATION ERROR: 'four.m' seems to be depricated, it is now 'multiscat.m'
! All z values are sequential, that is each whole basis function data is together,
! going from minimum z to maximum z, before progressing to the next FC.

subroutine loadfixedpot(nzfixed,nfc,vfcfixed,fourierfile)

  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  integer        i,j                               !loop indecies
  integer        nzfixed                             !number of z values in fixed fourier components
  integer        nfc                           !number of fourier components
  complex*16     vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX) !Fixed Fourier component data
  character*40   fourierfile                         !Fourier component data file                          
  common /const/rmlmda

  ! Initiates the fourier components as a matrix of zeros
  vfcfixed=0.0d0
  !open the data file and read in the fourier components (1 DC component + nvfcfixed fourier components)
  open(20,file=fourierfile)
  !discard the first 5 lines
  read(20,*)
  read(20,*)
  read(20,*)
  read(20,*)
  read(20,*)

  do i=1,nfc           !loop over fourier components
    do j=1,nzfixed     !loop over z values in fourier components
      read (20,*) vfcfixed(j,i)
    end do
  end do

  !Scale to the program units
  vfcfixed = vfcfixed * rmlmda
end subroutine loadfixedpot

!************************************************************************
! potent  - interpolates the data from Matlab to the data points requested by
!           the call to tshape/findmz
!
!           The whole of the requested vfc matrix is generated here
!           Also, have to set which is the zero fourier cmpt (nfc00)

subroutine potent(stepzmin,stepzmax,nzfixed,vfcfixed,nfc,vfc,m,z,ividx,ivflag)

  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  complex*16 vfc(m,nfc)                          ! requested fourier cmpts
  complex*16 vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX)
  dimension z(m)                                 ! requested z points (m=num z points)
  dimension ividx(nfc), ivflag(nfc)
  complex*16 atmp, btmp, vrealtmp


 !generate vfc matrix by interpolation of vfcfixed

  do i=1,nfc                                ! loop over fourier components

    do j=1,m                                ! loop over reqested points
      
      ! Locate what would be the index in the list of z points
      zindex = (z(j)-stepzmin)/(stepzmax-stepzmin)*(nzfixed-1)+1
      indexlow=int(zindex)               ! truncate to integer
      
      ! Pick out the value we are interested in
      if (zindex.eq.indexlow) then
!        ! have got exact value - no need to interpolate
        vfc(j,i) = vfcfixed(indexlow,i)
      else
        ! need to interpolate - interpolate real and imaginary parts separately (does it?)
        
        ! Interpolate potential
        atmp = vfcfixed(indexlow,i)
        btmp = vfcfixed(indexlow+1,i)
        vrealtmp = atmp + (btmp-atmp) * (zindex-dfloat(indexlow))
        
        ! Store for use
        vfc(j,i) = vrealtmp
      end if

    end do
  
  end do
return
end subroutine potent



!*****************************************************************************
!onedbasis - this is a version of the basis subroutine, modified for operation
!            in one dimension, rather than two.  Also, code is converted to
!            Fortran90 (cos fortran77 is crap.)
!
!c     
!c     calculate reciprocal lattice, read in channels from chanfile,
!c     calculate d (z-component of energy of outgoing wave) for 
!c     each channel
!c     
!c     
!c     area      = area of real space unit cell
!c     
!c     gax, gbx  = the unit vector of reciprocal lattice along symmetry 
!c     direction
!c     
!c     gay, gby  = y components of unit vector of reciprocal lattice
!c     along symmetry direction 2
!c     
!c     a1     length of real space lattice vector along axis
!c     a2     x coordinate of other real space lattice vector
!c     b2     y coordinate of other real space lattice vector
!c     2
!c                                    y|   _____a1_____ 
!c     note that the a1,               |  /|          /
!c     and b2 are the                  | / |b2       /
!c     paramenters of the              |/__|________/____ x
!c     unit cell of substrate           a2    
!c     
!c
!c     For each scattered channel:
!c     ered = square of incident wavevector for the channel
!c     eint = squre of surface component of the wavevector for the 
!c     channel
!c     
!c     -d(i) = square of the z component of the wavevector for the channel
!c     energy conservation: ered = eint + (-d(i))
!c     if d(i) < 0, channel open, possible diffraction spot
!c     if d(i) > 0, channel closed, no spot
!c     
!
! Note:
! The format of the command has been left unchanged, so the main program is
! virtually unaltered.



subroutine onedbasis(d,ix,iy,n,n00,dmax,imax,iwrite,itest)
      
  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  character*20   chanfile
  dimension      d(nmax), ix(nmax), iy(nmax)
      
  common /cells/ a1,a2,b2,ei,theta,phi,a0,gax,gay,gbx,gby
  common /const/ rmlmda
  DATA   Pi /3.141592653589793d0/
      
  if (itest.eq.1) print *, 'Calculating reciprocal lattice points (onedbasis).'
  
  !Old code, to allow for generally shaped 2d basis      
  ax=a1
  ay=0
  bx=a2     !basis vector b is zero in one dimension
  by=b2
  Auc=dabs(ax*by)
  RecUnit=2*Pi/Auc
  gax =  by*RecUnit
  gay = -bx*RecUnit
  gbx = -ay*RecUnit
  gby =  ax*RecUnit
    
  !Basis reciprocal lattice vectors are just 2pi/real space vectors, as
  !unit cell only 1D now.
  !ax=a1
  !ay=0
  !gax=2*Pi/ax
  !gay=0
      
   
  if (itest.eq.1) then
    write(iwrite,'(a)') '# unit cell:'
    write(iwrite,'(a,e14.6,a,e14.6,a)')'# real :(',ax,',',ay,')'
    write(iwrite,'(a,e14.6,a,e14.6,a)')'#       (',bx,',',by,')'
    write(iwrite,'(a,e14.6,a,e14.6,a)')'# recip:(',gax,',',gay,')'
    write(iwrite,'(a,e14.6,a,e14.6,a)')'#       (',gbx,',',gby,')'
  end if

  !Convert incidence angle to radians for use
  thetad = theta*pi/180.0d0
  
  !Calculate the momentum (?)
  ered   = rmlmda*ei                          !scale energy of incident He
  phid   = phi*pi/180.0d0 
  !pkx = sqrt(ered)*sin(thetad)*cos(phid)
  !pky = sqrt(ered)*sin(thetad)*sin(phid)
  !pkx = sqrt(ered)*sin(thetad)                !
  
      
  n=0
  do i1 = -imax,imax
    !do i2 = -imax,imax
      !gx = gax*i1 + gbx*i2
      !gy = gay*i1 + gby*i2
      !eint = (pkx+gx)**2 + (pky+gy)**2
      
      gx = gax*i1
      eint = (pkx+gx)**2
      di = eint-ered
      if (di.lt.dmax) then 
        n=n+1
        if (n.le.nmax) then
          ix(n)=i1
          !iy(n)=i2
          d(n)=di
          if (i1.eq.0) n00=n
        else
          stop 'ERROR: n too big! (subroutine onedbasis)'
        end if
      end if
    !end do
  end do
      
return
end subroutine onedbasis
         








!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!onedoutput - Outputs the data to the screen in a more convenient format that
!             the standard output (ie: suitable for sending to xmgr)
!
!             Outputs all diffractive scattering probabilities, along with
!             the parallel momentum change

subroutine onedoutput (ei,theta,phi,ix,iy,n,n00,d,p,nsf,iwrite,itest,a1,inumenergy)

  implicit double precision (a-h,o-z)
  dimension ix(n), iy(n), d(n), p(n)

  !Define a sum, to calculate unitarity with
  sum = 0.0d0
  deltaG = 2.0d0*3.141592654d0/a1

  if (itest.eq.1) then

    print *, 'Storing scattering amplitudes to scattering.dat'
    open (82,file='scattering.dat')

    !Loop over channels, outputting data and adding to sum
    do j = 1,n
      jx = ix(j)
      !jy = iy(j)
      deltaK = dfloat(jx)*deltaG
      if (d(j) .lt. 0.0d0) then
        sum = sum + p(j)
        write (82,*) jx,deltaK,p(j)
        print *, jx,deltaK,p(j)
      endif
    end do

    close(82)

    !Print out the unitarity indicator to check the calculation went OK
    print *, 'Unitarity = ',sum

  else
    
    vectork = (2.0*(ei/1000.0*1.602177d-19)*6.6465e-27/(1.054573d-34)**2.0)**0.5 / 1.0d10
    specularG = -2.0*vectork*sin(abs(theta*3.141592654/180.0))

    !Loop over channels, outputting data and adding to sum
    do j = 1,n
      jx = ix(j)
      jy = iy(j)
      deltaK = dfloat(jx)*deltaG
      if (d(j) .lt. 0.0d0) then
        sum = sum + p(j)
        !if (abs(deltaK-specularG).lt.1.0) print *, ei,deltaK,jx,p(j),specularG
        !if (abs(deltaK-specularG).lt.0.2) write (iwrite,'(F8.3,F9.5,I7,F9.5,F9.5)') ei,deltaK,jx,p(j),specularG
        if (jx.eq.2*inumenergy) write (83,'(F10.5,F10.5,I7,F10.5)') ei,deltaK,jx,p(j)
        if (jx.eq.2*inumenergy) write (6,'(F10.5,F10.5,I7,F10.5,F10.5)') ei,deltaK,jx,p(j)
      endif
    end do

    !Print out the results
    print *, 'Energy =',ei,' Unitarity = ',sum


  end if

  return

end subroutine onedoutput
