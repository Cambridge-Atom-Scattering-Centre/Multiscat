!Subroutines for linking multiscat core with periodic stepped potential data
!Andy Jardine, 12 January 1999 -->

!Modified fet Dec 07 and Dec 08
!Minor modifications F.Bello and E.Pierzchala Aug 2020

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
  common /const/hemass,rmlmda !modified by Boyao on 6 Dec 2020

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

subroutine potent(stepzmin,stepzmax,nzfixed,vfcfixed,nfc,vfc,m,z)

  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  complex*16 vfc(m,nfc)                          ! requested fourier cmpts
  complex*16 vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX)
  dimension z(m)                                 ! requested z points (m=num z points)
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
