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
! All z values are sequential, that is each whole basis function data is together,
! going from minimum z to maximum z, before progressing to the next FC.

subroutine loadfixedpot(nzfixed,nfc,vfcfixed,fourierfile,itest)

  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  integer        i,j                               !loop indecies
  integer        nzfixed                             !number of z values in fixed fourier components
  integer        nfc                           !number of fourier components
  integer        itest                               !Test mode 1=yes, 0=no
  complex*16     vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX) !Fixed Fourier component data
  character*40   fourierfile                         !Fourier component data file                          
  common /const/rmlmda


  !print up an indicator that this subroutine has been called
  !if (itest.eq.1) print *, 'Loading fcs (loadfixedpot) (',nzfixed,',',nvfc_xfix*nvfc_yfix,')'


  !open the data file and read in the fourier components (1 DC component + nvfcfixed fourier components)
  vfcfixed=0.0d0
  open(20,file=fourierfile)
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
  !complex*16 vfcfixed(nzfixed, nvfcfixed)        ! fixed fourier cmpts
  complex*16 vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX)
  dimension z(m)                                 ! requested z points (m=num z points)
  dimension ividx(nfc), ivflag(nfc)
  complex*16 atmp, btmp, vrealtmp


 !generate vfc matrix by interpolation of vfcfixed

  do i=1,nfc                                ! loop over fourier components
    
   ! k=ividx(i)                              ! get which fixed component
     k=i                                        ! k is the fixed point reference
    do j=1,m                                ! loop over reqested points
      
      ! Locate what would be the index in the list of z points
      zindex = (z(j)-stepzmin)/(stepzmax-stepzmin)*(nzfixed-1)+1
      indexlow=int(zindex)                       ! truncate to integer
      
      ! Pick out the value we are interested in
      if (zindex.eq.indexlow) then
!        ! have got exact value - no need to interpolate
        vfc(j,i) = vfcfixed(indexlow,k)
        !vfc(j,i)=vfcfixed(indexlow,i)
      else
        ! need to interpolate - interpolate real and imaginary parts separately
        
        ! Interpolate potential
        atmp = vfcfixed(indexlow,k)
        !atmp = vfcfixed(indexlow,i)
        btmp = vfcfixed(indexlow+1,k)
!btmp = vfcfixed(indexlow+1,i)
        vrealtmp = atmp + (btmp-atmp) * (zindex-dfloat(indexlow))
        
        ! Store for use
        vfc(j,i) = vrealtmp
!         vfc(j,i)=vfcfixed(j,i)
      end if
      
      ! Make sure the complex part is the 'right way around'
!      if (ivflag(i).eq.-1) then
!        vfc(j,i) = conjg(vfc(j,i))
!      end if

    end do
  
  end do
return
end subroutine potent
