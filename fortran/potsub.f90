!Subroutines for linking multiscat core with periodic stepped potential data
!Andy Jardine, 12 January 1999 -->

!Modified fet Dec 07 and Dec 08
!Minor modifications F.Bello and E.Pierzchala Aug 2020

! ********************************************************************************
! loadfixedpot - loads potential data from file produced by Matlab and stores data
!                in memory. Also loads the number of fourier components available
!                and number of z points.
!
! The data is to be stored as two numbers (a single complex value) per line.  
! Fortran representation is used - (a, b) - where a is the real part and b is the 
! imaginary part.  A total of nzfixed*nfcfixed lines should be present, produced 
! by the matlab script 'four.m'
! DOCUMENTATION ERROR: 'four.m' seems to be depricated, it is now 'multiscat.m'
! All z values are sequential, that is each whole basis function data is together,
! going from minimum z to maximum z, before progressing to the next FC.

subroutine loadfixedpot(nzfixed,nfc,ivx,ivy,nfc00,vfcfixed,fourierfile,ax,ay,bx,by,zmin,zmax)

  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  integer        i,j                               !loop indecies
  integer        kx,ky
  integer        nkx,nky
  integer        nfc_from_file
  integer        header_nzfixed
  integer        nzfixed                             !number of z values in fixed fourier components
  integer        nfc                           !number of fourier components
  integer        nfc00
  integer        ivx(nfcx), ivy(nfcx)
  complex*16     vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX) !Fixed Fourier component data
  character*40   fourierfile                         !Fourier component data file                          
  character*200  header_line
  common /const/hemass,rmlmda !modified by Boyao on 6 Dec 2020

  ! Initiates the fourier components as a matrix of zeros
  vfcfixed=0.0d0
  !open the data file and read in the fourier components
  open(20,file=fourierfile)
  read(20,'(A)') header_line
  read(20,'(A)') header_line
  read(20,*) nfc_from_file, nkx, nky, header_nzfixed
  read(20,'(A)') header_line
  read(20,*) ax, ay, bx, by
  read(20,'(A)') header_line
  read(20,*) zmin, zmax
  read(20,'(A)') header_line
  read(20,'(A)') header_line
  read(20,'(A)') header_line

  nfc = nfc_from_file
  nzfixed = header_nzfixed

  if (nfc.ne.(nkx*nky)) then
    print *, 'ERROR: inconsistent potential header; nfc must equal nkx*nky.'
    stop
  end if

  nfc00 = 1
  i = 0
  do kx=0,nkx-1
    do ky=0,nky-1
      i = i + 1
      ivx(i) = kx
      if (kx.gt.((nkx-1)/2)) ivx(i) = kx-nkx
      ivy(i) = ky
      if (ky.gt.((nky-1)/2)) ivy(i) = ky-nky
      if ((ivx(i).eq.0).and.(ivy(i).eq.0)) nfc00 = i
    end do
  end do

  do i=1,nfc           !loop over fourier components
    do j=1,nzfixed     !loop over z values in fourier components
      read (20,*) vfcfixed(j,i)
    end do
  end do
  close(20)

  !Scale to the program units
  vfcfixed = vfcfixed * rmlmda
end subroutine loadfixedpot

