! Multiscat: Fast Close Coupled Scattering Calculation Program
!
! This version modified by Andy Jardine, to perform 3d scattering calculation
! for a surface specified with a fourier transform of a square lattice
!
! Converted to f90 free format 9th May 2001
! Further modified by fay summer 2009

program multiscat
  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  !Define filenames
  character*40 inputfile,outfile,fourierfile
      
  !Arrays
  complex*16 x(mmax,nmax), y(mmax,nmax), vfc(mmax,nfcx)
  parameter (lmax=901)                          !gmres solver stograge
  complex*16 xx(nmax*mmax,lmax)
  complex*16 a(nmax), b(nmax), c(nmax), s(nmax)
      
  !More Arrays
  dimension icx(nmax,ncmx), icy(nmax,ncmx), ncn(ncmx), ncn00(ncmx)
  dimension ix(nmax), iy(nmax), ivx(nfcx), ivy(nfcx), ie(nmax,ncmx)
  dimension p(nmax), w(mmax), z(mmax)
  dimension ep(nmax,ncmx),eep(nmax,ncmx),cp(nmax,ncmx),dw(nmax,ncmx)
  dimension d(nmax), e(mmax), f(mmax,nmax), t(mmax,mmax)
  dimension ivsym(nfcx),ividx(nfcx),ivflag(nfcx)
  dimension ctheta(ncmx),cphi(ncmx),cei(ncmx), ne(ncmx)
  dimension xc(npmax), xc0(npmax),q(1000), ei_array(1000), phi_array(1000), theta_array(1000)
  dimension A0w(npwx),B0w(npwx),vmw(npwx),apw(npwx),bpw(npwx)
  parameter (hbarsq = 4.18020)
  integer startindex,endindex !start and ending indexes of the potential files to be used 
  !Variables for potential, represented as fourier data
  complex*16 vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX)   !FC's at the fixed points

  common /const/ hemass
  common /const/ rmlmda
  common /cells/ a1,a2,b2,ei,theta,phi,a0


  !===========================================================================


  !Begin the main program
  print *, ''
  print *, 'Multiscat: Close Coupled Scattering Program'
  print *, '============================================='
  print *, ''

  !get the name of the config file
  call getarg(1,inputfile)
  if (inputfile.eq.'') stop 'Error: you must supply a configuration file to run Multiscat.'
  print *, 'Reading parameters from input file: ',inputfile
  print *, ''

  !=====================read in parameters from config file==========================
 
  !read in parameters from the config file and make preliminary calculations
  open (80,file=inputfile) 
  read (80,*) itest
  print *, 'Output mode = ',itest
  read (80,*) ipc
  if (ipc.lt.0) ipc = 0   !only ipc = 0 and 1 are implemented in gmres:
  if (ipc.gt.1) ipc = 1
  print *, 'GMRES preconditioner flag = ',ipc
  read (80,*) nsf
  if (nsf.lt.2) nsf = 2   !place an upper and lower limit on the precision
  if (nsf.gt.5) nsf = 10
  eps = 0.5d0*(10.0d0**(-nsf))
  print *, 'Convergence sig. figures = ',nsf
  read (80,*) nfc
  print *, 'Total number of fourier components to use = ',nfc
  print *, ''
  read (80,*) zmin,zmax    !this is the required integration range; later we calculate how many points in the ramge are required and the potential is interpolated to those points
  print *, 'z integration range = (',zmin,',',zmax,')'
  read (80,*) vmin
  print *, 'Potential well depth = ',vmin
  read (80,*) dmax
  print *, 'Max energy of closed channels = ',dmax
  read (80,*) imax
  print *, 'Max index of channels = ',imax
  print *, ''

  ! read in the ranges of energy and scattering angles to be used
  
  read (80,*) eii,dei,nei 
  if (dei == 0 .AND. nei > 0) then
    print *, 'ERROR: Unsupported energy range found in .conf file.'
    print *, '', nei, ' steps requested with ', dei,' energy delta between them.'
    stop
  end if 

  READ(80,*)  (ei_array(i), i=1,0)
  read (80,*) thetai,dtheta,ntheta
  if (dtheta == 0 .AND. ntheta > 0) then
    print *, 'ERROR: Unsupported theta range found in .conf file.'
    print *, '', ntheta, ' steps requested with ', dtheta,' theta delta between them.'
    stop
  end if 

  READ(80,*)  (theta_array(i), i=1,0 )
  read (80,*) phii,dphi,nphi
  if (dphi == 0 .AND. nphi > 0) then
    print *, 'ERROR: Unsupported phi range found in .conf file.'
    print *, '', nphi, ' steps requested with ', dphi,' phi delta between them.'
    stop
  end if 

  READ(80,*)  (phi_array(i), i=1,0 )

  print *, 'Initial energy (meV) = ',eii
  print *, 'Energy change (meV) = ',dei
  print *, 'Num energy steps = ',nei
  print *, 'Energy array = ',(ei_array(i), i=1, nei)
  print *, 'Initial theta (deg) = ',thetai
  print *, 'Theta change (deg) = ',dtheta
  print *, 'Num theta steps = ',ntheta
  print *, 'Theta array = ',(theta_array(i), i=1, ntheta)
  print *, 'Initial phi (meV) = ',phii
  print *, 'Phi change (meV) = ',dphi
  print *, 'Num phi steps = ',nphi
  print *, 'Phi array = ',(phi_array(i), i=1, nphi)
  
  !read in and set shape of real space lattice; it is hexagonal lattice, but a1,a2 and b2 are
  ! its dimensions in cartesian coordinates 
  read (80,*) a1       !surface lattice constant in x direction (see basis)
  read (80,*) a2
  read (80,*) b2       !surface lattice constant in y direction
  print *, 'Unit cell (A) = ',a1,'x',b2

  read (80,*) nzfixed   ! number of z points in Fourier components of potential
  print *, 'Number of z points in fourier components (nzfixed) = ',nzfixed
  print *, ''           
  read (80,*) stepzmin  !maximum and minimumn values of z in the potential file read in
  read (80,*) stepzmax
  read(80,*) startindex !the start and end indices of the potential files to be used
  read(80,*) endindex
  print *, 'Calculating for potential input files between ',startindex,'.in and ',endindex,'.in'
  read(80,*) hemass
  
!===============preliminary calculation and setting up ===========================
  rmlmda = 2.0d0*hemass/hbarsq
  iread=5
  iwrite=6
  ireadp=10
  ireadc=10
  ireade=10
  iwritep=10
  iwritel=11
  ireadip=12
 
  !Label the fourier components-they are listed in 'FourierLabels' and appear in the same order as in the potential file
  open (98, file='FourierLabels.in', status='old')

  do i=1, nfc
     read (98,*)  ivx(i), ivy(i)
     if  ((ivx(i).eq.0) .and. (ivy(i).eq.0)) nfc00=i
  end do
  close (98)

! ============================================================================
!do loop for using different potential files
  do in=startindex,endindex
      write(fourierfile,599) in
599 format('pot',i5,'.in')
    if (itest.eq.1) write(outfile,598) in
598 format('diffrac',i5,'.out')
   ! diffrac will be the output file containing diffraction calculations;
  if (itest.eq.1) open(21,file=outfile,status='unknown')
  if (itest.eq.1) write(21,*) 'Diffraction intensities for potential:',fourierfile 
    
!========Initialize the potential================================================

  call loadfixedpot(stepzmin,stepzmax,nzfixed,nfc,vfcfixed,fourierfile,itest)
  !this will read in the potential Fourier components and convert to the program units

!========Do the scaterring calculations=========================================
  !Calculate scattering over the incident conditions required
  print *, ''
  print *, 'Calculating scattering for potential:',fourierfile
  print *, 'Energy / meV    Theta / deg    Phi / deg        I00         Sum ' 
  
 do iei=0,nei

    ei=eii+(iei*dei)
    do iphi=0,nphi
      phi=phii+(iphi*dphi)

      do itheta=0,ntheta
        theta=thetai+(itheta*dtheta)
      
 !find number of z values required
    call findmz (emax,vmin,nsf,zmin,zmax,m,itest)
    if (itest.eq.1) write(21,*) 'Required number of z grid points, m = ',m
    if (m.gt.mmax) stop 'ERROR: m too big!'
 !calculation kinetic enery matrix (possibly?)    
    call tshape (zmin,zmax,m,w,z,t,itest)
  
 !interpolate vfcs to required z positions
    call potent(stepzmin,stepzmax,nzfixed,vfcfixed,nfc,vfc,m,z,ividx,ivflag,itest,ivx,ivy)

    !get reciprocal lattice points    (also calculate how many channels are required for the calculation) 
    call basis(d,ix,iy,n,n00,dmax,imax,iwrite,itest)
    if (itest.eq.1) write(21,*) 'Number of diffraction channels, n =',n
    if (n.gt.nmax) stop 'ERROR: n too big!'

!routines for actually doing the calculation
    print *, 'n = ', n, '|| m = ', m
    do i = 1,n
      call waves (d(i),a(i),b(i),c(i),zmax)
      b(i) = b(i)/w(m)
      c(i) = c(i)/(w(m)**2)

    end do
    call precon (m,n,vfc,nfc,nfc00,d,e,f,t)
    ifail=0
    call gmres  (x,xx,y,m,ix,iy,n,n00,vfc,ivx,ivy,nfc,a,b,c,d,e,f,p,s,t,eps,ipc,ifail)

    !if failure, then put all intensity to -1
    if (ifail.eq.1) then
      p=-1
    end if
    
! write outputs 
    call output(ei,theta,phi,ix,iy,n,n00,d,p,nsf,outfile,itest)

  end do             
  end do
  end do
  if (itest.eq.1) close (21)
  end do

end program multiscat

