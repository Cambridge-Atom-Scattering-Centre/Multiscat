! Multiscat: Fast Close Coupled Scattering Calculation Program
!
! This version modified by Andy Jardine, to perform 3d scattering calculation
! for a surface specified with a fourier transform of a square lattice
!
! Converted to f90 free format 9th May 2001
! Further modified by fay summer 2009
! Modified by F.Bello and E. Pierzchala summer 2020

program multiscat
  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  !Define filenames
  character*40 inputfile,outfile,fourierfile, fourierLabelsFile, scattCondFile
      
  !Arrays
  complex*16 x(mmax,nmax), y(mmax,nmax), vfc(mmax,nfcx)
  parameter (lmax=901)                          !gmres solver stograge
  complex*16 xx(nmax*mmax,lmax)
  complex*16 a(nmax), b(nmax), c(nmax), s(nmax)
      
  !More Arrays
  dimension icx(nmax,ncmx), icy(nmax,ncmx)
  dimension ix(nmax), iy(nmax), ivx(nfcx), ivy(nfcx), ie(nmax,ncmx)
  dimension p(nmax), w(mmax), z(mmax)
  dimension ep(nmax,ncmx),eep(nmax,ncmx),cp(nmax,ncmx),dw(nmax,ncmx)
  dimension d(nmax), e(mmax), f(mmax,nmax), t(mmax,mmax)
  dimension ividx(nfcx),ivflag(nfcx)
  dimension ctheta(ncmx),cphi(ncmx),cei(ncmx)
  parameter (hbarsq = 4.18020)
  integer startindex,endindex !start and ending indexes of the potential files to be used 
  integer endOfFile
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
  
  !Load filenames for fourier labes and conditions
  read (80,*) fourierLabelsFile
  print *, 'Fourier labels file = ', fourierLabelsFile
  read (80,*) scattCondFile
  print *, 'Loading scattering conditions from ', scattCondFile

  open (81, file=scattCondFile)
  read (81, *)!Skip the first line of conditions file

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

  ! Checks that parameters don't clash
  if (nfc .gt. nfcx) then 
    print *, 'ERROR: the .conf file needs more fourier components', &
    ' than allowed by the .inc file (nfc>nfcx)'
    stop
  else if (nzfixed .gt. NZFIXED_MAX) then !not sure what something is!
    print *, 'ERROR: the .conf file needs more (something) than', &
    ' allowed by the .inc file (nzfixed>NZFIXED_MAX)'
    stop
  else if (nfc .gt. NVFCFIXED_MAX) then !not sure what something is!
    print *, 'ERROR: the .conf file needs more (something) than', &
    ' allowed by the .inc file (nfc>NVFCFIXED_MAX)'
    stop
  end if

  !Label the fourier components-they are listed in 'FourierLabels' and appear 
  !in the same order as in the potential file
  open (98, file=fourierlabelsfile, status='old')

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
  
    call loadfixedpot(nzfixed,nfc,vfcfixed,fourierfile)
    !this will read in the potential Fourier components and convert to the program units
  
  !========Do the scaterring calculations=========================================
    !Calculate scattering over the incident conditions required
    print *, ''
    print *, 'Calculating scattering for potential:',fourierfile
    print *, 'Energy / meV    Theta / deg    Phi / deg        I00         Sum ' 
   
    do
      read (81, *, iostat=endOfFile) ei, theta, phi !iostat checks for the end of the file
      if (endOfFile==0) then !Normal input
        !print *, ei, theta, phi
          
        !find number of z values required
        call findmz (emax,vmin,nsf,zmin,zmax,m)
        if (itest.eq.1) write(21,*) 'Required number of z grid points, m = ',m
        if (m.gt.mmax) stop 'ERROR: m too big!'
        !calculation kinetic enery matrix (possibly?)    
        call tshape (zmin,zmax,m,w,z,t)
      
        !interpolate vfcs to required z positions
        call potent(stepzmin,stepzmax,nzfixed,vfcfixed,nfc,vfc,m,z,ividx,ivflag,ivx)
    
        !get reciprocal lattice points    (also calculate how many channels are required for the calculation) 
        call basis(d,ix,iy,n,n00,dmax,imax)
        if (itest.eq.1) write(21,*) 'Number of diffraction channels, n =',n
        if (n.gt.nmax) stop 'ERROR: n too big!'
    
        !routines for actually doing the calculation
        do i = 1,n
          call waves (d(i),a(i),b(i),c(i),zmax)
          b(i) = b(i)/w(m)
          c(i) = c(i)/(w(m)**2)
          ! print *, b(i), c(i)
        end do
        call precon (m,n,vfc,nfc,nfc00,d,e,f,t)
        ifail=0
        call gmres  (x,xx,y,m,ix,iy,n,n00,vfc,ivx,ivy,nfc,a,b,c,d,e,f,p,s,t,eps,ipc,ifail)
    
        !if failure, then put all intensity to -1
        if (ifail.eq.1) then
          p=-1
        end if
        
        ! write outputs 
        call output(ei,theta,phi,ix,iy,n,n00,d,p,itest)
    
      else if (endOfFile<0) then !End of file
        print *, '-- End of scattering conditions file --'
        if (itest.eq.1) close (21)
        exit
      else !Unknown error
        print *, '#### ERROR: Invalid line found in input file ####'
        stop
      end if
    end do
  end do  
end program multiscat

