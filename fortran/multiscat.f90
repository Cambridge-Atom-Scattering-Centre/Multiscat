! Multiscat: Fast Close Coupled Scattering Calculation Program
!
! This version modified by Andy Jardine, to perform 3d scattering calculation
! for a surface specified with a fourier transform of a square lattice
!
! Converted to f90 free format 9th May 2001
! Further modified by fay summer 2009
! Modified by F.Bello and E. Pierzchala summer 2020

program multiscat
  use multiscat_io_loaders
  implicit double precision (a-h,o-z)
  include 'multiscat.inc'

  !Define filenames
  character*40 optimizationFile,outfile,fourierfile, scattCondFile
  character*80 arg
      
  !Arrays
  complex*16 x(mmax,nmax), y(mmax,nmax), vfc(mmax,nfcx)
  parameter (lmax=901)                          !gmres solver stograge
  complex*16 xx(nmax*mmax,lmax)
  complex*16 a(nmax), b(nmax), c(nmax), s(nmax)
      
  !More Arrays
  dimension ix(nmax), iy(nmax), ivx(nfcx), ivy(nfcx)
  dimension p(nmax), w(mmax), z(mmax)
  dimension d(nmax), e(mmax), f(mmax,nmax), t(mmax,mmax)
  parameter (hbarsq = 4.18020)
  integer argc, iarg
  integer icond

  type(OptimizationFileData) :: optimization_data
  type(ScatteringConditionsData) :: scatt_conditions_data
  type(FixedPotentialData) :: potential_data

  !Variables for potential, represented as fourier data
  complex*16 vfcfixed(NZFIXED_MAX,NVFCFIXED_MAX)   !FC's at the fixed points

  common /const/ hemass, rmlmda
  !common /const/ rmlmda !commented by Boyao on 6 Dec 2020
  common /cells/ ax,ay,bx,by,ei,theta,phi,a0


  !===========================================================================


  !Begin the main program
  print *, ''
  print *, 'Multiscat: Close Coupled Scattering Program'
  print *, '============================================='
  print *, ''

  optimizationFile = ''
  fourierfile = ''
  scattCondFile = ''

  ! Parse required CLI flags.
  argc = command_argument_count()
  iarg = 1
  do while (iarg.le.argc)
    call getarg(iarg,arg)
    if (trim(arg).eq.'--optimization') then
      iarg = iarg + 1
      if (iarg.gt.argc) stop 'Error: --optimization requires a file path.'
      call getarg(iarg,optimizationFile)
    else if (trim(arg).eq.'--potential') then
      iarg = iarg + 1
      if (iarg.gt.argc) stop 'Error: --potential requires a file path.'
      call getarg(iarg,fourierfile)
    else if (trim(arg).eq.'--condition') then
      iarg = iarg + 1
      if (iarg.gt.argc) stop 'Error: --condition requires a file path.'
      call getarg(iarg,scattCondFile)
    else
      stop 'Error: unrecognized argument. Use --potential, --condition, --optimization.'
    end if
    iarg = iarg + 1
  end do

  if (optimizationFile.eq.'') stop 'Error: you must supply --optimization <file>.'
  if (fourierfile.eq.'') stop 'Error: you must supply --potential <file>.'
  if (scattCondFile.eq.'') stop 'Error: you must supply --condition <file>.'

  print *, 'Reading optimization parameters from input file: ',optimizationFile
  print *, 'Loading scattering conditions from ', scattCondFile
  print *, 'Calculating for potential input file ',trim(fourierfile)
  print *, ''

  !=====================read in parameters from config file==========================
  call load_optimization_file(optimizationFile, optimization_data)
  call load_scattering_conditions_file(scattCondFile, scatt_conditions_data)
  hemass = scatt_conditions_data%helium_mass
  rmlmda = 2.0d0*hemass/hbarsq
  call load_fixed_potential(fourierfile, rmlmda, potential_data)

  itest = optimization_data%output_mode
  ipc = optimization_data%gmres_preconditioner_flag
  nsf = optimization_data%convergence_significant_figures
  eps = 0.5d0*(10.0d0**(-nsf))
  dmax = optimization_data%max_closed_channel_energy
  imax = optimization_data%max_channel_index
  
!===============preliminary calculation and setting up ===========================

  iread=5
  iwrite=6
  ireadp=10
  ireadc=10
  ireade=10
  iwritep=10
  iwritel=11
  ireadip=12

! ============================================================================
  if (itest.eq.1) then
    outfile='diffrac.out'
    ! diffrac will be the output file containing diffraction calculations;
    open(21,file=outfile,status='unknown')
    write(21,*) 'Diffraction intensities for potential:',fourierfile
  end if
      
  !========Initialize the potential================================================
  

    !this will read in the potential Fourier components and convert to the program units

    nfc = potential_data%fourier_component_count
    nzfixed = potential_data%z_point_count
    nfc00 = potential_data%specular_component_index
    ax = potential_data%unit_cell_ax
    ay = potential_data%unit_cell_ay
    bx = potential_data%unit_cell_bx
    by = potential_data%unit_cell_by
    zmin = potential_data%zmin
    zmax = potential_data%zmax

    ivx(1:nfc) = potential_data%fourier_indices_x(1:nfc)
    ivy(1:nfc) = potential_data%fourier_indices_y(1:nfc)
    vfcfixed(1:nzfixed,1:nfc) = potential_data%fixed_fourier_values(1:nzfixed,1:nfc)
  
  !========Do the scaterring calculations=========================================
    !Calculate scattering over the incident conditions required
    print *, ''
    print *, 'Calculating scattering for potential:',fourierfile
    print *, 'Energy / meV    Theta / deg    Phi / deg        I00         Sum ' 
   
    do icond = 1, scatt_conditions_data%condition_count
      ei = scatt_conditions_data%incident_energy_mev(icond)
      theta = scatt_conditions_data%theta_degrees(icond)
      phi = scatt_conditions_data%phi_degrees(icond)

        !Use the Lobatto z grid size from the loaded potential file.
        m = nzfixed
        if (itest.eq.1) write(21,*) 'Using z grid points from potential file, m = ',m
        if (m.gt.mmax) stop 'ERROR: m too big!'
           
        call tshape (zmin,zmax,m,w,z,t)

        !Input potential is already provided on this Lobatto z grid.
        do i=1,nfc
          do j=1,m
            vfc(j,i)=vfcfixed(j,i)
          end do
        end do
    
        !get reciprocal lattice points    (also calculate how many channels are required for the calculation) 
        call basis(d,ix,iy,n,n00,dmax,imax)
        if (itest.eq.1) write(21,*) 'Number of diffraction channels, n =',n
        if (n.gt.nmax) stop 'ERROR: n too big!'
    
        !routines for actually doing the calculation
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
        call output(ei,theta,phi,ix,iy,n,n00,d,p,itest)
    end do

    print *, '-- End of scattering conditions file --'
    if (itest.eq.1) close (21)
 
end program multiscat

