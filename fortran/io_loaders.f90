module multiscat_io_loaders
  implicit none

  type :: OptimizationFileData
    integer :: output_mode = 0
    integer :: gmres_preconditioner_flag = 0
    integer :: convergence_significant_figures = 2
    double precision :: max_closed_channel_energy = 0.0d0
    integer :: max_channel_index = 0
  end type OptimizationFileData

  type :: ScatteringConditionsData
    double precision :: helium_mass = 0.0d0
    integer :: condition_count = 0
    double precision, allocatable :: incident_energy_mev(:)
    double precision, allocatable :: theta_degrees(:)
    double precision, allocatable :: phi_degrees(:)
  end type ScatteringConditionsData

  type :: FixedPotentialData
    integer :: z_point_count = 0
    integer :: fourier_component_count = 0
    integer :: specular_component_index = 1
    integer :: fourier_grid_x_count = 0
    integer :: fourier_grid_y_count = 0
    double precision :: unit_cell_ax = 0.0d0
    double precision :: unit_cell_ay = 0.0d0
    double precision :: unit_cell_bx = 0.0d0
    double precision :: unit_cell_by = 0.0d0
    double precision :: zmin = 0.0d0
    double precision :: zmax = 0.0d0
    integer, allocatable :: fourier_indices_x(:)
    integer, allocatable :: fourier_indices_y(:)
    complex*16, allocatable :: fixed_fourier_values(:,:)
  end type FixedPotentialData

contains

  subroutine load_optimization_file(optimizationFile, data)
    implicit double precision (a-h,o-z)
    character(len=*), intent(in) :: optimizationFile
    type(OptimizationFileData), intent(out) :: data
    integer :: unit_number
    integer :: io_status

    open(newunit=unit_number, file=trim(optimizationFile), status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, 'ERROR: could not open optimization file: ', trim(optimizationFile)
      stop
    end if

    read(unit_number,*,iostat=io_status) data%output_mode
    if (io_status /= 0) stop 'ERROR: failed to read output mode from optimization file.'
    read(unit_number,*,iostat=io_status) data%gmres_preconditioner_flag
    if (io_status /= 0) stop 'ERROR: failed to read GMRES preconditioner flag from optimization file.'
    read(unit_number,*,iostat=io_status) data%convergence_significant_figures
    if (io_status /= 0) stop 'ERROR: failed to read convergence significant figures from optimization file.'
    read(unit_number,*,iostat=io_status) data%max_closed_channel_energy
    if (io_status /= 0) stop 'ERROR: failed to read max closed channel energy from optimization file.'
    read(unit_number,*,iostat=io_status) data%max_channel_index
    if (io_status /= 0) stop 'ERROR: failed to read max channel index from optimization file.'

    close(unit_number)

    if (data%gmres_preconditioner_flag < 0) data%gmres_preconditioner_flag = 0
    if (data%gmres_preconditioner_flag > 1) data%gmres_preconditioner_flag = 1
    if (data%convergence_significant_figures < 2) data%convergence_significant_figures = 2
    if (data%convergence_significant_figures > 5) data%convergence_significant_figures = 10

    print *, 'Output mode = ', data%output_mode
    print *, 'GMRES preconditioner flag = ', data%gmres_preconditioner_flag
    print *, 'Convergence sig. figures = ', data%convergence_significant_figures
    print *, ''
    print *, 'Max energy of closed channels = ', data%max_closed_channel_energy
    print *, 'Max index of channels = ', data%max_channel_index
    print *, ''
  end subroutine load_optimization_file

  subroutine load_scattering_conditions_file(scattCondFile, data)
    implicit none
    character(len=*), intent(in) :: scattCondFile
    type(ScatteringConditionsData), intent(out) :: data
    integer :: unit_number
    integer :: io_status
    integer :: i
    integer :: condition_count
    double precision :: incident_energy_mev
    double precision :: theta_degrees
    double precision :: phi_degrees
    character(len=200) :: header_line

    open(newunit=unit_number, file=trim(scattCondFile), status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, 'ERROR: could not open scattering conditions file: ', trim(scattCondFile)
      stop
    end if

    read(unit_number,'(A)',iostat=io_status) header_line
    if (io_status /= 0) stop 'ERROR: failed to read scattering conditions header line.'
    read(unit_number,*,iostat=io_status) data%helium_mass
    if (io_status /= 0) stop 'ERROR: failed to read helium mass from scattering conditions file.'

    print *, 'Helium mass = ', data%helium_mass

    condition_count = 0
    do
      read(unit_number,*,iostat=io_status) incident_energy_mev, theta_degrees, phi_degrees
      if (io_status < 0) exit
      if (io_status > 0) then
        print *, 'ERROR: invalid line in scattering conditions file (empty/malformed line).'
        stop
      end if
      condition_count = condition_count + 1
    end do

    data%condition_count = condition_count
    allocate(data%incident_energy_mev(condition_count), data%theta_degrees(condition_count), data%phi_degrees(condition_count))

    rewind(unit_number)
    read(unit_number,'(A)') header_line
    read(unit_number,*) data%helium_mass
    do i = 1, condition_count
      read(unit_number,*) data%incident_energy_mev(i), data%theta_degrees(i), data%phi_degrees(i)
    end do

    close(unit_number)
  end subroutine load_scattering_conditions_file

  subroutine load_fixed_potential(fourierfile, rmlmda, data)
    implicit none
    character(len=*), intent(in) :: fourierfile
    double precision, intent(in) :: rmlmda
    type(FixedPotentialData), intent(out) :: data
    integer :: i
    integer :: j
    integer :: unit_number
    integer :: io_status
    integer :: fourier_x_index
    integer :: fourier_y_index
    integer :: fourier_component_count_from_file
    integer :: z_point_count_from_header
    character(len=200) :: header_line
    integer, parameter :: MAX_Z_FIXED = 550
    integer, parameter :: MAX_FIXED_FOURIER_COMPONENTS = 4096
    integer, parameter :: MAX_FOURIER_COMPONENTS = 4096

    open(newunit=unit_number, file=trim(fourierfile), status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, 'ERROR: could not open potential file: ', trim(fourierfile)
      stop
    end if

    read(unit_number,'(A)') header_line
    read(unit_number,'(A)') header_line
    read(unit_number,*) fourier_component_count_from_file, data%fourier_grid_x_count, &
      data%fourier_grid_y_count, z_point_count_from_header
    read(unit_number,'(A)') header_line
    read(unit_number,*) data%unit_cell_ax, data%unit_cell_ay, data%unit_cell_bx, data%unit_cell_by
    read(unit_number,'(A)') header_line
    read(unit_number,*) data%zmin, data%zmax
    read(unit_number,'(A)') header_line
    read(unit_number,'(A)') header_line
    read(unit_number,'(A)') header_line

    data%fourier_component_count = fourier_component_count_from_file
    data%z_point_count = z_point_count_from_header

    if (data%fourier_component_count /= (data%fourier_grid_x_count * data%fourier_grid_y_count)) then
      print *, 'ERROR: inconsistent potential header; nfc must equal nkx*nky.'
      stop
    end if

      if (data%fourier_component_count > MAX_FOURIER_COMPONENTS) then
      print *, 'ERROR: the potential file needs more fourier components', &
      ' than allowed by the .inc file (nfc>nfcx)'
      stop
      else if (data%z_point_count > MAX_Z_FIXED) then
      print *, 'ERROR: the potential file needs more z points than', &
      ' allowed by the .inc file (nzfixed>NZFIXED_MAX)'
      stop
      else if (data%fourier_component_count > MAX_FIXED_FOURIER_COMPONENTS) then
      print *, 'ERROR: the potential file needs more fourier components than', &
      ' allowed by the .inc file (nfc>NVFCFIXED_MAX)'
      stop
    end if

    allocate(data%fourier_indices_x(data%fourier_component_count), &
      data%fourier_indices_y(data%fourier_component_count))
    allocate(data%fixed_fourier_values(data%z_point_count, data%fourier_component_count))
    data%fixed_fourier_values = (0.0d0, 0.0d0)

    data%specular_component_index = 1
    i = 0
    do fourier_x_index = 0, data%fourier_grid_x_count - 1
      do fourier_y_index = 0, data%fourier_grid_y_count - 1
        i = i + 1
        data%fourier_indices_x(i) = fourier_x_index
        if (fourier_x_index > ((data%fourier_grid_x_count - 1) / 2)) then
          data%fourier_indices_x(i) = fourier_x_index - data%fourier_grid_x_count
        end if
        data%fourier_indices_y(i) = fourier_y_index
        if (fourier_y_index > ((data%fourier_grid_y_count - 1) / 2)) then
          data%fourier_indices_y(i) = fourier_y_index - data%fourier_grid_y_count
        end if
        if ((data%fourier_indices_x(i) == 0) .and. (data%fourier_indices_y(i) == 0)) then
          data%specular_component_index = i
        end if
      end do
    end do

    do i = 1, data%fourier_component_count
      do j = 1, data%z_point_count
        read(unit_number,*) data%fixed_fourier_values(j,i)
      end do
    end do
    close(unit_number)

    data%fixed_fourier_values = data%fixed_fourier_values * rmlmda

    print *, 'Total number of fourier components from potential = ', data%fourier_component_count
    print *, 'Number of z points in fourier components (nzfixed) = ', data%z_point_count
    print *, 'Unit cell vectors (A): a = (', data%unit_cell_ax, ',', data%unit_cell_ay, &
      '), b = (', data%unit_cell_bx, ',', data%unit_cell_by, ')'
    print *, 'z integration range = (', data%zmin, ',', data%zmax, ')'
  end subroutine load_fixed_potential

end module multiscat_io_loaders