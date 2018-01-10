! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_box_model program

!> A simple box model for testing \ref phlex_chem "phlex-chem"
!> mechanisms.
program pmc_phlex_box_model

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind), parameter :: NUM_TIME_STEP = 100
  ! Time step (s)
  real(kind=dp), parameter :: TIME_STEP = 0.1
  ! Temperature (K)
  real(kind=dp), parameter :: ENV_TEMP = 298.0

  ! Config file list path
  character(len=300) :: file_list_arg
  character(len=:), allocatable :: file_list
  ! Output file name
  character(len=300) :: output_file_arg
  character(len=:), allocatable :: output_file

  ! Initialize mpi
  call pmc_mpi_init()

  if (pmc_mpi_rank() == 0) then
    if (command_argument_count() /= 2) then
      write(*,*) "Usage: phlex_box_model input_file_list.jsoni output_file.txt"
      call die_msg(695622653, "Incorrect number of command line arguments")
    end if
    call get_command_argument(1, file_list_arg)
    call get_command_argument(2, output_file_arg)
  end if

  file_list = trim(file_list_arg)
  output_file = trim(output_file_arg)

  ! Run the phlex-chem box model
  call run_phlex_box(file_list, output_file)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the phlex-chem box model
  subroutine run_phlex_box(file_list, output_file)

    use pmc_integration_data

    !> Input file list
    character(len=:), allocatable :: file_list
    !> Output file name
    character(len=:), allocatable :: output_file

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    type(integration_data_t), pointer :: integration_data

    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: buffer_size, pos
 
    type(string_t), allocatable :: species_names(:)
    type(property_t), pointer :: prop_set
    character(len=:), allocatable :: key
    real(kind=dp) :: real_val

    character(len=:), allocatable :: file_row
    real(kind=dp), allocatable :: model_conc(:,:)
    integer(kind=i_kind) :: i_time, i_spec, row_size, row_pos, row_end

    real(kind=dp) :: comp_start, comp_end

    ! Check for an available solver
    integration_data => integration_data_t()
    call assert_msg(831452409, integration_data%is_solver_available(), &
            "No solver available")

    ! Set up the phlex_core_t and phlex_state_t variables
    if (pmc_mpi_rank() == 0) then

      ! Start the computational timer
      call cpu_time(comp_start)

      ! Construct a phlex_core_t variable and load input files
      phlex_core => phlex_core_t(file_list)

      ! Initialize the phlex_core_t variable
      call phlex_core%initialize()

      ! Get a new state variable
      phlex_state = phlex_core%new_state()

      ! Set the environmental conditions
      phlex_state%env_state%temp = ENV_TEMP

      ! Set the initial conditions
      phlex_state%state_var(:) = real(0.0, kind=dp)
      key = "init conc"
      species_names = phlex_core%chem_spec_data%spec_names_by_type(GAS_SPEC)
      do i_spec = 1, size(species_names)
        prop_set = phlex_core%chem_spec_data%get_property_set( &
                species_names(i_spec)%string)   
        if (associated(prop_set)) then
          if (prop_set%get_real(key, real_val)) then 
            phlex_state%state_var( &
              phlex_core%chem_spec_data%gas_state_id( &
                  species_names(i_spec)%string)) = real_val
          end if
        end if
      end do

#ifdef PMC_USE_MPI
      ! Get the size of the mpi buffer
      buffer_size = phlex_core%pack_size() + &
              phlex_state%pack_size()
    end if

    ! Tell everyone the buffer size
    call pmc_mpi_bcast_integer(buffer_size)
    allocate(buffer(buffer_size))
    pos = 0
      
    if (pmc_mpi_rank() == 0) then
      ! Pack the initialized core and state to the other nodes
      call phlex_core%bin_pack(buffer, pos)
      call phlex_state%bin_pack(buffer, pos)
    end if

    ! Broadcast the data to everyone
    call pmc_mpi_bcast_packed(buffer)
    
    if (pmc_mpi_rank() /= 0) then

      ! Construct an empty phlex_core_t variable
      phlex_core => phlex_core_t()

      ! Unpack the data
      call phlex_core%bin_unpack(buffer, pos)
      call phlex_state%bin_unpack(buffer, pos)
#endif
    end if

    ! Set up the model concentration array
    allocate(model_conc(0:NUM_TIME_STEP,size(phlex_state%state_var)))
    model_conc(0,:) = phlex_state%state_var(:)

    ! Calculate the initialization time, and reset the timer for the model run
    if (pmc_mpi_rank() == 0) then
      call cpu_time(comp_end)
      write(*,*) "Initialization time: ", comp_end-comp_start, " s"
      call cpu_time(comp_start)
    end if

    write (*,*) "***** Phlex-chem configuration *******"
    call phlex_core%mechanism(1)%print()
    write (*,*) "Phlex state: ", phlex_state%state_var(:)
    write (*,*) "Temp: ", phlex_state%env_state%temp
    write (*,*) "***** end Phlex-chem configuration *******"

    ! Integrate the mechanism
    do i_time = 1, NUM_TIME_STEP

      ! Get the modeled conc
      call phlex_core%solve(phlex_state, TIME_STEP)

      ! Save the modeled conc
      model_conc(i_time,:) = phlex_state%state_var(:)

    end do

    ! Calculate the model run time
    if (pmc_mpi_rank() == 0) then
      call cpu_time(comp_end)
      write(*,*) "Model run time: ", comp_end-comp_start, " s"
    end if

    ! TODO Check child nodes to be sure results match parent node

    ! Save the results
    if (pmc_mpi_rank() == 0) then
      open(unit=7, file="out/"//output_file, status="replace", action="write")
      species_names = phlex_core%chem_spec_data%spec_names_by_type(GAS_SPEC)
      row_size = 0
      do i_spec = 1, size(species_names)
        row_size = row_size + len(species_names(i_spec)%string) + 1
      end do
      allocate(character(row_size) :: file_row)
      row_pos = 1
      do i_spec = 1, size(species_names)
        row_end = row_pos + len(species_names(i_spec)%string)
        file_row(row_pos:row_end) = species_names(i_spec)%string//" "
        row_pos = row_end + 1
      end do
      write (7,*) "time ", file_row
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, model_conc(i_time,:)
      end do
    end if

  end subroutine run_phlex_box 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_phlex_box_model
