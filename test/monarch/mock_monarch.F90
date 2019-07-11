! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The mock_monarch program

!> Mock version of the MONARCH model for testing integration with PartMC
program mock_monarch

  use pmc_util,                                 only : assert_msg
  use pmc_monarch_interface
  use pmc_mpi

  implicit none

  !> File unit for model run-time messages
  integer, parameter :: OUTPUT_FILE_UNIT = 6
  !> File unit for model results
  integer, parameter :: RESULTS_FILE_UNIT = 7
  !> File unit for script generation
  integer, parameter :: SCRIPTS_FILE_UNIT = 8

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters for mock MONARCH model !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of total species in mock MONARCH
  integer, parameter :: NUM_MONARCH_SPEC = 800
  !> Number of W-E cells in mock MONARCH
  integer, parameter :: NUM_WE_CELLS = 20
  !> Number of S-N cells in mock MONARCH
  integer, parameter :: NUM_SN_CELLS = 30
  !> Number of vertical cells in mock MONARCH
  integer, parameter :: NUM_VERT_CELLS = 1
  !> Starting W-E cell for phlex-chem call
  integer, parameter :: I_W = 1!9
  !> Ending W-E cell for phlex-chem call
  integer, parameter :: I_E = 15!15!11
  !> Starting S-N cell for phlex-chem call
  integer, parameter :: I_S = 1!14
  !> Ending S-N cell for phlex-chem call
  integer, parameter :: I_N = 15!15!16
  !> Starting index for phlex-chem species in tracer array
  integer, parameter :: START_PHLEX_ID = 100
  !> Ending index for phlex-chem species in tracer array
  integer, parameter :: END_PHLEX_ID = 650
  !> Time step (min)
  real, parameter :: TIME_STEP = 1.6
  !> Number of time steps to integrate over
  integer, parameter :: NUM_TIME_STEP = 20!100
  !> Index for water vapor in water_conc()
  integer, parameter :: WATER_VAPOR_ID = 5
  !> Start time
  real, parameter :: START_TIME = 360.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! State variables for mock MONARCH model !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> NMMB style arrays (W->E, S->N, top->bottom, ...)
  !> Temperature (K)
  real :: temperature(NUM_WE_CELLS, NUM_SN_CELLS, NUM_VERT_CELLS)
  !> Species conc (various units)
  real :: species_conc(NUM_WE_CELLS, NUM_SN_CELLS, NUM_VERT_CELLS, NUM_MONARCH_SPEC)
  !> Water concentrations (kg_H2O/kg_air)
  real :: water_conc(NUM_WE_CELLS, NUM_SN_CELLS, NUM_VERT_CELLS, WATER_VAPOR_ID)

  !> WRF-style arrays (W->E, bottom->top, N->S)
  !> Air density (kg_air/m^3)
  real :: air_density(NUM_WE_CELLS, NUM_VERT_CELLS, NUM_SN_CELLS)
  !> Air pressure (Pa)
  real :: pressure(NUM_WE_CELLS, NUM_VERT_CELLS, NUM_SN_CELLS)

  !> Starting time for mock model run (min since midnight) TODO check how time
  !! is tracked in MONARCH
  real :: curr_time = START_TIME

  !> Plot start time is after first call to solve chemistry, so initial
  !! concentrations do not affect y-axis scaling
  real :: plot_start_time = START_TIME + TIME_STEP

  !> !!! Add to MONARCH variables !!!
  type(monarch_interface_t), pointer :: pmc_interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Mock model setup and evaluation variables !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Phlex-chem input file file
  character(len=:), allocatable :: phlex_input_file
  !> PartMC-phlex <-> MONARCH interface configuration file
  character(len=:), allocatable :: interface_input_file
  !> Results file prefix
  character(len=:), allocatable :: output_file_prefix

  character(len=500) :: arg
  integer :: status_code, i_time
  integer :: num_cells = 1


  ! Check the command line arguments
  call assert_msg(129432506, command_argument_count().eq.3, "Usage: "// &
          "./mock_monarch phlex_input_file_list.json "// &
          "interface_input_file.json output_file_prefix")

  ! initialize mpi (to take the place of a similar MONARCH call)
  call pmc_mpi_init()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! **** Add to MONARCH during initialization **** !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize PartMC-phlex
  call get_command_argument(1, arg, status=status_code)
  call assert_msg(678165802, status_code.eq.0, "Error getting PartMC-phlex "//&
          "configuration file name")
  phlex_input_file = trim(arg)
  call get_command_argument(2, arg, status=status_code)
  call assert_msg(664104564, status_code.eq.0, "Error getting PartMC-phlex "//&
          "<-> MONARCH interface configuration file name")
  interface_input_file = trim(arg)

  !Cells to solve simultaneously
  num_cells = (I_E - I_W+1)*(I_N - I_S+1)*NUM_VERT_CELLS
  !num_cells = 1

  pmc_interface => monarch_interface_t(phlex_input_file, interface_input_file, &
          START_PHLEX_ID, END_PHLEX_ID, num_cells)
  !pmc_interface => monarch_interface_t(phlex_input_file, interface_input_file, &
  !        START_PHLEX_ID, END_PHLEX_ID)
  deallocate(phlex_input_file)
  deallocate(interface_input_file)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! **** end initialization modification **** !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize the mock model
  call get_command_argument(3, arg, status=status_code)
  call assert_msg(234156729, status_code.eq.0, "Error getting output file prefix")
  output_file_prefix = trim(arg)
  call model_initialize(output_file_prefix)
  call pmc_interface%get_init_conc(species_conc, water_conc, WATER_VAPOR_ID, &
          air_density)

  ! call pmc_interface%print( )

  ! Run the model
  do i_time=0, NUM_TIME_STEP

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! **** Add to MONARCH during runtime for each time step **** !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call output_results(curr_time)
    call pmc_interface%integrate(curr_time,         & ! Starting time (min)
                                 TIME_STEP,         & ! Time step (min)
                                 I_W,               & ! Starting W->E grid cell
                                 I_E,               & ! Ending W->E grid cell
                                 I_S,               & ! Starting S->N grid cell
                                 I_N,               & ! Ending S->N grid cell
                                 temperature,       & ! Temperature (K)
                                 species_conc,      & ! Tracer array
                                 water_conc,        & ! Water concentrations (kg_H2O/kg_air)
                                 WATER_VAPOR_ID,    & ! Index in water_conc() corresponding to water vapor
                                 air_density,       & ! Air density (kg_air/m^3)
                                 pressure)            ! Air pressure (Pa)
    curr_time = curr_time + TIME_STEP

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! **** end runtime modification **** !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end do

  write(*,*) "Model run time: ", comp_time, " s"

  ! Output results and scripts
  if (pmc_mpi_rank().eq.0) then
    call output_results(curr_time)
    call create_gnuplot_script(pmc_interface, output_file_prefix, &
            plot_start_time, curr_time)
  end if

  ! TODO evaluate results

  write(*,*) "MONARCH interface tests - PASS"

  ! close the output file
  close(RESULTS_FILE_UNIT)
  deallocate(output_file_prefix)

  ! free the interface
  deallocate(pmc_interface)

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the mock model
  subroutine model_initialize(file_prefix)

    !> File prefix for model results
    character(len=:), allocatable, intent(in) :: file_prefix

    integer :: i_spec
    character(len=:), allocatable :: file_name

    file_name = file_prefix//"_results.txt"

    ! Open the output file
    open(RESULTS_FILE_UNIT, file=file_name, status="replace", action="write")

    ! TODO refine initial model conditions
    temperature(:,:,:) = 300.614166259766
    species_conc(:,:,:,:) = 0.0
    water_conc(:,:,:,:) = 0.0
    water_conc(:,:,:,WATER_VAPOR_ID) = 0.01
    air_density(:,:,:) = 1.225
    pressure(:,:,:) = 94165.7187500000

    deallocate(file_name)

  end subroutine model_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the model results
  subroutine output_results(curr_time)

    !> Current model time (min since midnight)
    real, intent(in) :: curr_time

    write(RESULTS_FILE_UNIT, *) curr_time, &
            species_conc(2,2,1,START_PHLEX_ID:END_PHLEX_ID), &
            water_conc(2,2,1,WATER_VAPOR_ID)
    !            species_conc(10,15,1,START_PHLEX_ID:END_PHLEX_ID), &
    !        water_conc(10,15,1,WATER_VAPOR_ID)

  end subroutine output_results
!Tolerance -E4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a gnuplot script for viewing species concentrations
  subroutine create_gnuplot_script(pmc_interface, file_prefix, start_time, &
            end_time)

    !> PartMC-phlex <-> MONARCH interface
    type(monarch_interface_t), intent(in) :: pmc_interface
    !> File prefix for gnuplot script
    character(len=:), allocatable :: file_prefix
    !> Plot start time
    real :: start_time
    !> Plot end time
    real :: end_time

    type(string_t), allocatable :: species_names(:)
    integer(kind=i_kind), allocatable :: tracer_ids(:)
    character(len=:), allocatable :: file_name, spec_name
    integer(kind=i_kind) :: i_char, i_spec, tracer_id

    ! Get the species names and ids
    call pmc_interface%get_MONARCH_species(species_names, tracer_ids)

    ! Adjust the tracer ids to match the results file
    tracer_ids(:) = tracer_ids(:) - START_PHLEX_ID + 2

    ! Create the gnuplot script
    file_name = file_prefix//".conf"
    open(unit=SCRIPTS_FILE_UNIT, file=file_name, status="replace", action="write")
    write(SCRIPTS_FILE_UNIT,*) "# "//file_name
    write(SCRIPTS_FILE_UNIT,*) "# Run as: gnuplot "//file_name
    write(SCRIPTS_FILE_UNIT,*) "set terminal png truecolor"
    write(SCRIPTS_FILE_UNIT,*) "set autoscale"
    write(SCRIPTS_FILE_UNIT,*) "set xrange [", start_time, ":", end_time, "]"
    do i_spec = 1, size(species_names)
      spec_name = species_names(i_spec)%string
      forall (i_char = 1:len(spec_name), spec_name(i_char:i_char).eq.'/') &
                spec_name(i_char:i_char) = '_'
      write(SCRIPTS_FILE_UNIT,*) "set output '"//file_prefix//"_"// &
              spec_name//".png'"
      write(SCRIPTS_FILE_UNIT,*) "plot\"
      write(SCRIPTS_FILE_UNIT,*) " '"//file_prefix//"_results.txt'\"
      write(SCRIPTS_FILE_UNIT,*) " using 1:"// &
              trim(to_string(tracer_ids(i_spec)))//" title '"// &
              species_names(i_spec)%string//" (MONARCH)'"
    end do
    tracer_id = END_PHLEX_ID - START_PHLEX_ID + 3
    write(SCRIPTS_FILE_UNIT,*) "set output '"//file_prefix//"_H2O.png'"
    write(SCRIPTS_FILE_UNIT,*) "plot\"
    write(SCRIPTS_FILE_UNIT,*) " '"//file_prefix//"_results.txt'\"
    write(SCRIPTS_FILE_UNIT,*) " using 1:"// &
            trim(to_string(tracer_id))//" title 'H2O (MONARCH)'"
    close(SCRIPTS_FILE_UNIT)

    deallocate(species_names)
    deallocate(tracer_ids)
    deallocate(file_name)
    deallocate(spec_name)

  end subroutine create_gnuplot_script

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program mock_monarch
