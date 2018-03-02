! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_monarch module

!> Interface for the MONACH model and PartMC-phlex
module pmc_monarch

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : assert_msg, string_t
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_property
  use pmc_phlex_solver_data
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: pmc_initialize, pmc_integrate

  !> MPI node id from MONARCH
  integer(kind=i_kind) :: MONARCH_NODE = 0 ! TODO fill this in with MONARCH param
  !> Phlex-chem core
  type(phlex_core_t), pointer, save :: phlex_core 
  !> Phlex-chem state
  type(phlex_state_t), target, save :: phlex_state
  !> MONACH <-> PartMC species map
  integer(kind=i_kind), allocatable, save :: spec_map(:)
  !> Starting index for PartMC species on MONARCH tracter array
  integer(kind=i_kind) :: tracer_starting_id

  !> List of species names to map between MONARCH and PartMC.
  !! These should be in same order as they appear in the MONARCH
  !! tracer array.
  !! TODO Could read this from an input file in the future
  type(string_t), allocatable :: spec_names(:)
  !! TODO Set up similar mapping for aerosol species

  ! TEMPORARY
  real(kind=dp), public, save :: comp_time = 0.0d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the model. Call this routine once at the beginning of the 
  !! model run from each node. The master node should pass a string containing
  !! the path to the PartMC input file directory, and the starting and
  !! ending indices for chemical species in the tracer array.
  subroutine pmc_initialize(input_file_path, starting_id, ending_id)

    !> Path to the PartMC input file directory
    character(len=:), allocatable, optional :: input_file_path
    !> Starting index for chemical species in the MONARCH tracer array
    integer, optional :: starting_id
    !> Ending index for chemical species in the MONARCH tracer array
    integer, optional :: ending_id

    type(phlex_solver_data_t), pointer :: phlex_solver_data
    character(len=:), allocatable :: buffer
    integer(kind=i_kind) :: pos, pack_size
    integer(kind=i_kind) :: i_spec

    ! Computation time variable
    real(kind=dp) :: comp_start, comp_end

    ! Check for an available solver
    call assert_msg(332298164, phlex_solver_data%is_solver_available(), &
            "No solver available")

    ! Initialize the time-invariant model data on each node
    if (MONARCH_NODE.eq.0) then
      
      ! Start the computation timer on the primary node
      call cpu_time(comp_start)

      call assert_msg(304676624, present(input_file_path), &
              "Missing input file directory")
      call assert_msg(937567597, present(starting_id), &
              "Missing starting tracer index for chemical species")
      call assert_msg(593895016, present(ending_id), &
              "Missing ending tracer index for chemical species")
      call assert_msg(478217898, ending_id-starting_id+1 .ge. size(spec_names), &
              "Not enough room for chemical species on the tracer array")
      
      ! Initialize the phlex-chem core
      phlex_core => phlex_core_t(input_file_path)
      call phlex_core%initialize()

      ! Generate the MONARCH <-> PartMC species index map
      tracer_starting_id = starting_id
      allocate(spec_names(3))
      spec_names(1) = string_t("A")
      spec_names(2) = string_t("B")
      spec_names(3) = string_t("C")
      allocate(spec_map(size(spec_names)))
      do i_spec = 1, size(spec_names)
        spec_map(i_spec) = &
              phlex_core%chem_spec_data%gas_state_id(spec_names(i_spec)%string)
        call assert_msg(916977002, spec_map(i_spec).gt.0, &
                "Could not find index for species '"//spec_names(i_spec)%string)
      end do

#ifdef PMC_USE_MPI
      pack_size = phlex_core%pack_size() + &
              pmc_mpi_pack_size_string_array(spec_names)
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
      call pmc_mpi_pack_string_array(spec_names)
      ! TODO send buffer to all the other nodes
    else
      phlex_core => phlex_core_t()
      ! TODO get buffer from the primary node
      pos = 0
      call phlex_core%bin_unpack(buffer, pos)
      call pmc_mpi_unpack_string_array(spec_names)
#endif
    end if

    ! Create a state variable on each node
    phlex_state = phlex_core%new_state()

    ! Calculate the intialization time
    if (MONARCH_NODE.eq.0) then
      call cpu_time(comp_end)
      write(*,*) "Initialization time: ", comp_end-comp_start, " s"
    end if

  end subroutine pmc_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the PartMC mechanism for a particular set of cells and timestep
  subroutine pmc_integrate(start_time, time_step, i_start, i_end, j_start, &
                  j_end, temperature, MONARCH_conc, water_conc, &
                  water_vapor_index, air_density, pressure)

    !> Integration start time (min since midnight)
    real, intent(in) :: start_time
    !> Integration time step
    real, intent(in) :: time_step
    !> Grid-cell W->E starting index
    integer, intent(in) :: i_start
    !> Grid-cell W->E ending index
    integer, intent(in) :: i_end
    !> Grid-cell S->N starting index
    integer, intent(in) :: j_start
    !> Grid-cell S->N ending index
    integer, intent(in) :: j_end

    !> NMMB style arrays (W->E, S->N, top->bottom, ...)
    !> Temperature (K)
    real, intent(in) :: temperature(:,:,:)
    !> MONARCH species concentration (ppm or ug/m^3)
    real, intent(inout) :: MONARCH_conc(:,:,:,:)
    !> Atmospheric water concentrations (kg_H2O/kg_air)
    real, intent(in) :: water_conc(:,:,:,:)
    !> Index in water_conc corresponding to water vapor
    integer, intent(in) :: water_vapor_index

    !> WRF-style arrays (W->E, bottom->top, N->S)
    !> Air density (kg_air/m^3)
    real, intent(in) :: air_density(:,:,:)
    !> Pressure (Pa)
    real, intent(in) :: pressure(:,:,:)

    integer :: i, j, k, k_flip, i_spec

    ! Computation time variables
    real(kind=dp) :: comp_start, comp_end

    ! Loop through the grid cells
    do i=i_start, i_end
      do j=j_start, j_end
        do k=1, size(MONARCH_conc,3)

          ! Calculate the vertical index for NMMB-style arrays
          k_flip = size(MONARCH_conc,3) - k + 1

          ! Update the environmental state
          phlex_state%env_state%temp = temperature(i,j,k_flip)
          phlex_state%env_state%pressure = pressure(i,k,j)
          ! TODO finish environmental state setup

          ! Update species concentrations in PMC
          phlex_state%state_var(:) = 0.0
          phlex_state%state_var(spec_map(:)) = &
                  MONARCH_conc(i,j,k_flip,tracer_starting_id:&
                  tracer_starting_id+size(spec_map))

          ! Start the computation timer
          if (MONARCH_NODE.eq.0 .and. i.eq.i_start .and. j.eq.j_start &
                  .and. k.eq.1) then
            call cpu_time(comp_start)
          end if

          ! Integrate the PMC mechanism
          call phlex_core%solve(phlex_state, real(time_step, kind=dp))

          ! Calculate the computation time
          if (MONARCH_NODE.eq.0 .and. i.eq.i_start .and. j.eq.j_start &
                  .and. k.eq.1) then
            call cpu_time(comp_end)
            comp_time = comp_time + (comp_end-comp_start)
          end if

          ! Update the MONARCH tracer array with new species concentrations
          MONARCH_conc(i,j,k_flip,tracer_starting_id:&
                  tracer_starting_id+size(spec_map)) = &
                  phlex_state%state_var(spec_map(:))
 
        end do
      end do
    end do

  end subroutine pmc_integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_monarch
