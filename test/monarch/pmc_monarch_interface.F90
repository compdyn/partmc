! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The monarch_interface_t object and related functions

!> Interface for the MONACH model and PartMC-camp
module pmc_monarch_interface

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : assert_msg, string_t, &
                                             warn_assert_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_modal_binned_mass
  use pmc_chem_spec_data
  use pmc_property
  use pmc_camp_solver_data
  use pmc_mechanism_data,            only : mechanism_data_t
  use pmc_rxn_data,                  only : rxn_data_t
  use pmc_rxn_photolysis
  use pmc_solver_stats
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: monarch_interface_t

  !> PartMC <-> MONARCH interface
  !!
  !! Contains all data required to intialize and run PartMC from MONARCH data
  !! and map state variables between PartMC and MONARCH
  type :: monarch_interface_t
    !private
    !> CAMP-chem core
    type(camp_core_t), pointer :: camp_core
    !> CAMP-chem state
    type(camp_state_t), pointer :: camp_state
    !> MONARCH species names
    type(string_t), allocatable :: monarch_species_names(:)
    !> MONARCH <-> PartMC species map
    integer(kind=i_kind), allocatable :: map_monarch_id(:), map_camp_id(:)
    !> PartMC-camp ids for initial concentrations
    integer(kind=i_kind), allocatable :: init_conc_camp_id(:)
    !> Initial species concentrations
    real(kind=dp), allocatable :: init_conc(:)
    !> Number of cells to compute simultaneously
    integer(kind=i_kind) :: n_cells = 1
    !> Starting index for PartMC species on the MONARCH tracer array
    integer(kind=i_kind) :: tracer_starting_id
    !> Ending index for PartMC species on the MONARCH tracer array
    integer(kind=i_kind) :: tracer_ending_id
    !> PartMC-camp <-> MONARCH species map input data
    type(property_t), pointer :: species_map_data
    !> Gas-phase water id in PartMC-camp
    integer(kind=i_kind) :: gas_phase_water_id
    !> Initial concentration data
    type(property_t), pointer :: init_conc_data
    !> Interface input data
    type(property_t), pointer :: property_set
    type(rxn_update_data_photolysis_t), allocatable :: rate_update(:)
    real(kind=dp), allocatable :: photo_rates(:)
    integer :: n_photo_rxn
    !> Solve multiple grid cells at once?
    logical :: solve_multiple_cells = .false.
  contains
    !> Integrate PartMC for the current MONARCH state over a specified time step
    procedure :: integrate
    !> Get initial concentrations (for testing only)
    procedure :: get_init_conc
    !> Get monarch species names and ids (for testing only)
    procedure :: get_MONARCH_species
    !> Print the PartMC-camp data
    procedure :: print => do_print
    !> Load interface data from a set of input files
    procedure, private :: load
    !> Create the PartMC <-> MONARCH species map
    procedure, private :: create_map
    !> Load the initial concentrations
    procedure, private :: load_init_conc
    !> Finalize the interface
    final :: finalize
  end type monarch_interface_t

  !> PartMC <-> MONARCH interface constructor
  interface monarch_interface_t
    procedure :: constructor
  end interface monarch_interface_t

  !> MPI node id from MONARCH
  integer(kind=i_kind) :: MONARCH_PROCESS ! TODO replace with MONARCH param
  ! TEMPORARY
  real(kind=dp), public, save :: comp_time = 0.0d0

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create and initialize a new monarch_interface_t object
  !!
  !! Create a monarch_interface_t object at the beginning of the  model run
  !! for each node. The master node should pass a string containing the path
  !! to the PartMC confirguration file list, the path to the interface
  !! configuration file and the starting and ending indices for chemical
  !! species in the tracer array.
  function constructor(camp_config_file, interface_config_file, &
                       starting_id, ending_id, n_cells, mpi_comm) result (this)

    !> A new MONARCH interface
    type(monarch_interface_t), pointer :: this
    !> Path to the PartMC-camp configuration file list
    character(len=:), allocatable, optional :: camp_config_file
    !> Path to the PartMC-camp <-> MONARCH interface input file
    character(len=:), allocatable, optional :: interface_config_file
    !> Starting index for chemical species in the MONARCH tracer array
    integer, optional :: starting_id
    !> Ending index for chemical species in the MONARCH tracer array
    integer, optional :: ending_id
    !> MPI communicator
    integer, intent(in), optional :: mpi_comm
    !> Num cells to compute simulatenously
    integer, optional :: n_cells

    type(camp_solver_data_t), pointer :: camp_solver_data
    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: pos, pack_size
    integer(kind=i_kind) :: i_spec, i_photo_rxn
    type(string_t), allocatable :: unique_names(:)

    class(aero_rep_data_t), pointer :: aero_rep
    integer(kind=i_kind) :: i_sect_om, i_sect_bc, i_sect_sulf, i_sect_opm
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    ! Computation time variable
    real(kind=dp) :: comp_start, comp_end

#ifdef PMC_USE_MPI
    integer :: local_comm

    if (present(mpi_comm)) then
      local_comm = mpi_comm
    else
      local_comm = MPI_COMM_WORLD
    endif
#endif
    ! Set the MPI rank (TODO replace with MONARCH param)
    MONARCH_PROCESS = pmc_mpi_rank()

    ! Create a new interface object
    allocate(this)

    if (.not.present(n_cells).or.n_cells.eq.1) then
      this%solve_multiple_cells = .false.
    else
      this%solve_multiple_cells = .true.
      this%n_cells=n_cells
    end if

    ! Check for an available solver
    camp_solver_data => camp_solver_data_t()
    call assert_msg(332298164, camp_solver_data%is_solver_available(), &
            "No solver available")
    deallocate(camp_solver_data)

    ! Initialize the time-invariant model data on each node
    if (MONARCH_PROCESS.eq.0) then

      ! Start the computation timer on the primary node
      call cpu_time(comp_start)

      call assert_msg(304676624, present(camp_config_file), &
              "Missing PartMC-camp configuration file list")
      call assert_msg(194027509, present(interface_config_file), &
              "Missing MartMC-camp <-> MONARCH interface configuration file")
      call assert_msg(937567597, present(starting_id), &
              "Missing starting tracer index for chemical species")
      call assert_msg(593895016, present(ending_id), &
              "Missing ending tracer index for chemical species")

      ! Load the interface data
      call this%load(interface_config_file)

      ! Initialize the camp-chem core
      this%camp_core => camp_core_t(camp_config_file, this%n_cells)
      call this%camp_core%initialize()

      ! Set the aerosol representation id
      if (this%camp_core%get_aero_rep("MONARCH mass-based", aero_rep)) then
        select type (aero_rep)
          type is (aero_rep_modal_binned_mass_t)
            call this%camp_core%initialize_update_object( aero_rep, &
                                                             update_data_GMD)
            call this%camp_core%initialize_update_object( aero_rep, &
                                                             update_data_GSD)
            call assert(889473105, &
                        aero_rep%get_section_id("organic_matter", i_sect_om))
            call assert(648042550, &
                        aero_rep%get_section_id("black_carbon", i_sect_bc))
            !call assert(760360895, &
            !            aero_rep%get_section_id("sulfate", i_sect_sulf))
            call assert(307728742, &
                        aero_rep%get_section_id("other_PM", i_sect_opm))
          class default
            call die_msg(351392791, &
                         "Wrong type for aerosol representation "// &
                         "'MONARCH mass-based'")
        end select
      else
        i_sect_om = -1
        i_sect_bc = -1
        i_sect_sulf = -1
        i_sect_opm = -1
      end if

      ! Set the MONARCH tracer array bounds
      this%tracer_starting_id = starting_id
      this%tracer_ending_id = ending_id

      ! Generate the PartMC-camp <-> MONARCH species map
      call this%create_map()

      ! Load the initial concentrations
      call this%load_init_conc()

#ifdef PMC_USE_MPI

      ! Change a bit init_conc to denote different initial values
      !this%init_conc(:) = this%init_conc(:) + 0.1*MONARCH_PROCESS

      pack_size = this%camp_core%pack_size() + &
              update_data_GMD%pack_size() + &
              update_data_GSD%pack_size() + &
              pmc_mpi_pack_size_integer_array(this%map_monarch_id) + &
              pmc_mpi_pack_size_integer_array(this%map_camp_id) + &
              pmc_mpi_pack_size_integer_array(this%init_conc_camp_id) + &
              pmc_mpi_pack_size_real_array(this%init_conc) + &
              pmc_mpi_pack_size_integer(this%gas_phase_water_id) + &
              pmc_mpi_pack_size_integer(i_sect_om) + &
              pmc_mpi_pack_size_integer(i_sect_bc) + &
              pmc_mpi_pack_size_integer(i_sect_sulf) + &
              pmc_mpi_pack_size_integer(i_sect_opm)
      allocate(buffer(pack_size))
      pos = 0
      call this%camp_core%bin_pack(buffer, pos)
      call update_data_GMD%bin_pack(buffer, pos)
      call update_data_GSD%bin_pack(buffer, pos)
      call pmc_mpi_pack_integer_array(buffer, pos, this%map_monarch_id)
      call pmc_mpi_pack_integer_array(buffer, pos, this%map_camp_id)
      call pmc_mpi_pack_integer_array(buffer, pos, this%init_conc_camp_id)
      call pmc_mpi_pack_real_array(buffer, pos, this%init_conc)
      call pmc_mpi_pack_integer(buffer, pos, this%gas_phase_water_id)
      call pmc_mpi_pack_integer(buffer, pos, i_sect_om)
      call pmc_mpi_pack_integer(buffer, pos, i_sect_bc)
      call pmc_mpi_pack_integer(buffer, pos, i_sect_sulf)
      call pmc_mpi_pack_integer(buffer, pos, i_sect_opm)
    endif

    ! broadcast the buffer size
    call pmc_mpi_bcast_integer(pack_size, local_comm)

    if (MONARCH_PROCESS.ne.0) then
      ! allocate the buffer to receive data
      allocate(buffer(pack_size))
    end if

    ! boradcast the buffer
    call pmc_mpi_bcast_packed(buffer, local_comm)

    if (MONARCH_PROCESS.ne.0) then
      ! unpack the data
      this%camp_core => camp_core_t()
      pos = 0
      call this%camp_core%bin_unpack(buffer, pos)
      call update_data_GMD%bin_unpack(buffer, pos)
      call update_data_GSD%bin_unpack(buffer, pos)
      call pmc_mpi_unpack_integer_array(buffer, pos, this%map_monarch_id)
      call pmc_mpi_unpack_integer_array(buffer, pos, this%map_camp_id)
      call pmc_mpi_unpack_integer_array(buffer, pos, this%init_conc_camp_id)
      call pmc_mpi_unpack_real_array(buffer, pos, this%init_conc)
      call pmc_mpi_unpack_integer(buffer, pos, this%gas_phase_water_id)
      call pmc_mpi_unpack_integer(buffer, pos, i_sect_om)
      call pmc_mpi_unpack_integer(buffer, pos, i_sect_bc)
      call pmc_mpi_unpack_integer(buffer, pos, i_sect_sulf)
      call pmc_mpi_unpack_integer(buffer, pos, i_sect_opm)
#endif
    end if

#ifdef PMC_USE_MPI
    deallocate(buffer)
#endif

    ! Initialize the solver on all nodes
    call this%camp_core%solver_initialize()

    ! Create a state variable on each node
    this%camp_state => this%camp_core%new_state()

#ifndef ENABLE_CB05_SOA
    do i_photo_rxn = 1, this%n_photo_rxn
      call this%rate_update(i_photo_rxn)%set_rate(real(this%photo_rates(i_photo_rxn), kind=dp))
      call this%camp_core%update_data(this%rate_update(i_photo_rxn))
    end do
#endif

    ! Set the aerosol mode dimensions

    ! organic matter
    if (i_sect_om.gt.0) then
      call update_data_GMD%set_GMD(i_sect_om, 2.12d-8)
      call update_data_GSD%set_GSD(i_sect_om, 2.24d0)
      call this%camp_core%update_data(update_data_GMD)
      call this%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_bc.gt.0) then
    ! black carbon
      call update_data_GMD%set_GMD(i_sect_bc, 1.18d-8)
      call update_data_GSD%set_GSD(i_sect_bc, 2.00d0)
      call this%camp_core%update_data(update_data_GMD)
      call this%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_sulf.gt.0) then
    ! sulfate
      call update_data_GMD%set_GMD(i_sect_sulf, 6.95d-8)
      call update_data_GSD%set_GSD(i_sect_sulf, 2.12d0)
      call this%camp_core%update_data(update_data_GMD)
      call this%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_opm.gt.0) then
    ! other PM
      call update_data_GMD%set_GMD(i_sect_opm, 2.12d-8)
      call update_data_GSD%set_GSD(i_sect_opm, 2.24d0)
      call this%camp_core%update_data(update_data_GMD)
      call this%camp_core%update_data(update_data_GSD)
    end if

    ! Calculate the intialization time
    if (MONARCH_PROCESS.eq.0) then
      call cpu_time(comp_end)
      write(*,*) "Initialization time: ", comp_end-comp_start, " s"
      !call this%camp_core%print()
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the PartMC mechanism for a particular set of cells and timestep
  subroutine integrate(this, start_time, time_step, i_start, i_end, j_start, &
                  j_end, temperature, MONARCH_conc, water_conc, &
                  water_vapor_index, air_density, pressure, conv, i_hour)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this
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
    real, intent(in) :: conv
    integer, intent(inout) :: i_hour

#ifndef ENABLE_CB05_SOA

    type(chem_spec_data_t), pointer :: chem_spec_data

    real, dimension(30) :: SO2_emi
    real, dimension(30) :: NO2_emi
    real, dimension(30) :: NO_emi
    real, dimension(30) :: NH3_emi
    real, dimension(30) :: CO_emi
    real, dimension(30) :: ALD2_emi
    real, dimension(30) :: FORM_emi
    real, dimension(30) :: ETH_emi
    real, dimension(30) :: IOLE_emi
    real, dimension(30) :: OLE_emi
    real, dimension(30) :: TOL_emi
    real, dimension(30) :: XYL_emi
    real, dimension(30) :: PAR_emi
    real, dimension(30) :: ISOP_emi
    real, dimension(30) :: MEOH_emi
    real, dimension(30) :: rate_emi

    !todo Should be a better way than this
    integer :: SO2_id_camp
    integer :: NO2_id_camp
    integer :: NO_id_camp
    integer :: NH3_id_camp
    integer :: CO_id_camp
    integer :: ALD2_id_camp
    integer :: FORM_id_camp
    integer :: ETH_id_camp
    integer :: IOLE_id_camp
    integer :: OLE_id_camp
    integer :: TOL_id_camp
    integer :: XYL_id_camp
    integer :: PAR_id_camp
    integer :: ISOP_id_camp
    integer :: MEOH_id_camp

#endif

    ! MPI
    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: pos, pack_size
    integer :: local_comm
    real(kind=dp), allocatable :: mpi_conc(:)

    integer :: i, j, k, k_flip, i_spec, z, o, i2, i_cell, i_photo_rxn
    integer :: k_end

    ! Computation time variables
    real(kind=dp) :: comp_start, comp_end

    type(solver_stats_t), target :: solver_stats
    integer :: state_size_per_cell, n_cell_check

    if(this%n_cells.eq.1) then
      state_size_per_cell = 0
    else
      state_size_per_cell = this%camp_core%state_size_per_cell()
    end if

    k_end = size(MONARCH_conc,3)

#ifndef ENABLE_CB05_SOA
    call assert_msg(731700229, &
            this%camp_core%get_chem_spec_data(chem_spec_data), &
            "No chemical species data in camp_core.")

    SO2_emi = (/ 4.234E-09, 5.481E-09, 5.089E-09, 5.199E-09, 5.221E-09, 5.284E-09, 5.244E-09, 5.280E-09, 5.560E-09, 5.343E-09, 4.480E-09, 3.858E-09, 3.823E-09, 3.607E-09, 3.533E-09, 3.438E-09, 2.866E-09, 2.667E-09, 2.636E-09, 2.573E-09, 2.558E-09, 2.573E-09, 2.715E-09, 3.170E-09, 4.234E-09, 5.481E-09, 5.089E-09, 5.199E-09, 5.221E-09, 5.284E-09 /)
    NO2_emi = (/ 3.024E-09, 3.334E-09, 3.063E-09, 3.281E-09, 3.372E-09, 3.523E-09, 3.402E-09, 3.551E-09, 3.413E-09, 3.985E-09, 3.308E-09, 2.933E-09, 2.380E-09, 1.935E-09, 1.798E-09, 1.537E-09, 9.633E-10, 8.873E-10, 7.968E-10, 6.156E-10, 5.920E-10, 6.320E-10, 9.871E-10, 1.901E-09, .024E-09, 3.334E-09, 3.063E-09, 3.281E-09, 3.372E-09, 3.523E-09 /)
    NO_emi = (/ 5.749E-08, 6.338E-08, 5.825E-08, 6.237E-08, 6.411E-08, 6.699E-08, 6.468E-08, 6.753E-08, 6.488E-08, 7.575E-08, 6.291E-08, 5.576E-08, 4.524E-08, 3.679E-08, 3.419E-08, 2.924E-08, 1.832E-08, 1.687E-08, 1.515E-08, 1.171E-08, 1.125E-08, 1.202E-08, 1.877E-08, 3.615E-08, 5.749E-08, 6.338E-08, 5.825E-08, 6.237E-08, 6.411E-08, 6.699E-08 /)
    NH3_emi = (/ 8.93E-09, 8.705E-09, 1.639E-08, 1.466E-08, 1.6405E-08, 1.8805E-08, 1.65E-08, 1.8045E-08, 1.347E-08, 6.745E-09, 5.415E-09, 2.553E-09, 2.087E-09, 2.2885E-09, 2.7265E-09, 2.738E-09, 9.96E-10, 2.707E-09, 9.84E-10, 9.675E-10, 9.905E-10, 1.0345E-09, 1.0825E-09, 2.7465E-09, 8.93E-09, 8.705E-09, 1.639E-08, 1.466E-08, 1.6405E-08, 1.8805E-08 /)
    CO_emi = (/ 7.839E-07, 5.837E-07, 4.154E-07, 4.458E-07, 4.657E-07, 4.912E-07, 4.651E-07, 4.907E-07, 6.938E-07, 8.850E-07, 8.135E-07, 4.573E-07, 3.349E-07, 2.437E-07, 2.148E-07, 1.662E-07, 8.037E-08, 7.841E-08, 6.411E-08, 2.551E-08, 2.056E-08, 3.058E-08, 1.083E-07, 3.938E-07, 7.839E-07, 5.837E-07, 4.154E-07, 4.458E-07, 4.657E-07, 4.912E-07 /)
    ALD2_emi = (/ 1.702E-09, 1.283E-09, 9.397E-10, 1.024E-09, 1.076E-09, 1.132E-09, 1.068E-09, 1.130E-09, 1.651E-09, 2.132E-09, 1.985E-09, 1.081E-09, 7.847E-10, 5.676E-10, 5.003E-10, 3.838E-10, 1.784E-10, 1.766E-10, 1.430E-10, 5.173E-11, 4.028E-11, 6.349E-11, 2.428E-10, 8.716E-10, 1.702E-09, 1.283E-09, 9.397E-10, 1.024E-09, 1.076E-09, 1.132E-09 /)
    FORM_emi = (/ 4.061E-09, 3.225E-09, 2.440E-09, 2.639E-09, 2.754E-09, 2.888E-09, 2.741E-09, 2.885E-09, 4.088E-09, 5.186E-09, 4.702E-09, 2.601E-09, 1.923E-09, 1.412E-09, 1.252E-09, 9.776E-10, 4.687E-10, 4.657E-10, 3.836E-10, 1.717E-10, 1.448E-10, 1.976E-10, 6.193E-10, 2.090E-09, 4.061E-09, 3.225E-09, 2.440E-09, 2.639E-09, 2.754E-09, 2.888E-09 /)
    ETH_emi = (/ 1.849E-08, 1.391E-08, 1.010E-08, 1.095E-08, 1.148E-08, 1.209E-08, 1.142E-08, 1.205E-08, 1.806E-08, 2.320E-08, 2.149E-08, 1.146E-08, 8.384E-09, 6.124E-09, 5.414E-09, 4.119E-09, 1.953E-09, 1.927E-09, 1.575E-09, 6.164E-10, 4.973E-10, 7.420E-10, 2.653E-09, 9.477E-09, 1.849E-08, 1.391E-08, 1.010E-08, 1.095E-08, 1.148E-08, 1.209E-08 /)
    IOLE_emi = (/ 5.948E-09, 4.573E-09, 3.374E-09, 3.668E-09, 3.851E-09, 4.050E-09, 3.841E-09, 4.052E-09, 6.094E-09, 7.795E-09, 7.215E-09, 3.738E-09, 2.718E-09, 1.973E-09, 1.729E-09, 1.338E-09, 6.333E-10, 6.394E-10, 5.126E-10, 2.089E-10, 1.708E-10, 2.480E-10, 8.947E-10, 3.057E-09, 5.948E-09, 4.573E-09, 3.374E-09, 3.668E-09, 3.851E-09, 4.050E-09 /)
    OLE_emi = (/ 5.948E-09, 4.573E-09, 3.374E-09, 3.668E-09, 3.851E-09, 4.050E-09, 3.841E-09, 4.052E-09, 6.094E-09, 7.795E-09, 7.215E-09, 3.738E-09, 2.718E-09, 1.973E-09, 1.729E-09, 1.338E-09, 6.333E-10, 6.394E-10, 5.126E-10, 2.089E-10, 1.708E-10, 2.480E-10, 8.947E-10, 3.057E-09, 5.948E-09, 4.573E-09, 3.374E-09, 3.668E-09, 3.851E-09, 4.050E-09 /)
    TOL_emi = (/ 6.101E-09, 8.706E-09, 7.755E-09, 8.024E-09, 8.202E-09, 8.410E-09, 8.218E-09, 8.407E-09, 1.020E-08, 1.139E-08, 7.338E-09, 4.184E-09, 3.078E-09, 2.283E-09, 2.010E-09, 1.575E-09, 8.966E-10, 6.705E-10, 5.395E-10, 2.462E-10, 2.106E-10, 2.852E-10, 9.300E-10, 3.144E-09, 6.101E-09, 8.706E-09, 7.755E-09, 8.024E-09, 8.202E-09, 8.410E-09 /)
    XYL_emi = (/ 5.599E-09, 4.774E-09, 3.660E-09, 3.909E-09, 4.060E-09, 4.239E-09, 4.060E-09, 4.257E-09, 6.036E-09, 7.448E-09, 6.452E-09, 3.435E-09, 2.525E-09, 1.859E-09, 1.650E-09, 1.302E-09, 6.852E-10, 6.773E-10, 5.437E-10, 2.697E-10, 2.358E-10, 3.059E-10, 8.552E-10, 2.861E-10, 5.599E-09, 4.774E-09, 3.660E-09, 3.909E-09, 4.060E-09, 4.239E-09 /)
    PAR_emi = (/ 1.709E-07, 1.953E-07, 1.698E-07, 1.761E-07, 1.808E-07, 1.865E-07, 1.822E-07, 1.859E-07, 2.412E-07, 2.728E-07, 2.174E-07, 1.243E-07, 9.741E-08, 7.744E-08, 6.931E-08, 5.805E-08, 3.900E-08, 3.317E-08, 2.956E-08, 2.306E-08, 2.231E-08, 2.395E-08, 4.284E-08, 9.655E-08, 1.709E-07, 1.953E-07, 1.698E-07, 1.761E-07, 1.808E-07, 1.865E-07 /)
    ISOP_emi = (/ 2.412E-10, 2.814E-10, 3.147E-10, 4.358E-10, 5.907E-10, 6.766E-10, 6.594E-10, 5.879E-10, 5.435E-10, 6.402E-10, 5.097E-10, 9.990E-11, 7.691E-11, 5.939E-11, 5.198E-11, 4.498E-11, 3.358E-11, 2.946E-11, 2.728E-11, 2.183E-11, 1.953E-11, 1.890E-11, 2.948E-11, 1.635E-10, 2.412E-10, 2.814E-10, 3.147E-10, 4.358E-10, 5.907E-10, 6.766E-10 /)
    MEOH_emi = (/ 2.368E-10, 6.107E-10, 6.890E-10, 6.890E-10, 6.890E-10, 6.889E-10, 6.886E-10, 6.890E-10, 6.890E-10, 5.414E-10, 3.701E-10, 2.554E-10, 1.423E-10, 6.699E-11, 2.912E-11, 2.877E-11, 2.825E-11, 2.056E-12, 2.056E-12, 2.056E-12, 2.435E-12, 2.435E-12, 4.030E-11, 1.168E-10, 2.368E-10, 6.107E-10, 6.890E-10, 6.890E-10, 6.890E-10, 6.889E-10 /)
    rate_emi = (/0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

#endif

    call cpu_time(comp_start)

    if(.not.this%solve_multiple_cells) then
      do i=i_start, i_end
        do j=j_start, j_end
          do k=1, k_end

            ! Calculate the vertical index for NMMB-style arrays
            k_flip = size(MONARCH_conc,3) - k + 1

            ! Update the environmental state
            call this%camp_state%env_states(1)%set_temperature_K( &
              real( temperature(i,j,k_flip), kind=dp ) )
            call this%camp_state%env_states(1)%set_pressure_Pa(   &
              real( pressure(i,k,j), kind=dp ) )

            this%camp_state%state_var(:) = 0.0

            this%camp_state%state_var(this%map_camp_id(:)) = &
                    this%camp_state%state_var(this%map_camp_id(:)) + &
                            MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))
            this%camp_state%state_var(this%gas_phase_water_id) = &
                    water_conc(i,j,k_flip,water_vapor_index) * &
                    air_density(i,k,j) * 1.0d9

            if(mod(int(start_time),60).eq.0) then
              write(*,*) "i_hour loop", i_hour
            end if
#ifndef ENABLE_CB05_SOA
              !Add emissions
              this%camp_state%state_var(chem_spec_data%gas_state_id("SO2"))=this%camp_state%state_var(chem_spec_data%gas_state_id("SO2"))+SO2_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("NO2"))=this%camp_state%state_var(chem_spec_data%gas_state_id("NO2"))+NO2_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("NO"))=this%camp_state%state_var(chem_spec_data%gas_state_id("NO"))+NO_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("NH3"))=this%camp_state%state_var(chem_spec_data%gas_state_id("NH3"))+NH3_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("CO"))=this%camp_state%state_var(chem_spec_data%gas_state_id("CO"))+CO_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("ALD2"))=this%camp_state%state_var(chem_spec_data%gas_state_id("ALD2"))+ALD2_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("FORM"))=this%camp_state%state_var(chem_spec_data%gas_state_id("FORM"))+FORM_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("ETH"))=this%camp_state%state_var(chem_spec_data%gas_state_id("ETH"))+ETH_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("IOLE"))=this%camp_state%state_var(chem_spec_data%gas_state_id("IOLE"))+IOLE_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("OLE"))=this%camp_state%state_var(chem_spec_data%gas_state_id("OLE"))+OLE_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("TOL"))=this%camp_state%state_var(chem_spec_data%gas_state_id("TOL"))+TOL_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("XYL"))=this%camp_state%state_var(chem_spec_data%gas_state_id("XYL"))+XYL_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("PAR"))=this%camp_state%state_var(chem_spec_data%gas_state_id("PAR"))+PAR_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("ISOP"))=this%camp_state%state_var(chem_spec_data%gas_state_id("ISOP"))+ISOP_emi(i_hour)*rate_emi(i_hour)*conv
              this%camp_state%state_var(chem_spec_data%gas_state_id("MEOH"))=this%camp_state%state_var(chem_spec_data%gas_state_id("MEOH"))+MEOH_emi(i_hour)*rate_emi(i_hour)*conv
#endif


            !write(*,*) "State_var input",this%camp_state%state_var(this%map_camp_id(:)+(z*state_size_per_cell))

            ! Start the computation timer
            if (MONARCH_PROCESS.eq.0 .and. i.eq.i_start .and. j.eq.j_start &
                    .and. k.eq.1) then
              !solver_stats%debug_out = .false.
            else
              !solver_stats%debug_out = .false.
            end if

            ! Integrate the PMC mechanism
            call this%camp_core%solve(this%camp_state, &
                    real(time_step, kind=dp), solver_stats = solver_stats)

            call assert_msg(376450931, solver_stats%status_code.eq.0, &
                            "Solver failed with code "// &
                            to_string(solver_stats%solver_flag))

#ifdef PMC_DEBUG
            ! Check the Jacobian evaluations
            call assert_msg(611569150, solver_stats%Jac_eval_fails.eq.0,&
                          trim( to_string( solver_stats%Jac_eval_fails ) )// &
                          " Jacobian evaluation failures at time "// &
                          trim( to_string( start_time ) ) )

            ! Only evaluate the Jacobian for the first cell because it is
            ! time consuming
            solver_stats%eval_Jac = .false.
#endif

            ! Update the MONARCH tracer array with new species concentrations
            MONARCH_conc(i,j,k_flip,this%map_monarch_id(:)) = &
                    this%camp_state%state_var(this%map_camp_id(:))

          end do
        end do
      end do

    else

      ! solve multiple grid cells at once
      !  FIXME this only works if this%n_cells ==
      !       (i_end - i_start + 1) * (j_end - j_start + 1 ) * k_end
      !n_cell_check = (i_end - i_start + 1) * (j_end - j_start + 1 ) * k_end
      !call assert_msg(559245176, this%n_cells .eq. n_cell_check, &
      !        "Grid cell number mismatch, got "// &
      !                trim(to_string(n_cell_check))//", expected "// &
      !                trim(to_string(this%n_cells)))

      ! Set initial conditions and environmental parameters for each grid cell
      do i=i_start, i_end
        do j=j_start, j_end
          do k=1, k_end
            !Remember fortran read matrix in inverse order for optimization!
            ! TODO add descriptions for o and z, or preferably use descriptive
            !      variable names
            o = (j-1)*(i_end) + (i-1) !Index to 3D
            z = (k-1)*(i_end*j_end) + o !Index for 2D

            ! Calculate the vertical index for NMMB-style arrays
            k_flip = size(MONARCH_conc,3) - k + 1

            ! Update the environmental state
            !call this%camp_state%env_states(1)%set_temperature_K(real(temperature(i,j,k_flip),kind=dp))
            !call this%camp_state%env_states(1)%set_pressure_Pa(real(pressure(i,k,j),kind=dp))
            call this%camp_state%env_states(z+1)%set_temperature_K(real(temperature(i,j,k_flip),kind=dp))
            call this%camp_state%env_states(z+1)%set_pressure_Pa(real(pressure(i,k,j),kind=dp))

            !write(*,*) "State_var input",this%camp_state%state_var(this%map_camp_id(:)+(z*state_size_per_cell))
            !write(*,*) "Monarch_conc input", MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))
            !todo fix better Nan values
            this%camp_state%state_var(this%map_camp_id(:)+(z*state_size_per_cell))=&
                    this%camp_state%state_var(this%map_camp_id(:))
            MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))=MONARCH_conc(1,1,1,this%map_monarch_id(:))
            !write(*,*) "Monarch_conc input", MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))
            !write(*,*) "State_var input",this%camp_state%state_var(this%map_camp_id(:)+(z*state_size_per_cell))

            !Reset state conc
            this%camp_state%state_var(this%map_camp_id(:) + &
                                       (z*state_size_per_cell)) = 0.0

            this%camp_state%state_var(this%map_camp_id(:) + &
                                       (z*state_size_per_cell)) = &
                    this%camp_state%state_var(this%map_camp_id(:) + &
                                               (z*state_size_per_cell)) + &
                    MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))
            this%camp_state%state_var(this%gas_phase_water_id + &
                                       (z*state_size_per_cell)) = &
                    water_conc(i,j,k_flip,water_vapor_index) * &
                          air_density(i,k,j) * 1.0d9

          end do
        end do
      end do

      ! Integrate the PMC mechanism
      call this%camp_core%solve(this%camp_state, &
              real(time_step, kind=dp), solver_stats = solver_stats)

      do i=i_start, i_end
        do j=j_start, j_end
          do k=1, k_end
            o = (j-1)*(i_end) + (i-1) !Index to 3D
            z = (k-1)*(i_end*j_end) + o !Index for 2D

            k_flip = size(MONARCH_conc,3) - k + 1
            MONARCH_conc(i,j,k_flip,this%map_monarch_id(:)) = &
                    this%camp_state%state_var(this%map_camp_id(:) + &
                                               (z*state_size_per_cell))
          end do
        end do
      end do

    end if

#ifndef ENABLE_CB05_SOA
    if(mod(int(start_time),60).eq.0) then
      i_hour=i_hour+1
    end if
#endif

    !W8 until all process to send data and measure correctly times
#ifdef PMC_USE_MPI
    call pmc_mpi_barrier()
#endif

    !todo move this comp_time to only the solve part, like the other tests
    call cpu_time(comp_end)
    comp_time = comp_time + (comp_end-comp_start)


#ifdef PMC_USE_MPI

if (pmc_mpi_rank().eq.0) then
  !call solver_stats%print( )
end if

#endif

  end subroutine integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load the MONARCH <-> PartMC-camp interface input data
  subroutine load(this, config_file)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this
    !> Interface configuration file path
    character(len=:), allocatable :: config_file

#ifdef PMC_USE_JSON

    type(json_core), pointer :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next, j_child
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val

    character(len=:), allocatable :: str_val
    integer(kind=i_kind) :: var_type
    logical :: found

    ! Initialize the property sets
    this%species_map_data => property_t()
    this%init_conc_data => property_t()
    this%property_set => property_t()

    ! Get a new json core
    allocate(json)

    ! Initialize the json objects
    j_obj => null()
    j_next => null()

    ! Initialize the json file
    call j_file%initialize()
    call j_file%get_core(json)
    call assert_msg(207035903, allocated(config_file), &
              "Received non-allocated string for file path")
    call assert_msg(368569727, trim(config_file).ne."", &
              "Received empty string for file path")
    inquire( file=config_file, exist=found )
    call assert_msg(134309013, found, "Cannot find file: "// &
              config_file)
    call j_file%load_file(filename = config_file)

    ! Find the interface data
    call j_file%get('monarch-data(1)', j_obj)

    ! Load the data to the property_set
    do while (associated(j_obj))

      ! Find the object type
      call json%get(j_obj, 'type', unicode_str_val, found)
      call assert_msg(236838162, found, "Missing type in json input file "// &
              config_file)
      str_val = unicode_str_val

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Load property sets according to type !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Species Map data
      if (str_val.eq."SPECIES_MAP") then
        call json%get_child(j_obj, j_child)
        do while (associated(j_child))
          call json%info(j_child, name=key, var_type=var_type)
          if (key.ne."type".and.key.ne."name") then
            call this%species_map_data%load(json, j_child, .false., key)
          end if
          j_next => j_child
          call json%get_next(j_next, j_child)
        end do

      ! Initial concentration data
      else if (str_val.eq."INIT_CONC") then
        call json%get_child(j_obj, j_child)
        do while (associated(j_child))
          call json%info(j_child, name=key, var_type=var_type)
          if (key.ne."type".and.key.ne."name") then
            call this%init_conc_data%load(json, j_child, .false., key)
          end if
          j_next => j_child
          call json%get_next(j_next, j_child)
        end do

      ! Data of unknown type
      else
        call this%property_set%load(json, j_obj, .false., str_val)
      end if

      j_next => j_obj
      call json%get_next(j_next, j_obj)
    end do

    ! Clean up the json objects
    call j_file%destroy()
    call json%destroy()
    deallocate(json)

#else
    call die_msg(635417227, "PartMC-camp <-> MONARCH interface requires "// &
                  "JSON file support.")
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create the PartMC-camp <-> MONARCH species map
  subroutine create_map(this)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(property_t), pointer :: gas_species_list, aero_species_list, species_data
    character(len=:), allocatable :: key_name, spec_name, rep_name
    integer(kind=i_kind) :: i_spec, num_spec

    integer :: i_rxn, i_photo_rxn, i_base_rate
    type(mechanism_data_t), pointer :: mechanism
    class(rxn_data_t), pointer :: rxn
    character(len=:), allocatable :: key, str_val
    real(kind=dp) :: base_rate
    !real(kind=dp), allocatable :: photo_rates(:)

    key = "MONARCH mod37"

    !mechanism => this%camp_core%mechanism( 1 ) %val
    call assert(418262750, this%camp_core%get_mechanism(key, mechanism))

    !key = "Fast-J id"
    key = "base rate"

    this%n_photo_rxn = 0
    do i_rxn = 1, mechanism%size( )
      rxn => mechanism%get_rxn( i_rxn )
      select type( rxn_photo => rxn )
      class is( rxn_photolysis_t )
        !if( rxn%property_set%get_string( key, str_val ) ) then
        if( rxn%property_set%get_real( key, base_rate ) ) then
          this%n_photo_rxn = this%n_photo_rxn + 1
        end if
      end select
    end do

    allocate(this%rate_update(this%n_photo_rxn))
    allocate(this%photo_rates(this%n_photo_rxn))

    ! Set the PMC and Fast-J ids
    i_photo_rxn = 1
    do i_rxn = 1, mechanism%size( )
      rxn => mechanism%get_rxn( i_rxn )
      select type( rxn_photo => rxn )
      class is( rxn_photolysis_t ) !type is (rxn_photolysis_t)
        !if( rxn%property_set%get_string( key, str_val ) ) then
        if( rxn%property_set%get_real( key, base_rate ) ) then

          ! Initialize the update data object for this reaction
          call this%camp_core%initialize_update_object( rxn,               &
                  this%rate_update(i_photo_rxn))

          !Update photo_rate
          this%photo_rates(i_photo_rxn)=base_rate
          !call this%rate_update(i_photo_rxn)%set_rate(base_rate)
          !call this%camp_core%update_data(this%rate_update(i_photo_rxn))

          !write(*,*) "base rate", base_rate

          i_photo_rxn = i_photo_rxn+1
        end if
      end select
    end do

    !do i_photo_rxn = 1, this%n_photo_rxn
    !  call this%rate_update(i_photo_rxn)%set_rate(real(this%photo_rates(i_photo_rxn), kind=dp))
    !  call this%camp_core%update_data(this%rate_update(i_photo_rxn))
    !end do

    ! Make sure nothing strange happened
    call assert( 410333867, this%n_photo_rxn .eq. i_photo_rxn - 1 )

    ! Get the gas-phase species ids
    key_name = "gas-phase species"
    call assert_msg(939097252, &
            this%species_map_data%get_property_t(key_name, gas_species_list), &
            "Missing set of gas-phase species MONARCH ids")
    num_spec = gas_species_list%size()

    ! Get the aerosol-phase species ids
    key_name = "aerosol-phase species"
    if (this%species_map_data%get_property_t(key_name, &
            aero_species_list)) then
      num_spec = num_spec + aero_species_list%size()
    end if

    ! Set up the species map and MONARCH names array
    allocate(this%monarch_species_names(num_spec))
    allocate(this%map_monarch_id(num_spec))
    allocate(this%map_camp_id(num_spec))

    ! Get the chemical species data
    call assert_msg(731700229, &
            this%camp_core%get_chem_spec_data(chem_spec_data), &
            "No chemical species data in camp_core.")

    ! Set the gas-phase water id
    key_name = "gas-phase water"
    call assert_msg(413656652, &
            this%species_map_data%get_string(key_name, spec_name), &
            "Missing gas-phase water species for MONARCH interface.")
    this%gas_phase_water_id = chem_spec_data%gas_state_id(spec_name)
    call assert_msg(910692272, this%gas_phase_water_id.gt.0, &
            "Could not find gas-phase water species '"//spec_name//"'.")

    ! Loop through the gas-phase species and set up the map
    call gas_species_list%iter_reset()
    i_spec = 1
    do while (gas_species_list%get_key(spec_name))

      this%monarch_species_names(i_spec)%string = spec_name

      call assert_msg(599522862, &
              gas_species_list%get_property_t(val=species_data), &
              "Missing species data for '"//spec_name//"' in PartMC-camp "// &
              "<-> MONARCH species map.")

      key_name = "monarch id"
      call assert_msg(643926329, &
              species_data%get_int(key_name, this%map_monarch_id(i_spec)), &
              "Missing monarch id for species '"//spec_name//" in "// &
              "PartMC-camp <-> MONARCH species map.")
      this%map_monarch_id(i_spec) = this%map_monarch_id(i_spec) + &
              this%tracer_starting_id - 1
      call assert_msg(450258014, &
              this%map_monarch_id(i_spec).le.this%tracer_ending_id, &
              "Monarch id for species '"//spec_name//"' out of specified "// &
              "tracer array bounds.")

      this%map_camp_id(i_spec) = chem_spec_data%gas_state_id(spec_name)
      call assert_msg(916977002, this%map_camp_id(i_spec).gt.0, &
                "Could not find species '"//spec_name//"' in PartMC-camp.")

      call gas_species_list%iter_next()
      i_spec = i_spec + 1
    end do

    ! Loop through the aerosol-phase species and add them to the map
    if (associated(aero_species_list)) then

      call aero_species_list%iter_reset()
      do while(aero_species_list%get_key(spec_name))

        this%monarch_species_names(i_spec)%string = spec_name

        call assert_msg(567689501, &
                aero_species_list%get_property_t(val=species_data), &
                "Missing species data for '"//spec_name//"' in " //&
                "PartMC-camp <-> MONARCH species map.")

        key_name = "monarch id"
        call assert_msg(615451741, &
                species_data%get_int(key_name, this%map_monarch_id(i_spec)), &
                "Missing monarch id for species '"//spec_name//"' in "// &
                "PartMC-camp <-> MONARCH species map.")
        this%map_monarch_id(i_spec) = this%map_monarch_id(i_spec) + &
                this%tracer_starting_id - 1
        call assert_msg(382644266, &
                this%map_monarch_id(i_spec).le.this%tracer_ending_id, &
                "Monarch id for species '"//spec_name//"' out of "// &
                "specified tracer array bounds.")

        key_name = "aerosol representation name"
        call assert_msg(963222513, &
                species_data%get_string(key_name, rep_name), &
                "Missing aerosol representation name for species '"// &
                spec_name//"' in PartMC-camp <-> MONARCH species map.")

        ! Find the species PartMC id
        this%map_camp_id(i_spec) = 0
        call assert_msg(377850668, &
                this%camp_core%get_aero_rep(rep_name, aero_rep_ptr), &
                "Could not find aerosol representation '"//rep_name//"'")
        this%map_camp_id(i_spec) = aero_rep_ptr%spec_state_id(spec_name)
        call assert_msg(887136850, this%map_camp_id(i_spec) .gt. 0, &
                "Could not find aerosol species '"//spec_name//"' in "// &
                "aerosol representation '"//rep_name//"'.")

        !Conver from [ug m-3] to [kg m-3] (not needed)
        !this%map_monarch_id(i_spec)=this%map_monarch_id(i_spec)*1.0e-9

        call aero_species_list%iter_next()
        i_spec = i_spec + 1
      end do
    end if

  end subroutine create_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load initial concentrations
  subroutine load_init_conc(this)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(property_t), pointer :: gas_species_list, aero_species_list, species_data
    character(len=:), allocatable :: key_name, spec_name, rep_name
    integer(kind=i_kind) :: i_spec, num_spec

    num_spec = 0

    ! Get the gas-phase species
    key_name = "gas-phase species"
    if (this%init_conc_data%get_property_t(key_name, gas_species_list)) then
      num_spec = num_spec + gas_species_list%size()
    end if

    ! Get the aerosol-phase species
    key_name = "aerosol-phase species"
    if (this%init_conc_data%get_property_t(key_name, aero_species_list)) then
      num_spec = num_spec + aero_species_list%size()
    end if

    ! Get the chemical species data
    call assert_msg(885063268, &
            this%camp_core%get_chem_spec_data(chem_spec_data), &
            "No chemical species data in camp_core.")

    ! Allocate space for the initial concentrations and indices
    allocate(this%init_conc_camp_id(num_spec))
    allocate(this%init_conc(num_spec))

    ! Add the gas-phase initial concentrations
    if (associated(gas_species_list)) then

      ! Loop through the gas-phase species and load the initial concentrations
      call gas_species_list%iter_reset()
      i_spec = 1
      do while (gas_species_list%get_key(spec_name))

        call assert_msg(325582312, &
                gas_species_list%get_property_t(val=species_data), &
                "Missing species data for '"//spec_name//"' for "// &
                "PartMC-camp initial concentrations.")

        key_name = "init conc"
        call assert_msg(445070498, &
                species_data%get_real(key_name, this%init_conc(i_spec)), &
                "Missing 'init conc' for species '"//spec_name//" for "// &
                "PartMC-camp initial concentrations.")

        this%init_conc_camp_id(i_spec) = &
                chem_spec_data%gas_state_id(spec_name)
        call assert_msg(940200584, this%init_conc_camp_id(i_spec).gt.0, &
                "Could not find species '"//spec_name//"' in PartMC-camp.")

        call gas_species_list%iter_next()
        i_spec = i_spec + 1
      end do

    end if

    ! Add the aerosol-phase species initial concentrations
    if (associated(aero_species_list)) then

      call aero_species_list%iter_reset()
      do while(aero_species_list%get_key(spec_name))

        call assert_msg(331096555, &
                aero_species_list%get_property_t(val=species_data), &
                "Missing species data for '"//spec_name//"' for " //&
                "PartMC-camp initial concentrations.")

        key_name = "init conc"
        call assert_msg(782275469, &
                species_data%get_real(key_name, this%init_conc(i_spec)), &
                "Missing 'init conc' for species '"//spec_name//"' for "// &
                "PartMC-camp initial concentrations.")

        key_name = "aerosol representation name"
        call assert_msg(150863332, &
                species_data%get_string(key_name, rep_name), &
                "Missing aerosol representation name for species '"// &
                spec_name//"' for PartMC-camp initial concentrations.")

        ! Find the species PartMC id
        this%init_conc_camp_id(i_spec) = 0
        call assert_msg(258814777, &
                this%camp_core%get_aero_rep(rep_name, aero_rep_ptr), &
                "Could not find aerosol representation '"//rep_name//"'")
        this%init_conc_camp_id(i_spec) = &
                aero_rep_ptr%spec_state_id(spec_name)
        call assert_msg(437149649, this%init_conc_camp_id(i_spec) .gt. 0, &
                "Could not find aerosol species '"//spec_name//"' in "// &
                "aerosol representation '"//rep_name//"'.")

        call aero_species_list%iter_next()
        i_spec = i_spec + 1
      end do
    end if

    !write(*,*) "Init conc",this%init_conc(:)

  end subroutine load_init_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get initial concentrations for the mock MONARCH model (for testing only)
  subroutine get_init_conc(this, MONARCH_conc, MONARCH_water_conc, &
      WATER_VAPOR_ID, MONARCH_air_density)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this
    !> MONARCH species concentrations to update
    real, intent(inout) :: MONARCH_conc(:,:,:,:)
    !> Atmospheric water concentrations (kg_H2O/kg_air)
    real, intent(out) :: MONARCH_water_conc(:,:,:,:)
    !> Index in water_conc corresponding to water vapor
    integer, intent(in) :: WATER_VAPOR_ID
    !> Air density (kg_air/m^3)
    real, intent(out) :: MONARCH_air_density(:,:,:)

    integer(kind=i_kind) :: i_spec, water_id, i,j,k
    real :: factor_ppb_to_ppm

#ifndef ENABLE_CB05_SOA
    factor_ppb_to_ppm=1.0E-3
#else
    factor_ppb_to_ppm=1.0
#endif

    ! Reset the species concentrations in PMC and MONARCH
    this%camp_state%state_var(:) = 0.0
    MONARCH_conc(:,:,:,:) = 0.0
    MONARCH_water_conc(:,:,:,WATER_VAPOR_ID) = 0.0

    ! Set the air density to a nominal value
    MONARCH_air_density(:,:,:) = 1.225

    ! Set initial concentrations in PMC
    this%init_conc(:) = this%init_conc(:) * factor_ppb_to_ppm
    this%camp_state%state_var(this%init_conc_camp_id(:)) = this%init_conc(:)

    ! Copy species concentrations to MONARCH array
    forall (i_spec = 1:size(this%map_monarch_id))
      MONARCH_conc(:,:,:,this%map_monarch_id(i_spec)) = &
              this%camp_state%state_var(this%map_camp_id(i_spec))
    end forall

    !write(*,*) "Init conc",this%init_conc(this%map_camp_id(:))
    !write(*,*) "State_var init",this%camp_state%state_var(this%map_camp_id(:))

    ! Set the relative humidity
    MONARCH_water_conc(:,:,:,WATER_VAPOR_ID) = &
            this%camp_state%state_var(this%gas_phase_water_id) * &
            1.0d-9 / 1.225d0

  end subroutine get_init_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the MONARCH species names and indices (for testing only)
  subroutine get_MONARCH_species(this, species_names, MONARCH_ids)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this
    !> Set of MONARCH species names
    type(string_t), allocatable, intent(out) :: species_names(:)
    !> MONARCH tracer ids
    integer(kind=i_kind), allocatable, intent(out) :: MONARCH_ids(:)

    species_names = this%monarch_species_names
    MONARCH_ids = this%map_monarch_id

  end subroutine get_MONARCH_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the PartMC-camp data
  subroutine do_print(this)

    !> PartMC-camp <-> MONARCH interface
    class(monarch_interface_t) :: this

    call this%camp_core%print()

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the interface
  elemental subroutine finalize(this)

    !> PartMC-camp <-> MONARCH interface
    type(monarch_interface_t), intent(inout) :: this

    if (associated(this%camp_core)) &
            deallocate(this%camp_core)
    if (associated(this%camp_state)) &
            deallocate(this%camp_state)
    if (allocated(this%monarch_species_names)) &
            deallocate(this%monarch_species_names)
    if (allocated(this%map_monarch_id)) &
            deallocate(this%map_monarch_id)
    if (allocated(this%map_camp_id)) &
            deallocate(this%map_camp_id)
    if (allocated(this%init_conc_camp_id)) &
            deallocate(this%init_conc_camp_id)
    if (allocated(this%init_conc)) &
            deallocate(this%init_conc)
    if (associated(this%species_map_data)) &
            deallocate(this%species_map_data)
    if (associated(this%init_conc_data)) &
            deallocate(this%init_conc_data)
    if (associated(this%property_set)) &
            deallocate(this%property_set)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_monarch_interface
