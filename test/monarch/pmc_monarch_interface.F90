! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The monarch_interface_t object and related functions

!> Interface for the MONACH model and PartMC-camp
module pmc_monarch_interface

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : assert_msg, string_t
  use pmc_camp_core
  use pmc_camp_state
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_modal_binned_mass
  use pmc_chem_spec_data
  use pmc_property
  use pmc_camp_solver_data
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
    private
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
                       starting_id, ending_id, n_cells, mpi_comm) result (new_obj)

    !> A new MONARCH interface
    type(monarch_interface_t), pointer :: new_obj
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
    integer(kind=i_kind) :: i_spec
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
    allocate(new_obj)

    if (.not.present(n_cells).or.n_cells.eq.1) then
      new_obj%solve_multiple_cells = .false.
    else
      new_obj%solve_multiple_cells = .true.
      new_obj%n_cells=n_cells
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
      call new_obj%load(interface_config_file)

      ! Initialize the camp-chem core
      new_obj%camp_core => camp_core_t(camp_config_file, new_obj%n_cells)
      call new_obj%camp_core%initialize()

      ! Set the aerosol representation id
      if (new_obj%camp_core%get_aero_rep("MONARCH mass-based", aero_rep)) then
        select type (aero_rep)
          type is (aero_rep_modal_binned_mass_t)
            call new_obj%camp_core%initialize_update_object( aero_rep, &
                                                             update_data_GMD)
            call new_obj%camp_core%initialize_update_object( aero_rep, &
                                                             update_data_GSD)
            call assert(889473105, &
                        aero_rep%get_section_id("organic matter", i_sect_om))
            call assert(648042550, &
                        aero_rep%get_section_id("black carbon", i_sect_bc))
            call assert(760360895, &
                        aero_rep%get_section_id("sulfate", i_sect_sulf))
            call assert(307728742, &
                        aero_rep%get_section_id("other PM", i_sect_opm))
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
      new_obj%tracer_starting_id = starting_id
      new_obj%tracer_ending_id = ending_id

      ! Generate the PartMC-camp <-> MONARCH species map
      call new_obj%create_map()

      ! Load the initial concentrations
      call new_obj%load_init_conc()

#ifdef PMC_USE_MPI
      pack_size = new_obj%camp_core%pack_size() + &
              update_data_GMD%pack_size() + &
              update_data_GSD%pack_size() + &
              pmc_mpi_pack_size_integer_array(new_obj%map_monarch_id) + &
              pmc_mpi_pack_size_integer_array(new_obj%map_camp_id) + &
              pmc_mpi_pack_size_integer_array(new_obj%init_conc_camp_id) + &
              pmc_mpi_pack_size_real_array(new_obj%init_conc) + &
              pmc_mpi_pack_size_integer(new_obj%gas_phase_water_id) + &
              pmc_mpi_pack_size_integer(i_sect_om) + &
              pmc_mpi_pack_size_integer(i_sect_bc) + &
              pmc_mpi_pack_size_integer(i_sect_sulf) + &
              pmc_mpi_pack_size_integer(i_sect_opm)
      allocate(buffer(pack_size))
      pos = 0
      call new_obj%camp_core%bin_pack(buffer, pos)
      call update_data_GMD%bin_pack(buffer, pos)
      call update_data_GSD%bin_pack(buffer, pos)
      call pmc_mpi_pack_integer_array(buffer, pos, new_obj%map_monarch_id)
      call pmc_mpi_pack_integer_array(buffer, pos, new_obj%map_camp_id)
      call pmc_mpi_pack_integer_array(buffer, pos, new_obj%init_conc_camp_id)
      call pmc_mpi_pack_real_array(buffer, pos, new_obj%init_conc)
      call pmc_mpi_pack_integer(buffer, pos, new_obj%gas_phase_water_id)
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
      new_obj%camp_core => camp_core_t()
      pos = 0
      call new_obj%camp_core%bin_unpack(buffer, pos)
      call update_data_GMD%bin_unpack(buffer, pos)
      call update_data_GSD%bin_unpack(buffer, pos)
      call pmc_mpi_unpack_integer_array(buffer, pos, new_obj%map_monarch_id)
      call pmc_mpi_unpack_integer_array(buffer, pos, new_obj%map_camp_id)
      call pmc_mpi_unpack_integer_array(buffer, pos, new_obj%init_conc_camp_id)
      call pmc_mpi_unpack_real_array(buffer, pos, new_obj%init_conc)
      call pmc_mpi_unpack_integer(buffer, pos, new_obj%gas_phase_water_id)
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
    call new_obj%camp_core%solver_initialize()

    ! Create a state variable on each node
    new_obj%camp_state => new_obj%camp_core%new_state()

    ! Set the aerosol mode dimensions

    ! organic matter
    if (i_sect_om.gt.0) then
      call update_data_GMD%set_GMD(i_sect_om, 2.12d-8)
      call update_data_GSD%set_GSD(i_sect_om, 2.24d0)
      call new_obj%camp_core%update_data(update_data_GMD)
      call new_obj%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_bc.gt.0) then
    ! black carbon
      call update_data_GMD%set_GMD(i_sect_bc, 1.18d-8)
      call update_data_GSD%set_GSD(i_sect_bc, 2.00d0)
      call new_obj%camp_core%update_data(update_data_GMD)
      call new_obj%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_sulf.gt.0) then
    ! sulfate
      call update_data_GMD%set_GMD(i_sect_sulf, 6.95d-8)
      call update_data_GSD%set_GSD(i_sect_sulf, 2.12d0)
      call new_obj%camp_core%update_data(update_data_GMD)
      call new_obj%camp_core%update_data(update_data_GSD)
    end if
    if (i_sect_opm.gt.0) then
    ! other PM
      call update_data_GMD%set_GMD(i_sect_opm, 2.12d-8)
      call update_data_GSD%set_GSD(i_sect_opm, 2.24d0)
      call new_obj%camp_core%update_data(update_data_GMD)
      call new_obj%camp_core%update_data(update_data_GSD)
    end if

    ! Calculate the intialization time
    if (MONARCH_PROCESS.eq.0) then
      call cpu_time(comp_end)
      write(*,*) "Initialization time: ", comp_end-comp_start, " s"
      !call new_obj%camp_core%print()
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the PartMC mechanism for a particular set of cells and timestep
  subroutine integrate(this, start_time, time_step, i_start, i_end, j_start, &
                  j_end, temperature, MONARCH_conc, water_conc, &
                  water_vapor_index, air_density, pressure)

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

    integer :: i, j, k, k_flip, i_spec, z, o, i2
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


#if 0
#ifdef PMC_DEBUG
    ! Evaluate the Jacobian during solving
    solver_stats%eval_Jac = .true.
#endif
#endif


    k_end = size(MONARCH_conc,3)

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

#if 0
#ifdef PMC_DEBUG
              ! Check the Jacobian evaluations
              call warn_assert_msg(611569150, solver_stats%Jac_eval_fails.eq.0,&
                          trim( to_string( solver_stats%Jac_eval_fails ) )// &
                          " Jacobian evaluation failures at time "// &
                          trim( to_string( start_time ) ) )
#endif
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
      n_cell_check = (i_end - i_start + 1) * (j_end - j_start + 1 ) * k_end
      call assert_msg(559245176, this%n_cells .eq. n_cell_check, &
              "Grid cell number mismatch, got "// &
                      trim(to_string(n_cell_check))//", expected "// &
                      trim(to_string(this%n_cells)))

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
            call this%camp_state%env_states(1)%set_temperature_K( &
              real( temperature(i,j,k_flip), kind=dp ) )
            call this%camp_state%env_states(1)%set_pressure_Pa(   &
              real( pressure(i,k,j), kind=dp ) )

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

    call cpu_time(comp_end)
    comp_time = comp_time + (comp_end-comp_start)

    ! call solver_stats%print( )

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

    ! Reset the species concentrations in PMC and MONARCH
    this%camp_state%state_var(:) = 0.0
    MONARCH_conc(:,:,:,:) = 0.0
    MONARCH_water_conc(:,:,:,WATER_VAPOR_ID) = 0.0

    ! Set the air density to a nominal value
    MONARCH_air_density(:,:,:) = 1.225

    ! Set initial concentrations in PMC
    this%camp_state%state_var(this%init_conc_camp_id(:)) = &
            this%init_conc(:)

    ! Copy species concentrations to MONARCH array
    forall (i_spec = 1:size(this%map_monarch_id))
      MONARCH_conc(:,:,:,this%map_monarch_id(i_spec)) = &
              this%camp_state%state_var(this%map_camp_id(i_spec))
    end forall

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
