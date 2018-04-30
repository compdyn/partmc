! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The monarch_interface_t object and related functions

!> Interface for the MONACH model and PartMC-phlex
module pmc_monarch_interface

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

  public :: monarch_interface_t

  !> PartMC <-> MONARCH interface
  !!
  !! Contains all data required to intialize and run PartMC from MONARCH data
  !! and map state variables between PartMC and MONARCH
  type :: monarch_interface_t
    private
    !> Phlex-chem core
    type(phlex_core_t), pointer :: phlex_core 
    !> Phlex-chem state
    type(phlex_state_t), pointer :: phlex_state
    !> MONACH <-> PartMC species map
    integer(kind=i_kind), allocatable :: map_monarch_id(:), map_phlex_id(:)
    !> PartMC-phlex ids for initial concentrations
    integer(kind=i_kind), allocatable :: init_conc_phlex_id(:)
    !> Initial species concentrations
    real(kind=dp), allocatable :: init_conc(:)
    !> Starting index for PartMC species on the MONARCH tracer array
    integer(kind=i_kind) :: tracer_starting_id
    !> Ending index for PartMC species on the MONARCH tracer array
    integer(kind=i_kind) :: tracer_ending_id
    !> PartMC-phlex <-> MONARCH species map input data
    type(property_t), pointer :: species_map_data
    !> Initial concentration data
    type(property_t), pointer :: init_conc_data
    !> Interface input data
    type(property_t), pointer :: property_set
  contains
    !> Integrate PartMC for the current MONARCH state over a specified time step
    procedure :: integrate
    !> Get initial concentrations (for testing only)
    procedure :: get_init_conc
    !> Load interface data from a set of input files
    procedure, private :: load
    !> Create the PartMC <-> MONARCH species map
    procedure, private :: create_map
    !> Load the initial concentrations
    procedure, private :: load_init_conc
  end type monarch_interface_t

  !> PartMC <-> MONARCH interface constructor
  interface monarch_interface_t
    procedure :: constructor
  end interface monarch_interface_t

  !> MPI node id from MONARCH
  integer(kind=i_kind) :: MONARCH_NODE ! TODO replace with MONARCH param
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
  function constructor(phlex_config_file, interface_config_file, &
                       starting_id, ending_id) result (new_obj)

    !> A new MONARCH interface
    type(monarch_interface_t), pointer :: new_obj
    !> Path to the PartMC-phlex configuration file list
    character(len=:), allocatable, optional :: phlex_config_file
    !> Path to the PartMC-phlex <-> MONARCH interface input file
    character(len=:), allocatable, optional :: interface_config_file
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

    ! Set the MPI rank (TODO replace with MONARCH param)
    MONARCH_NODE = pmc_mpi_rank()

    ! Create a new interface object
    allocate(new_obj)

    ! Check for an available solver
    call assert_msg(332298164, phlex_solver_data%is_solver_available(), &
            "No solver available")

    ! Initialize the time-invariant model data on each node
    if (MONARCH_NODE.eq.0) then
      
      ! Start the computation timer on the primary node
      call cpu_time(comp_start)

      call assert_msg(304676624, present(phlex_config_file), &
              "Missing PartMC-phlex configuration file list")
      call assert_msg(194027509, present(interface_config_file), &
              "Missing MartMC-phlex <-> MONARCH interface configuration file")
      call assert_msg(937567597, present(starting_id), &
              "Missing starting tracer index for chemical species")
      call assert_msg(593895016, present(ending_id), &
              "Missing ending tracer index for chemical species")

      ! Load the interface data
      call new_obj%load(interface_config_file)

      ! Initialize the phlex-chem core
      new_obj%phlex_core => phlex_core_t(phlex_config_file)
      call new_obj%phlex_core%initialize()

      ! Set the MONARCH tracer array bounds
      new_obj%tracer_starting_id = starting_id
      new_obj%tracer_ending_id = ending_id

      ! Generate the PartMC-phlex <-> MONARCH species map
      call new_obj%create_map()      

      ! Load the initial concentrations
      call new_obj%load_init_conc()
#ifdef PMC_USE_MPI
      pack_size = new_obj%phlex_core%pack_size() + &
              pmc_mpi_pack_size_int_array(new_obj%map_monarch_id) + &
              pmc_mpi_pack_size_int_array(new_obj%map_phlex_id) + &
              pmc_mpi_pack_size_int_array(new_obj%init_conc_phlex_id) + &
              pmc_mpi_pack_size_real_array(new_obj%init_conc)
      allocate(buffer(pack_size))
      pos = 0
      call new_obj%phlex_core%bin_pack(buffer, pos)
      call pmc_mpi_pack_int_array(new_obj%map_monarch_id)
      call pmc_mpi_pack_int_array(new_obj%map_phlex_id)
      call pmc_mpi_pack_int_array(new_obj%init_conc_phlex_id)
      call pmc_mpi_pack_real_array(new_obj%init_conc)
      ! TODO send buffer to all the other nodes
    else
      new_obj%phlex_core => phlex_core_t()
      ! TODO get buffer from the primary node
      pos = 0
      call new_obj%phlex_core%bin_unpack(buffer, pos)
      call pmc_mpi_unpack_int_array(new_obj%map_monarch_id)
      call pmc_mpi_unpack_int_array(new_obj%map_phlex_id)
      call pmc_mpi_unpack_int_array(new_obj%init_conc_phlex_id)
      call pmc_mpi_unpack_real_array(new_obj%init_conc)
#endif
    end if

    ! Initialize the solver on all nodes
    call new_obj%phlex_core%solver_initialize()

    ! Create a state variable on each node
    allocate(new_obj%phlex_state)
    new_obj%phlex_state = new_obj%phlex_core%new_state()

    ! Calculate the intialization time
    if (MONARCH_NODE.eq.0) then
      call cpu_time(comp_end)
      write(*,*) "Initialization time: ", comp_end-comp_start, " s"
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the PartMC mechanism for a particular set of cells and timestep
  subroutine integrate(this, start_time, time_step, i_start, i_end, j_start, &
                  j_end, temperature, MONARCH_conc, water_conc, &
                  water_vapor_index, air_density, pressure)

    !> PartMC-phlex <-> MONARCH interface
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
          this%phlex_state%env_state%temp = temperature(i,j,k_flip)
          this%phlex_state%env_state%pressure = pressure(i,k,j)
          ! TODO finish environmental state setup

          ! Update species concentrations in PMC
          this%phlex_state%state_var(:) = 0.0
          this%phlex_state%state_var(this%map_phlex_id(:)) = &
                  MONARCH_conc(i,j,k_flip,this%map_monarch_id(:))

          ! Start the computation timer
          if (MONARCH_NODE.eq.0 .and. i.eq.i_start .and. j.eq.j_start &
                  .and. k.eq.1) then
            call cpu_time(comp_start)
          end if

          ! Integrate the PMC mechanism
          call this%phlex_core%solve(this%phlex_state, &
                  real(time_step, kind=dp))

          ! Calculate the computation time
          if (MONARCH_NODE.eq.0 .and. i.eq.i_start .and. j.eq.j_start &
                  .and. k.eq.1) then
            call cpu_time(comp_end)
            comp_time = comp_time + (comp_end-comp_start)
          end if

          ! Update the MONARCH tracer array with new species concentrations
          MONARCH_conc(i,j,k_flip,this%map_monarch_id(:)) = &
                  this%phlex_state%state_var(this%map_phlex_id(:))
 
        end do
      end do
    end do

  end subroutine integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load the MONARCH <-> PartMC-phlex interface input data
  subroutine load(this, config_file)

    !> PartMC-phlex <-> MONARCH interface
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
    allocate(this%species_map_data)
    allocate(this%init_conc_data)
    allocate(this%property_set)
    this%species_map_data = property_t()
    this%init_conc_data = property_t()
    this%property_set = property_t()

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
            call this%species_map_data%load(json, j_child, .false.)
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
            call this%init_conc_data%load(json, j_child, .false.)
          end if
          j_next => j_child
          call json%get_next(j_next, j_child)
        end do

      ! Data of unknown type
      else
        call this%property_set%load(json, j_obj, .false.)
      end if

      j_next => j_obj
      call json%get_next(j_next, j_obj)
    end do

    ! Clean up the json objects
    call j_file%destroy()
    call json%destroy()

#else
    call die_msg(635417227, "PartMC-phlex <-> MONARCH interface requires "// &
                  "JSON file support.")
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create the PartMC-phlex <-> MONARCH species map
  subroutine create_map(this)

    !> PartMC-phlex <-> MONARCH interface
    class(monarch_interface_t) :: this

    type(property_t), pointer :: species_list, species_data
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec

    ! Get the gas-phase species ids
    key_name = "gas-phase species"
    call assert_msg(939097252, &
            this%species_map_data%get_property_t(key_name, species_list), &
            "Missing set of gas-phase species MONARCH ids")
    
    ! Set up the species map
    allocate(this%map_monarch_id(species_list%size()))
    allocate(this%map_phlex_id(species_list%size()))

    ! Loop through the gas-phase species and set up the map      
    call species_list%iter_reset()
    i_spec = 1
    do while (species_list%get_key(spec_name))

      call assert_msg(599522862, &
              species_list%get_property_t(val=species_data), &
              "Missing species data for '"//spec_name//"' in PartMC-phlex "// &
              "<-> MONARCH species map.")

      key_name = "monarch id"
      call assert_msg(643926329, &
              species_data%get_int(key_name, this%map_monarch_id(i_spec)), &
              "Missing monarch id for species '"//spec_name//" in "// &
              "PartMC-phlex <-> MONARCH species map.")
      this%map_monarch_id(i_spec) = this%map_monarch_id(i_spec) + &
              this%tracer_starting_id - 1
      call assert_msg(450258014, &
              this%map_monarch_id(i_spec).le.this%tracer_ending_id, &
              "Monarch id for species '"//spec_name//"' out of specified "// &
              "tracer array bounds.")
        
      this%map_phlex_id(i_spec) = &
              this%phlex_core%chem_spec_data%gas_state_id(spec_name)
      call assert_msg(916977002, this%map_phlex_id(i_spec).gt.0, &
                "Could not find species '"//spec_name//"' in PartMC-phlex.")
      
      call species_list%iter_next()
      i_spec = i_spec + 1
    end do

    ! FIXME add aerosol species map creation

  end subroutine create_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load initial concentrations
  subroutine load_init_conc(this)

    !> PartMC-phlex <-> MONARCH interface
    class(monarch_interface_t) :: this

    type(property_t), pointer :: species_list, species_data
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec

    ! Get the gas-phase species ids
    key_name = "gas-phase species"
    if (this%init_conc_data%get_property_t(key_name, species_list)) then
    
      ! Set up the species map
      allocate(this%init_conc_phlex_id(species_list%size()))
      allocate(this%init_conc(species_list%size()))

      ! Loop through the gas-phase species and load the initial concentrations 
      call species_list%iter_reset()
      i_spec = 1
      do while (species_list%get_key(spec_name))

        call assert_msg(325582312, &
                species_list%get_property_t(val=species_data), &
                "Missing species data for '"//spec_name//"' for "// &
                "PartMC-phlex initial concentrations.")

        key_name = "init conc"
        call assert_msg(445070498, &
                species_data%get_real(key_name, this%init_conc(i_spec)), &
                "Missing 'init conc' for species '"//spec_name//" for "// &
                "PartMC-phlex initial concentrations.")
        
        this%init_conc_phlex_id(i_spec) = &
                this%phlex_core%chem_spec_data%gas_state_id(spec_name)
        call assert_msg(940200584, this%init_conc_phlex_id(i_spec).gt.0, &
                "Could not find species '"//spec_name//"' in PartMC-phlex.")
      
        call species_list%iter_next()
        i_spec = i_spec + 1
      end do

    else

      ! No gas-phase species intial concentrations specified
      allocate(this%init_conc_phlex_id(0))
      allocate(this%init_conc(0))

    end if

  end subroutine load_init_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get initial concentrations for the moch MONARCH model (for testing only)
  subroutine get_init_conc(this, MONARCH_conc)

    !> PartMC-phlex <-> MONARCH interface
    class(monarch_interface_t) :: this
    !> MONARCH species concentrations to update
    real, intent(inout) :: MONARCH_conc(:,:,:,:)

    integer(kind=i_kind) :: i_spec

    ! Reset the species concentrations in PMC and MONARCH
    this%phlex_state%state_var(:) = 0.0
    MONARCH_conc(:,:,:,:) = 0.0

    ! Set initial concentrations in PMC
    this%phlex_state%state_var(this%init_conc_phlex_id(:)) = &
            this%init_conc(:)

    ! Copy species concentrations to MONARCH array
    forall (i_spec = 1:size(this%map_monarch_id))
      MONARCH_conc(:,:,:,this%map_monarch_id(i_spec)) = &
              this%phlex_state%state_var(this%map_phlex_id(i_spec))
    end forall

  end subroutine get_init_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_monarch_interface
