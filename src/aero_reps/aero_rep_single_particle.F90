! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle module.

!> \page camp_aero_rep_single_particle CAMP: Single Particle Aerosol Representation
!!
!! The single particle aerosol representation is for use with a PartMC
!! particle-resolved run. The \c json object for this \ref camp_aero_rep
!! "aerosol representation" has the following format:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "name" : "my single particle aero rep",
!!      "type" : "AERO_REP_SINGLE_PARTICLE"
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be \b
!! AERO_REP_SINGLE_PARTICLE. This representation assumes that every \ref
!! input_format_aero_phase "aerosol phase" available will be present
!! once in each particle, and that the \ref input_format_mechanism
!! "chemical mechanisms" will be solved at each time step first for the
!! gas-phase then for phase-transfer and aerosol-phase chemistry for each
!! single particle in the \c  pmc_aero_particle_array::aero_particle_array_t
!! variable sequentially. This may be changed in the future to solve for all
!! particles simultaneously.
!!
!! The number concentration for each particle must be
!! set from an external model using
!! \c pmc_aero_rep_single_particle::aero_rep_update_data_single_particle_number_t
!! objects.

!> The aero_rep_single_particle_t type and associated subroutines.
module pmc_aero_rep_single_particle

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_camp_state
  use pmc_mpi
  use pmc_property
  use pmc_util,                                  only: dp, i_kind, &
                                                       string_t, assert_msg, &
                                                       die_msg, to_string, &
                                                       assert

  use iso_c_binding

  implicit none
  private

#define NUM_PHASE_ this%condensed_data_int(1)
#define AERO_REP_ID_ this%condensed_data_int(2)
#define MAX_PARTICLES_ this%condensed_data_int(3)
#define PARTICLE_STATE_SIZE_ this%condensed_data_int(4)
#define NUM_INT_PROP_ 4
#define NUM_REAL_PROP_ 0
#define NUM_ENV_PARAM_PER_PARTICLE_ 1
#define PHASE_STATE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PHASE_MODEL_DATA_ID_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+x)
#define PHASE_NUM_JAC_ELEM_(x) this%condensed_data_int(NUM_INT_PROP_+2*NUM_PHASE_+x)

  ! Update types (These must match values in aero_rep_single_particle.c)
  integer(kind=i_kind), parameter, public :: UPDATE_NUMBER_CONC = 0

  public :: aero_rep_single_particle_t, &
            aero_rep_update_data_single_particle_number_t

  !> Single particle aerosol representation
  !!
  !! Time-invariant data related to a single particle aerosol representation.
  type, extends(aero_rep_data_t) :: aero_rep_single_particle_t
    !> Unique names for each instance of every chemical species in the
    !! aerosol representaiton
    type(string_t), allocatable, private :: unique_names_(:)
    !> First state id for the representation (only used during initialization)
    integer(kind=i_kind) :: state_id_start = -99999
  contains
    !> Initialize the aerosol representation data, validating component data and
    !! loading any required information from the \c
    !! aero_rep_data_t::property_set. This routine should be called once for
    !! each aerosol representation at the beginning of a model run after all
    !! the input files have been read in. It ensures all data required during
    !! the model run are included in the condensed data arrays.
    procedure :: initialize
    !> Returns the maximum number of computational particles
    procedure :: maximum_computational_particles
    !> Initialize an update data number object
    procedure :: update_data_initialize_number => update_data_init_number
    !> Get the size of the section of the
    !! \c pmc_camp_state::camp_state_t::state_var array required for this
    !! aerosol representation.
    !!
    !! For a single particle representation, the size will correspond to the
    !! the sum of the sizes of a single instance of each aerosol phase
    !! provided to \c aero_rep_single_particle::initialize()
    procedure :: size => get_size
    !> Get the number of state variables per-particle
    !!
    !! Calling functions can assume each particle has the same size on the
    !! state array, and that individual particle states are contiguous and
    !! arranged sequentially
    procedure :: per_particle_size
    !> Get a list of unique names for each element on the
    !! \c pmc_camp_state::camp_state_t::state_var array for this aerosol
    !! representation. The list may be restricted to a particular phase and/or
    !! aerosol species by including the phase_name and spec_name arguments.
    !!
    !! For a single particle representation, the unique names will be the
    !! phase name with the species name separated by a '.'
    procedure :: unique_names
    !> Get a species id on the \c pmc_camp_state::camp_state_t::state_var
    !! array by its unique name. These are unique ids for each element on the
    !! state array for this \ref camp_aero_rep "aerosol representation" and
    !! are numbered:
    !!
    !!   \f[x_u \in x_f ... (x_f+n-1)\f]
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to the species
    !! with unique name \f$u\f$ on the \c
    !! pmc_camp_state::camp_state_t::state_var array, \f$x_f\f$ is the index
    !! of the first element for this aerosol representation on the state array
    !! and \f$n\f$ is the total number of variables on the state array from
    !! this aerosol representation.
    procedure :: spec_state_id
    !> Get the non-unique name of a species by its unique name
    procedure :: spec_name
    !> Get the number of instances of an aerosol phase in this representation
    procedure :: num_phase_instances
    !> Get the number of Jacobian elements used in calculations of aerosol
    !! mass, volume, number, etc. for a particular phase
    procedure :: num_jac_elem
    !> Finalize the aerosol representation
    final :: finalize

  end type aero_rep_single_particle_t

  ! Constructor for aero_rep_single_particle_t
  interface aero_rep_single_particle_t
    procedure :: constructor
  end interface aero_rep_single_particle_t

  !> Single particle update number concentration object
  type, extends(aero_rep_update_data_t) :: &
            aero_rep_update_data_single_particle_number_t
  private
    !> Flag indicating whether the update data is allocated
    logical :: is_malloced = .false.
    !> Unique id for finding aerosol representations during initialization
    integer(kind=i_kind) :: aero_rep_unique_id = 0
    !> Maximum number of computational particles
    integer(kind=i_kind) :: maximum_computational_particles = 0
  contains
    !> Update the number
    procedure :: set_number__n_m3 => update_data_set_number__n_m3
    !> Determine the pack size of the local update data
    procedure :: internal_pack_size => internal_pack_size_number
    !> Pack the local update data to a binary
    procedure :: internal_bin_pack => internal_bin_pack_number
    !> Unpack the local update data from a binary
    procedure :: internal_bin_unpack => internal_bin_unpack_number
    !> Finalize the number update data
    final :: update_data_number_finalize
  end type aero_rep_update_data_single_particle_number_t

  !> Interface to c aerosol representation functions
  interface

    !> Allocate space for a number update
    function aero_rep_single_particle_create_number_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function aero_rep_single_particle_create_number_update_data

    !> Set a new particle number concentration
    subroutine aero_rep_single_particle_set_number_update_data__n_m3( &
              update_data, aero_rep_unique_id, particle_id, number_conc) &
              bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Aerosol representation unique id
      integer(kind=c_int), value :: aero_rep_unique_id
      !> Computational particle index
      integer(kind=c_int), value :: particle_id
      !> New number (m)
      real(kind=c_double), value :: number_conc
    end subroutine aero_rep_single_particle_set_number_update_data__n_m3

    !> Free an update data object
    pure subroutine aero_rep_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine aero_rep_free_update_data

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_single_particle_t
  function constructor() result (new_obj)

    !> New aerosol representation
    type(aero_rep_single_particle_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the \c
  !! aero_rep_data_t::property_set. This routine should be called once for
  !! each aerosol representation at the beginning of a model run after all
  !! the input files have been read in. It ensures all data required during
  !! the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_phase_set, spec_state_id)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representationin the model species
    !! state array
    integer(kind=i_kind), intent(in) :: spec_state_id

    character(len=:), allocatable :: key_name
    integer(kind=i_kind) :: i_particle, i_phase, curr_id
    integer(kind=i_kind) :: num_int_param, num_float_param, num_particles

    ! Start off the counters
    num_int_param = NUM_INT_PROP_ + 3*size(aero_phase_set)
    num_float_param = NUM_REAL_PROP_

    ! Get the maximum number of computational particles
    key_name = "maximum computational particles"
    call assert_msg(857908074, &
                    this%property_set%get_int(key_name, num_particles), &
                    "Missing maximum number of computational particles")

    ! Assume all phases will be applied once to each particle
    allocate(this%aero_phase(size(aero_phase_set)*num_particles))
    do i_particle = 1, num_particles
      do i_phase = 1, size(aero_phase_set)
        this%aero_phase((i_particle-1)*size(aero_phase_set)+i_phase) = &
            aero_phase_set(i_phase)
      end do
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(num_int_param))
    allocate(this%condensed_data_real(num_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save the number of computational particles
    MAX_PARTICLES_ = num_particles

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_PER_PARTICLE_ * num_particles

    ! Set phase state and model data ids
    NUM_PHASE_ = size(aero_phase_set)
    this%state_id_start = spec_state_id
    curr_id = spec_state_id
    do i_phase = 1, NUM_PHASE_
      PHASE_STATE_ID_(i_phase) = curr_id
      PHASE_MODEL_DATA_ID_(i_phase) = i_phase
      curr_id = curr_id + aero_phase_set(i_phase)%val%size()
    end do
    PARTICLE_STATE_SIZE_ = curr_id - spec_state_id

    ! Initialize the aerosol representation id
    AERO_REP_ID_ = -1

    ! Set the unique names for the chemical species
    this%unique_names_ = this%unique_names( )

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the maximum nunmber of computational particles
  integer(kind=i_kind) function maximum_computational_particles(this)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    maximum_computational_particles = MAX_PARTICLES_

  end function maximum_computational_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the section of the
  !! \c pmc_camp_state::camp_state_t::state_var array required for this
  !! aerosol representation.
  !!
  !! For a single particle representation, the size will correspond to the
  !! the sum of the sizes of a single instance of each aerosol phase
  !! provided to \c aero_rep_single_particle::initialize()
  function get_size(this) result (state_size)

    !> Size on the state array
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    state_size = MAX_PARTICLES_ * PARTICLE_STATE_SIZE_

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of state variables per-particle
  !!
  !! Calling functions can assume each particle has the same size on the
  !! state array, and that individual particle states are contiguous and
  !! arranged sequentially
  function per_particle_size(this) result(state_size)

    !> Size on the state array per particle
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    state_size = PARTICLE_STATE_SIZE_

  end function per_particle_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element on the
  !! \c pmc_camp_state::camp_state_t::state_var array for this aerosol
  !! representation. The list may be restricted to a particular phase and/or
  !! aerosol species by including the phase_name and spec_name arguments.
  !!
  !! For a single particle representation, the unique names will be a 'P'
  !! followed by the computational particle number, a '.', the phase name,
  !! another '.', and the species name.
  function unique_names(this, phase_name, tracer_type, spec_name)

    use pmc_util,                      only : integer_to_string

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase name
    character(len=*), optional, intent(in) :: phase_name
    !> Aerosol-phase species tracer type
    integer(kind=i_kind), optional, intent(in) :: tracer_type
    !> Aerosol-phase species name
    character(len=*), optional, intent(in) :: spec_name

    integer(kind=i_kind) :: num_spec, i_part, i_spec, j_spec, i_phase
    integer(kind=i_kind) :: curr_tracer_type
    character(len=:), allocatable :: curr_phase_name
    type(string_t), allocatable :: spec_names(:)

    ! Copy saved unique names when no filters are included
    if (.not. present(phase_name) .and. &
        .not. present(tracer_type) .and. &
        .not. present(spec_name) .and. &
        allocated(this%unique_names_)) then
      unique_names = this%unique_names_
      return
    end if

    ! Count the number of unique names
    num_spec = 0
    do i_phase = 1, NUM_PHASE_
      curr_phase_name = this%aero_phase(i_phase)%val%name()
      if (present(phase_name)) then
        if (phase_name.ne.curr_phase_name) cycle
      end if
      if (present(spec_name).or.present(tracer_type)) then
        spec_names = this%aero_phase(i_phase)%val%get_species_names()
        do j_spec = 1, size(spec_names)
          curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type( &
                  spec_names(j_spec)%string)
          if (present(spec_name)) then
            if (spec_name.ne.spec_names(j_spec)%string) cycle
          end if
          if (present(tracer_type)) then
            if (tracer_type.ne.curr_tracer_type) cycle
          end if
          num_spec = num_spec + 1
        end do
        deallocate(spec_names)
      else
        num_spec = num_spec + this%aero_phase(i_phase)%val%size()
      end if
    end do

    ! Allocate space for the unique names and assign them
    allocate(unique_names(num_spec*MAX_PARTICLES_))
    i_spec = 1
    do i_part = 1, MAX_PARTICLES_
      do i_phase = 1, NUM_PHASE_
        curr_phase_name = this%aero_phase(i_phase)%val%name()
        if (present(phase_name)) then
          if (phase_name.ne.curr_phase_name) cycle
        end if
        spec_names = this%aero_phase(i_phase)%val%get_species_names()
        num_spec = this%aero_phase(i_phase)%val%size()
        do j_spec = 1, num_spec
          curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type( &
                  spec_names(j_spec)%string)
          if (present(spec_name)) then
            if (spec_name.ne.spec_names(j_spec)%string) cycle
          end if
          if (present(tracer_type)) then
            if (tracer_type.ne.curr_tracer_type) cycle
          end if
          unique_names(i_spec)%string = 'P'//trim(integer_to_string(i_part))//&
                  "."//curr_phase_name//'.'//spec_names(j_spec)%string
          i_spec = i_spec + 1
        end do
      end do
      deallocate(spec_names)
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_camp_state::camp_state_t::state_var
  !! array by its unique name. These are unique ids for each element on the
  !! state array for this \ref camp_aero_rep "aerosol representation" and
  !! are numbered:
  !!
  !!   \f[x_u \in x_f ... (x_f+n-1)\f]
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to the species
  !! with unique name \f$u\f$ on the \c
  !! pmc_camp_state::camp_state_t::state_var array, \f$x_f\f$ is the index
  !! of the first element for this aerosol representation on the state array
  !! and \f$n\f$ is the total number of variables on the state array from
  !! this aerosol representation.
  function spec_state_id(this, unique_name) result (spec_id)

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Unique name
    character(len=*), intent(in) :: unique_name

    integer(kind=i_kind) :: i_spec

    spec_id = this%state_id_start
    do i_spec = 1, size(this%unique_names_)
      if (this%unique_names_(i_spec)%string .eq. unique_name) then
        return
      end if
      spec_id = spec_id + 1
    end do
    call die_msg( 449087541, "Cannot find species '"//unique_name//"'" )

  end function spec_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the non-unique name of a species in this aerosol representation by
  !! id.
  function spec_name(this, unique_name)

    !> Chemical species name
    character(len=:), allocatable :: spec_name
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Unique name of the species in this aerosol representation
    character(len=*), intent(in) :: unique_name

    ! Indices for iterators
    integer(kind=i_kind) :: i_part, i_spec, j_spec, i_phase

    ! Species in aerosol phase
    type(string_t), allocatable :: spec_names(:)

    call assert( 124916561, allocated( this%unique_names_ ) )
    i_spec = 1
    do i_phase = 1, size(this%aero_phase) ! each phase in each partice
      spec_names = this%aero_phase(i_phase)%val%get_species_names()
      do j_spec = 1, size(spec_names)
        if (unique_name .eq. this%unique_names_(i_spec)%string) then
          spec_name = spec_names(j_spec)%string
          return
        end if
        i_spec = i_spec + 1
      end do
    end do
    call die_msg(101731871, "Could not find unique name '"//unique_name//"'")

  end function spec_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of instances of a specified aerosol phase. In the single
  !! particle representation, if an aerosol phase is present, is exists once
  !! in each computational particle.
  function num_phase_instances(this, phase_name)

    !> Number of instances of the aerosol phase
    integer(kind=i_kind) :: num_phase_instances
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase name
    character(len=*), intent(in) :: phase_name

    integer(kind=i_kind) :: i_phase

    num_phase_instances = 0
    do i_phase = 1, NUM_PHASE_
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        num_phase_instances = MAX_PARTICLES_
        return
      end if
    end do

  end function num_phase_instances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of Jacobian elements used in calculations of aerosol mass,
  !! volume, number, etc. for a particular phase
  function num_jac_elem(this, phase_id)

    !> Number of Jacobian elements used
    integer(kind=i_kind) :: num_jac_elem
    !> Aerosol respresentation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase id
    integer(kind=i_kind), intent(in) :: phase_id

    integer(kind=i_kind) :: i_phase

    call assert_msg( 401502046, phase_id .ge. 1 .and. &
                                phase_id .le. size( this%aero_phase ), &
                     "Aerosol phase index out of range. Got "// &
                     trim( integer_to_string( phase_id ) )//", expected 1:"// &
                     trim( integer_to_string( size( this%aero_phase ) ) ) )
    num_jac_elem = 0
    do i_phase = 1, NUM_PHASE_
      num_jac_elem = num_jac_elem + &
                     this%aero_phase(i_phase)%val%num_jac_elem()
    end do

  end function num_jac_elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the aerosol representation
  elemental subroutine finalize(this)

    !> Aerosol representation data
    type(aero_rep_single_particle_t), intent(inout) :: this

    if (allocated(this%rep_name)) deallocate(this%rep_name)
    if (allocated(this%aero_phase)) then
      ! The core will deallocate the aerosol phases
      call this%aero_phase(:)%dereference()
      deallocate(this%aero_phase)
    end if
    if (associated(this%property_set)) deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize an update data object
  subroutine update_data_init_number(this, update_data, aero_rep_type)

    use pmc_rand,                                only : generate_int_id

    !> Aerosol representation to update
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> Update data object
    class(aero_rep_update_data_single_particle_number_t), intent(out) :: &
        update_data
    !> Aerosol representaiton id
    integer(kind=i_kind), intent(in) :: aero_rep_type

    ! If an aerosol representation id has not been generated, do it now
    if (AERO_REP_ID_.eq.-1) then
      AERO_REP_ID_ = generate_int_id()
    end if

    update_data%aero_rep_unique_id = AERO_REP_ID_
    update_data%maximum_computational_particles = &
        this%maximum_computational_particles( )
    update_data%aero_rep_type = int(aero_rep_type, kind=c_int)
    update_data%update_data = &
      aero_rep_single_particle_create_number_update_data()
    update_data%is_malloced = .true.

  end subroutine update_data_init_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for particle number (#/m3) for a particular
  !! computational particle.
  subroutine update_data_set_number__n_m3(this, particle_id, number_conc)

    !> Update data
    class(aero_rep_update_data_single_particle_number_t), intent(inout) :: &
            this
    !> Computational particle index
    integer(kind=i_kind), intent(in) :: particle_id
    !> Updated number
    real(kind=dp), intent(in) :: number_conc

    call assert_msg(897092373, this%is_malloced, &
            "Trying to set number of uninitialized update object.")
    call assert_msg(357138177, particle_id .ge. 1 .and. &
                    particle_id .le. this%maximum_computational_particles, &
                    "Invalid computational particle index: "// &
                    trim(integer_to_string(particle_id)))
    call aero_rep_single_particle_set_number_update_data__n_m3( &
            this%get_data(), this%aero_rep_unique_id, particle_id-1, &
            number_conc)

  end subroutine update_data_set_number__n_m3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the reaction data
  integer(kind=i_kind) function internal_pack_size_number(this, comm) &
      result(pack_size)

    !> Aerosol representation update data
    class(aero_rep_update_data_single_particle_number_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size = &
      pmc_mpi_pack_size_logical(this%is_malloced, comm) + &
      pmc_mpi_pack_size_integer(this%maximum_computational_particles, comm) + &
      pmc_mpi_pack_size_integer(this%aero_rep_unique_id, comm)

  end function internal_pack_size_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine internal_bin_pack_number(this, buffer, pos, comm)

    !> Aerosol representation update data
    class(aero_rep_update_data_single_particle_number_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_logical(buffer, pos, this%is_malloced, comm)
    call pmc_mpi_pack_integer(buffer, pos, &
                              this%maximum_computational_particles, comm)
    call pmc_mpi_pack_integer(buffer, pos, this%aero_rep_unique_id, comm)
    call assert(964639022, &
         pos - prev_position <= this%pack_size(comm))
#endif

  end subroutine internal_bin_pack_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine internal_bin_unpack_number(this, buffer, pos, comm)

    !> Aerosol representation update data
    class(aero_rep_update_data_single_particle_number_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_logical(buffer, pos, this%is_malloced, comm)
    call pmc_mpi_unpack_integer(buffer, pos, &
                                this%maximum_computational_particles, comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%aero_rep_unique_id, comm)
    call assert(459432617, &
         pos - prev_position <= this%pack_size(comm))
    this%update_data = aero_rep_single_particle_create_number_update_data()
#endif

  end subroutine internal_bin_unpack_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a number update data object
  elemental subroutine update_data_number_finalize(this)

    !> Update data object to free
    type(aero_rep_update_data_single_particle_number_t), intent(inout) :: this

    if (this%is_malloced) call aero_rep_free_update_data(this%update_data)

  end subroutine update_data_number_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_single_particle
