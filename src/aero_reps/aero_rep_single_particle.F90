! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle module.

!> \page phlex_aero_rep_single_particle Phlexible Module for Chemistry: Single Particle Aerosol Representation
!!
!! The single particle aerosol representation is for use with a PartMC
!! particle-resolved run. The \c json object for this \ref phlex_aero_rep
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

!> The aero_rep_single_particle_t type and associated subroutines.
module pmc_aero_rep_single_particle

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_phlex_state
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
#define RADIUS_ this%condensed_data_real(1)
#define NUMBER_CONC_ this%condensed_data_real(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 2
#define PHASE_STATE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PHASE_MODEL_DATA_ID_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+x)
#define PHASE_NUM_JAC_ELEM_(x) this%condensed_data_int(NUM_INT_PROP_+2*NUM_PHASE_+x)
#define PHASE_MASS_(x) this%condensed_data_real(NUM_REAL_PROP_+x)
#define PHASE_AVG_MW_(x) this%condensed_data_real(NUM_REAL_PROP_+NUM_PHASE_+x)

  ! Update types (These must match values in aero_rep_single_particle.c)
  integer(kind=i_kind), parameter, public :: UPDATE_RADIUS = 0
  integer(kind=i_kind), parameter, public :: UPDATE_NUMBER_CONC = 1

  public :: aero_rep_single_particle_t, &
            aero_rep_update_data_single_particle_radius_t, &
            aero_rep_update_data_single_particle_number_t

  !> Single particle aerosol representation
  !!
  !! Time-invariant data related to a single particle aerosol representation.
  type, extends(aero_rep_data_t) :: aero_rep_single_particle_t
    !> Phase state id (only used during initialization
    integer(kind=i_kind), allocatable :: phase_state_id(:)
  contains
    !> Initialize the aerosol representation data, validating component data and
    !! loading any required information from the \c
    !! aero_rep_data_t::property_set. This routine should be called once for
    !! each aerosol representation at the beginning of a model run after all
    !! the input files have been read in. It ensures all data required during
    !! the model run are included in the condensed data arrays.
    procedure :: initialize
    !> Set an id for this aerosol representation for use with updates from
    !! external modules
    procedure :: set_id
    !> Get the size of the section of the
    !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
    !! aerosol representation.
    !!
    !! For a single particle representation, the size will correspond to the
    !! the sum of the sizes of a single instance of each aerosol phase
    !! provided to \c aero_rep_single_particle::initialize()
    procedure :: size => get_size
    !> Get a list of unique names for each element on the
    !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
    !! representation. The list may be restricted to a particular phase and/or
    !! aerosol species by including the phase_name and spec_name arguments.
    !!
    !! For a single particle representation, the unique names will be the
    !! phase name with the species name separated by a '.'
    procedure :: unique_names
    !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array by its unique name. These are unique ids for each element on the
    !! state array for this \ref phlex_aero_rep "aerosol representation" and
    !! are numbered:
    !!
    !!   \f[x_u \in x_f ... (x_f+n-1)\f]
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to the species
    !! with unique name \f$u\f$ on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
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

  !> Single particle update radius object
  type, extends(aero_rep_update_data_t) :: &
            aero_rep_update_data_single_particle_radius_t
  private
    logical :: is_malloced = .false.
  contains
    !> Initialize the update data
    procedure :: initialize => update_data_init_radius
    !> Update the radius
    procedure :: set_radius => update_data_set_radius
    !> Finalize the radius update data
    final :: update_data_radius_finalize
  end type aero_rep_update_data_single_particle_radius_t

  !> Single particle update number concentration object
  type, extends(aero_rep_update_data_t) :: &
            aero_rep_update_data_single_particle_number_t
  private
    logical :: is_malloced = .false.
  contains
    !> Initialize the update data
    procedure :: initialize => update_data_init_number
    !> Update the number
    procedure :: set_number => update_data_set_number
    !> Finalize the number update data
    final :: update_data_number_finalize
  end type aero_rep_update_data_single_particle_number_t

  !> Interface to c aerosol representation functions
  interface

    !> Allocate space for a radius update
    function aero_rep_single_particle_create_radius_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function aero_rep_single_particle_create_radius_update_data

    !> Set a new particle radius
    subroutine aero_rep_single_particle_set_radius_update_data(update_data, &
              aero_rep_id, radius) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Aerosol representation id from
      !! pmc_aero_rep_single_particle::aero_rep_single_particle_t::set_id
      integer(kind=c_int), value :: aero_rep_id
      !> New radius (m)
      real(kind=c_double), value :: radius
    end subroutine aero_rep_single_particle_set_radius_update_data

    !> Allocate space for a number update
    function aero_rep_single_particle_create_number_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function aero_rep_single_particle_create_number_update_data

    !> Set a new particle number concentration
    subroutine aero_rep_single_particle_set_number_update_data(update_data, &
              aero_rep_id, number_conc) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Aerosol representation id from
      !! pmc_aero_rep_single_particle::aero_rep_single_particle_t::set_id
      integer(kind=c_int), value :: aero_rep_id
      !> New number (m)
      real(kind=c_double), value :: number_conc
    end subroutine aero_rep_single_particle_set_number_update_data

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

    integer(kind=i_kind) :: i_phase, curr_id
    integer(kind=i_kind) :: num_int_param, num_float_param

    ! Start off the counters
    num_int_param = NUM_INT_PROP_ + 3*size(aero_phase_set)
    num_float_param = NUM_REAL_PROP_ + 2*size(aero_phase_set)

    ! Assume all phases will be applied once to each particle
    allocate(this%aero_phase(size(aero_phase_set)))
    do i_phase = 1, size(aero_phase_set)
      this%aero_phase(i_phase) = aero_phase_set(i_phase)
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(num_int_param))
    allocate(this%condensed_data_real(num_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set phase state and model data ids
    NUM_PHASE_ = size(this%aero_phase)
    allocate(this%phase_state_id(NUM_PHASE_))
    curr_id = spec_state_id
    do i_phase = 1, NUM_PHASE_
      this%phase_state_id(i_phase) = curr_id
      curr_id = curr_id + aero_phase_set(i_phase)%val%size()
      PHASE_STATE_ID_(i_phase) = this%phase_state_id(i_phase)
      PHASE_MODEL_DATA_ID_(i_phase) = i_phase
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this aerosol representation that can be used by external
  !! modules to update the particle radius and number concentration
  subroutine set_id(this, new_id)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> New id
    integer(kind=i_kind), intent(in) :: new_id

    AERO_REP_ID_ = new_id

  end subroutine set_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the section of the
  !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
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

    integer(kind=i_kind) :: i_phase

    ! Get the total number of species across all phases
    state_size = 0
    do i_phase = 1, size(this%aero_phase)
      state_size = state_size + this%aero_phase(i_phase)%val%size()
    end do

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element on the
  !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
  !! representation. The list may be restricted to a particular phase and/or
  !! aerosol species by including the phase_name and spec_name arguments.
  !!
  !! For a single particle representation, the unique names will be the
  !! phase name with the species name separated by a '.'
  function unique_names(this, phase_name, tracer_type, spec_name)

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, optional, intent(in) :: phase_name
    !> Aerosol-phase species tracer type
    integer(kind=i_kind), optional, intent(in) :: tracer_type
    !> Aerosol-phase species name
    character(len=:), allocatable, optional, intent(in) :: spec_name

    integer(kind=i_kind) :: num_spec, i_spec, j_spec, i_phase
    integer(kind=i_kind) :: curr_tracer_type
    character(len=:), allocatable :: curr_phase_name
    type(string_t), allocatable :: spec_names(:)

    ! Count the number of unique names
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
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
    allocate(unique_names(num_spec))
    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
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
        unique_names(i_spec)%string = &
                curr_phase_name//'.'//spec_names(j_spec)%string
        i_spec = i_spec + 1
      end do
      deallocate(spec_names)
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by its unique name. These are unique ids for each element on the
  !! state array for this \ref phlex_aero_rep "aerosol representation" and
  !! are numbered:
  !!
  !!   \f[x_u \in x_f ... (x_f+n-1)\f]
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to the species
  !! with unique name \f$u\f$ on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
  !! of the first element for this aerosol representation on the state array
  !! and \f$n\f$ is the total number of variables on the state array from
  !! this aerosol representation.
  function spec_state_id(this, unique_name) result (spec_id)

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Unique name
    character(len=:), allocatable, intent(in) :: unique_name

    type(string_t), allocatable :: unique_names(:)
    integer(kind=i_kind) :: i_spec

    spec_id = 0
    unique_names = this%unique_names()
    do i_spec = 1, size(unique_names)
      if (unique_names(i_spec)%string .eq. unique_name) then
        spec_id = this%phase_state_id(1) + i_spec - 1
        return
      end if
    end do

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
    integer(kind=i_kind) :: i_spec, j_spec, i_phase

    ! Species in aerosol phase
    type(string_t), allocatable :: spec_names(:)

    ! Unique name list
    type(string_t), allocatable :: unique_names(:)

    unique_names = this%unique_names()

    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      spec_names = this%aero_phase(i_phase)%val%get_species_names()
      do j_spec = 1, size(spec_names)
        if (unique_name.eq.unique_names(i_spec)%string) then
          spec_name = spec_names(j_spec)%string
        end if
        i_spec = i_spec + 1
      end do
      deallocate(spec_names)
    end do

    deallocate(unique_names)

  end function spec_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of instances of a specified aerosol phase. In the single
  !! particle representation, if an aerosol phase is present, it has only
  !! one instance.
  function num_phase_instances(this, phase_name)

    !> Number of instances of the aerosol phase
    integer(kind=i_kind) :: num_phase_instances
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, intent(in) :: phase_name

    integer(kind=i_kind) :: i_phase

    num_phase_instances = 0
    do i_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        num_phase_instances = 1
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

    num_jac_elem = 0
    do i_phase = 1, size(this%aero_phase)
      num_jac_elem = num_jac_elem +                                          &
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

  !> Initialize radius update data object
  subroutine update_data_init_radius(this, aero_rep_type)

    !> Update data object
    class(aero_rep_update_data_single_particle_radius_t) :: this
    !> Aerosol representaiton id
    integer(kind=i_kind), intent(in) :: aero_rep_type

    call assert_msg(954504234, .not.this%is_malloced, &
            "Trying to re-initializea radius update object.")

    this%aero_rep_type = int(aero_rep_type, kind=c_int)
    this%update_data = aero_rep_single_particle_create_radius_update_data()
    this%is_malloced = .true.

  end subroutine update_data_init_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for particle radius
  subroutine update_data_set_radius(this, aero_rep_id, radius)

    !> Update data
    class(aero_rep_update_data_single_particle_radius_t), intent(inout) :: &
            this
    !> Aerosol representation id from
    !! pmc_aero_rep_single_particle::aero_rep_single_particle_t::set_id
    integer(kind=i_kind), intent(in) :: aero_rep_id
    !> Updated radius
    real(kind=dp), intent(in) :: radius

    call assert_msg(228446610, this%is_malloced, &
            "Trying to set radius of uninitialized update object.")
    call aero_rep_single_particle_set_radius_update_data(this%get_data(), &
            aero_rep_id, radius)

  end subroutine update_data_set_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a radius update data object
  elemental subroutine update_data_radius_finalize(this)

    !> Update data object to free
    type(aero_rep_update_data_single_particle_radius_t), intent(inout) :: this

    if (this%is_malloced) call aero_rep_free_update_data(this%update_data)

  end subroutine update_data_radius_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize an update data object
  subroutine update_data_init_number(this, aero_rep_type)

    !> Update data object
    class(aero_rep_update_data_single_particle_number_t) :: this
    !> Aerosol representaiton id
    integer(kind=i_kind), intent(in) :: aero_rep_type

    call assert_msg(886589356, .not.this%is_malloced, &
            "Trying to re-initializea number update object.")

    this%aero_rep_type = int(aero_rep_type, kind=c_int)
    this%update_data = aero_rep_single_particle_create_number_update_data()
    this%is_malloced = .true.

  end subroutine update_data_init_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for particle number
  subroutine update_data_set_number(this, aero_rep_id, number_conc)

    !> Update data
    class(aero_rep_update_data_single_particle_number_t), intent(inout) :: &
            this
    !> Aerosol representation id from
    !! pmc_aero_rep_single_particle::aero_rep_single_particle_t::set_id
    integer(kind=i_kind), intent(in) :: aero_rep_id
    !> Updated number
    real(kind=dp), intent(in) :: number_conc

    call assert_msg(897092373, this%is_malloced, &
            "Trying to set number of uninitialized update object.")
    call aero_rep_single_particle_set_number_update_data(this%get_data(), &
            aero_rep_id, number_conc)

  end subroutine update_data_set_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a number update data object
  elemental subroutine update_data_number_finalize(this)

    !> Update data object to free
    type(aero_rep_update_data_single_particle_number_t), intent(inout) :: this

    if (this%is_malloced) call aero_rep_free_update_data(this%update_data)

  end subroutine update_data_number_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef NUM_PHASE_
#undef AERO_REP_ID_
#undef RADIUS_
#undef NUMBER_CONC_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef PHASE_STATE_ID_
#undef PHASE_MODEL_DATA_ID_
#undef PHASE_MASS_
#undef PHASE_AVG_MW_

end module pmc_aero_rep_single_particle
