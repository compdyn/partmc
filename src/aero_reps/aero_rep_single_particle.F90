! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle module.

!> \page phlex_aero_rep_single_particle Phlexible Module for Chemistry: Single Particle Aerosol Representation
!!
!! The single particle aerosol representation is for use with a PartMC
!! particle-resolved run. The \c json object for this \ref phlex_aero_rep
!! "aerosol representation" is of the form:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "type" : "AERO_REP_SINGLE_PARTICLE"
!!    },
!!    ...
!!  ]}
!! \endcode
!! The only key-value pair required is the \b type, which must be \b
!! AERO_REP_SINGLE_PARTICLE. This representation assumes that every
!! \ref input_format_aero_phase "aerosol phase object" available will be
!! present once in each particle, and that integration of the
!! \ref input_format_mechanism "chemical mechanisms" at each time step will
!! solve the gas-phase first then do phase-transfer and aerosol-phase
!! chemistry for each single particle in the \c pmc_aero_particle_array::aero_particle_array_t variable
!! sequentially.

!> The aero_rep_single_particle_t type and associated subroutines.
module pmc_aero_rep_single_particle

  use pmc_util,                                  only: dp, i_kind, &
                                                       string_t, assert_msg, &
                                                       die_msg, to_string
  use pmc_property
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_phase_data
  use pmc_phlex_state

  implicit none
  private

#define _NUM_PHASE_ this%condensed_data_int(1)
#define _TOTAL_INT_PARAM_ this%condensed_data_int(2)
#define _RADIUS_ this%condensed_data_real(1)
#define _NUMBER_CONC_ this%condensed_data_real(2)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 2
#define _PHASE_INT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _NUM_SPEC_(x) this%condensed_data_int(_PHASE_INT_PARAM_LOC_(x))
#define _SPEC_STATE_ID_(x,y) this%condensed_data_int(_PHASE_INT_PARAM_LOC_(x)+x)
#define _AERO_PHASE_MASS_(x) this%condensed_data_real(_NUM_FLOAT_PROP_+x)
  
  ! Update types (These must match values in aero_rep_single_particle.c)
  integer(kind=i_kind), parameter, public :: UPDATE_RADIUS = 0
  integer(kind=i_kind), parameter, public :: UPDATE_NUMBER_CONC = 1

  public :: aero_rep_single_particle_t

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
    !! array by unique name. These are unique ids for each element on the
    !! state array for this \ref phlex_aero_rep  "aerosol representation" and
    !! are numbered:
    !!
    !!   \f$x_u = x_f ... (x_f+n-1)\f$
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to concentration
    !! of the species with unique name \f$u\f$ on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
    !! of the first element for this aerosol representation on the state array
    !! and \f$n\f$ is the total number of variables on the state array from
    !! this aerosol representation.
    procedure :: spec_state_id
    !> Get the non-unique name of a species in this aerosol representation by
    !! id.
    procedure :: spec_name_by_id
    !> Get the number of instances of an aerosol phase in this representation
    procedure :: num_phase_instances

  end type aero_rep_single_particle_t

  ! Constructor for aero_rep_single_particle_t
  interface aero_rep_single_particle_t
    procedure :: constructor
  end interface aero_rep_single_particle_t

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
  subroutine initialize(this, aero_phase_set, &
                  spec_state_id, chem_spec_data)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representationin the model species
    !! state array
    integer(kind=i_kind), intent(in) :: spec_state_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    integer(kind=i_kind) :: i_phase, i_spec, curr_id, curr_spec_id, num_spec
    type(property_t), pointer :: spec_props
    real(kind=dp) :: density
    character(len=:), allocatable :: key, spec_name
    integer(kind=i_kind) :: num_int_param, num_float_param

    ! Start off the counters
    num_int_param = _NUM_INT_PROP_
    num_float_param = _NUM_REAL_PROP_

    ! Assume all phases will be applied once to each particle
    allocate(this%aero_phase(size(aero_phase_set)))
    do i_phase = 1, size(aero_phase_set)
      allocate(this%aero_phase(i_phase)%val)
      this%aero_phase(i_phase)%val = aero_phase_set(i_phase)%val

      ! Add space for phase parameter location, number of species and
      ! species state ids
      num_int_param = num_int_param + 2 + &
              this%aero_phase(i_phase)%val%size()

      ! Add space for aerosol phase mass
      num_float_param = num_float_param + 1
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(num_int_param))
    allocate(this%condensed_data_real(num_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set the number of int parameters
    _TOTAL_INT_PARAM_ = num_int_param

    ! Set phase parameter locations, size and ids
    allocate(this%phase_state_id(size(this%aero_phase)))
    curr_id = spec_state_id
    num_int_param = _NUM_INT_PROP_ + 1
    do i_phase = 1, size(this%aero_phase)
      _PHASE_INT_PARAM_LOC_(i_phase) = num_int_param
      _NUM_SPEC_(i_phase) = this%aero_phase(i_phase)%val%size()
      this%phase_state_id(i_phase) = curr_id
      curr_id = curr_id + _NUM_SPEC_(i_phase)
      do i_spec = 1, _NUM_SPEC_(i_phase)
        _SPEC_STATE_ID_(i_phase, i_spec) = &
                this%phase_state_id(i_phase) + i_spec - 1
      end do
      num_int_param = num_int_param + 1 + _NUM_SPEC_(i_phase)
    end do

  end subroutine initialize

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
    character(len=:), allocatable :: curr_phase_name, curr_spec_name
    
    ! Count the number of unique names
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
      curr_phase_name = this%aero_phase(i_phase)%val%name()
      if (present(phase_name)) then
        if (phase_name.ne.curr_phase_name) cycle
      end if
      if (present(spec_name).or.present(tracer_type)) then
        do j_spec = 1, this%aero_phase(i_phase)%val%size()
          curr_spec_name = &
                  this%aero_phase(i_phase)%val%get_species_name(j_spec)
          curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type(j_spec)
          if (present(spec_name)) then
            if (spec_name.ne.curr_spec_name) cycle
          end if
          if (present(tracer_type)) then
            if (tracer_type.ne.curr_tracer_type) cycle
          end if
          num_spec = num_spec + 1
        end do
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
      num_spec = this%aero_phase(i_phase)%val%size()
      do j_spec = 1, num_spec
        curr_spec_name = &
                this%aero_phase(i_phase)%val%get_species_name(j_spec)
        curr_tracer_type = &
                this%aero_phase(i_phase)%val%get_species_type(j_spec)
        if (present(spec_name)) then
          if (spec_name.ne.curr_spec_name) cycle
        end if
        if (present(tracer_type)) then
          if (tracer_type.ne.curr_tracer_type) cycle
        end if
        unique_names(i_spec)%string = &
                curr_phase_name//'.'//curr_spec_name
        i_spec = i_spec + 1
      end do
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by unique name. These are unique ids for each element on the
  !! state array for this \ref phlex_aero_rep  "aerosol representation" and
  !! are numbered:
  !!
  !!   \f$x_u = x_f ... (x_f+n-1)\f$
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to concentration
  !! of the species with unique name \f$u\f$ on the \c
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
  function spec_name_by_id(this, aero_rep_spec_id)

    !> Chemical species name
    character(len=:), allocatable :: spec_name_by_id
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Indoex of species in this aerosol representation
    integer(kind=i_kind) :: aero_rep_spec_id

    ! Indices for iterators
    integer(kind=i_kind) :: i_spec, j_spec, i_phase
    
    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      do j_spec = 1, this%aero_phase(i_phase)%val%size()
        if (i_spec.eq.aero_rep_spec_id) then
          spec_name_by_id = this%aero_phase(i_phase)%val%get_species_name(j_spec)
        end if
        i_spec = i_spec + 1
      end do
    end do

  end function spec_name_by_id

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

#undef _NUM_PHASE_
#undef _TOTAL_INT_PARAM_
#undef _RADIUS_
#undef _NUMBER_CONC_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_INT_PARAM_LOC_
#undef _NUM_SPEC_
#undef _SPEC_STATE_ID_
#undef _AERO_PHASE_MASS_

end module pmc_aero_rep_single_particle
