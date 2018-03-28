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

!> The abstract aero_rep_single_particle_t structure and associated subroutines.
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

#define _RADIUS_ this%condensed_data_real(1)
#define _NUMBER_CONC_ this%condensed_data_real(2)
#define _NUM_INT_PROP_ 0
#define _NUM_REAL_PROP_ 2
  
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
    !! representation.
    !! 
    !! For a single particle representation, the unique names will be the
    !! phase name with the species name separated by a '.'
    procedure :: unique_names
    !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array by unique name. These are unique ids for each variable on the
    !! state array for this \ref phlex_aero_rep  "aerosol representation" and
    !! are numbered:
    !!
    !!   \f$x_u = x_f ... (x_f+n-1)\f$
    !!
    !! where \f$x_u\f$ is the id of the variable corresponding to unique name
    !! \f$u\f$ on the \c pmc_phlex_state::phlex_state_t::state_var array,
    !! \f$x_u\f$ is the index of the first variable for this aerosol
    !! representation on the state array and \f$n\f$ is the total number of
    !! variables on the state array from this aerosol representation.
    procedure :: state_id_by_unique_name
    !> Get a set of ids on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array for a particular aerosol species. These are unique ids for each
    !! variable on the state array that correspond to a particular species in
    !! this \ref phlex_aero_rep "aerosol representation" and are numbered:
    !!
    !!   \f$x_{si} \in x_f ... (x_f+n-1)\f$
    !!
    !! where \f$x_{si}\f$ is the index of species \f$s\f$ in phase instance
    !! \f$i\f$, \f$x_f\f$ is the index of the first variable for this aerosol
    !! representation on the \c pmc_phlex_state:phlex_state_t::state_var array
    !! and \f$n\f$ is the total number of variables on the state array for
    !! this aerosol representation. The size of the returned array will be:
    !!
    !!   \f$\text{size} = \sum_p n_p\f$
    !!
    !! where \f$n_p\f$ is the number of instances of phase \f$p\f$ in this
    !! aerosol representation. If species \f$s\f$ is not present in phase
    !! \f$p\f$ then \f$x_{si} \equiv 0 \quad \forall i \in p\f$.
    !! 
    !! This function should only be called during initialization
    procedure :: species_state_id
    !> Get the non-unique name of a species in this aerosol representation by
    !! id.
    procedure :: spec_name_by_id

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

    ! Assume all phases will be applied once to each particle
    allocate(this%aero_phase(size(aero_phase_set)))
    do i_phase = 1, size(aero_phase_set)
      allocate(this%aero_phase(i_phase)%val)
      this%aero_phase(i_phase)%val = aero_phase_set(i_phase)%val
    end do

    ! Get the total number of species across all phases
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
      num_spec = num_spec + this%aero_phase(i_phase)%val%size()
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(_NUM_INT_PROP_))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_))

    ! Set indexes
    allocate(this%phase_state_id(size(this%aero_phase)))
    curr_id = spec_state_id
    do i_phase = 1, size(this%aero_phase)
      this%phase_state_id(i_phase) = curr_id
      curr_id = curr_id + this%aero_phase(i_phase)%val%size()
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
  !! representation.
  !! 
  !! For a single particle representation, the unique names will be the
  !! phase name with the species name separated by a '.'
  function unique_names(this)

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    integer(kind=i_kind) :: num_spec, i_spec, j_spec, i_phase
    character(len=:), allocatable :: phase_name, spec_name
    
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
      num_spec = num_spec + this%aero_phase(i_phase)%val%size()
    end do
    allocate(unique_names(num_spec))
    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      phase_name = this%aero_phase(i_phase)%val%name()
      num_spec = this%aero_phase(i_phase)%val%size()
      do j_spec = 1, num_spec
        unique_names(i_spec + j_spec - 1)%string = &
                phase_name//'.'// &
                this%aero_phase(i_phase)%val%get_species_name(j_spec)
      end do
      i_spec = i_spec + num_spec
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by unique name. These are unique ids for each variable on the
  !! state array for this \ref phlex_aero_rep  "aerosol representation" and
  !! are numbered:
  !!
  !!   \f$x_u = x_f ... (x_f+n-1)\f$
  !!
  !! where \f$x_u\f$ is the id of the variable corresponding to unique name
  !! \f$u\f$ on the \c pmc_phlex_state::phlex_state_t::state_var array,
  !! \f$x_u\f$ is the index of the first variable for this aerosol
  !! representation on the state array and \f$n\f$ is the total number of
  !! variables on the state array from this aerosol representation.
  function state_id_by_unique_name(this, unique_name) result (spec_id)

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

  end function state_id_by_unique_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a set of ids on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array for a particular aerosol species. These are unique ids for each
  !! variable on the state array that correspond to a particular species in
  !! this \ref phlex_aero_rep "aerosol representation" and are numbered:
  !!
  !!   \f$x_{si} \in x_f ... (x_f+n-1)\f$
  !!
  !! where \f$x_{si}\f$ is the index of species \f$s\f$ in phase instance
  !! \f$i\f$, \f$x_f\f$ is the index of the first variable for this aerosol
  !! representation on the \c pmc_phlex_state:phlex_state_t::state_var array
  !! and \f$n\f$ is the total number of variables on the state array for
  !! this aerosol representation. The size of the returned array will be:
  !!
  !!   \f$\text{size} = \sum_p n_p\f$
  !!
  !! where \f$n_p\f$ is the number of instances of phase \f$p\f$ in this
  !! aerosol representation. If species \f$s\f$ is not present in phase
  !! \f$p\f$ then \f$x_{si} \equiv 0 \quad \forall i \in p\f$.
  !! 
  !! This function should only be called during initialization
  function species_state_id(this, phase_name, species_name) result(spec_index)
    use pmc_util,                                     only : i_kind

    !> Species index array
    integer(kind=i_kind), allocatable :: spec_index(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase id
    character(len=:), allocatable, intent(in) :: phase_name
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

    integer(kind=i_kind) :: i_phase, i_spec

    allocate(spec_index(1))
    i_phase = this%phase_id(phase_name)
    call assert_msg(820720954, i_phase.gt.0, "Invalid phase requested: "// &
            phase_name)
    i_spec = this%aero_phase(i_phase)%val%spec_id(species_name)
    call assert_msg(413945507, i_spec.gt.0, "Invalid species requested: "// &
            species_name//" for phase: "//phase_name)
    spec_index(1) = this%phase_state_id(i_phase) + i_spec - 1

  end function species_state_id

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

#undef _RADIUS_
#undef _NUMBER_CONC_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_

end module pmc_aero_rep_single_particle
