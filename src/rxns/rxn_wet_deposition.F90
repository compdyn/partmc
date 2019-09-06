! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_wet_deposition module.

!> \page phlex_rxn_wet_deposition Phlexible Module for Chemistry: Wet Deposition
!!
!! Wet Deposition reactions take the form:
!!
!! \f[
!!   \mbox{X} \rightarrow
!! \f]
!!
!! where \f$\mbox{\ch{X}}\f$ is a set of species in an aerosol phase
!! undergoing wet deposition at a given rate.
!!
!! Wet deposition rate constants can be constant or set from an external
!! module using the
!! \c pmc_rxn_wet_deposition::rxn_update_data_wet_deposition_rate_t object.
!! External modules can use the
!! \c pmc_rxn_wet_deposition::rxn_wet_deposition_t::get_property_set()
!! function during initilialization to access any needed reaction parameters
!! to identify certain wet deposition reactions. The
!! \c pmc_rxn_wet_deposition::rxn_wet_deposition_t::set_rxn_id() function
!! can be used during initialization to set an integer id for a particular
!! reaction that can be used during solving to update the wet deposition
!! rate from an external module.
!!
!! Input data for wet deposition reactions have the following format :
!! \code{.json}
!!   {
!!     "type" : "WET_DEPOSITION",
!!     "aerosol phase" : "my aero phase",
!!     "rate const" : 12.5,
!!     "scaling factor" : 1.2,
!!     ...
!!   }
!! \endcode
!! The key-value pair \b aerosol \b phase is required and its value must be
!! the name of the aerosol phase undergoing wet deposition. All species within
!! the aerosol phase in all instances of the aerosol phase will be removed
!! according the first-order loss rate constant. The key-value pair \b rate
!! \b const is optional and can be used to set a rate constant that remains
!! constant throughout the model run. The \b scaling \b factor is also
!! optional, and can be used to set a constant scaling factor for the rate
!! constant. When a \b scaling \b factor is not provided, it is assumed to be
!! 1.0. All other data is optional and will be available to external modules
!! during initialization. Rate constants are in units of \f$s^{-1}\f$.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_wet_deposition_t type and associated functions.
module pmc_rxn_wet_deposition

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, string_t, &
                                                  to_string, assert, &
                                                  assert_msg, die_msg

  use iso_c_binding

  implicit none
  private

#define RXN_ID_ this%condensed_data_int(1)
#define NUM_SPEC_ this%condensed_data_int(2)
#define BASE_RATE_ this%condensed_data_real(1)
#define SCALING_ this%condensed_data_real(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 2
#define NUM_ENV_PARAM_ 1
#define REACT_(s) this%condensed_data_int(NUM_INT_PROP_+s)
#define DERIV_ID_(s) this%condensed_data_int(NUM_INT_PROP_+NUM_SPEC_+s)
#define JAC_ID_(s) this%condensed_data_int(NUM_INT_PROP_+2*NUM_SPEC_+s))

public :: rxn_wet_deposition_t, rxn_update_data_wet_deposition_rate_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_wet_deposition_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Set the reaction id for this reaction
    procedure :: set_rxn_id
    !> Get the reaction property set
    procedure :: get_property_set
    !> Finalize the reaction
    final :: finalize
  end type rxn_wet_deposition_t

  !> Constructor for rxn_wet_deposition_t
  interface rxn_wet_deposition_t
    procedure :: constructor
  end interface rxn_wet_deposition_t

  !> Wet Deposition rate update object
  type, extends(rxn_update_data_t) :: rxn_update_data_wet_deposition_rate_t
  private
    logical :: is_malloced = .false.
  contains
    !> Initialize update data
    procedure :: initialize => update_data_rate_initialize
    !> Update the rate data
    procedure :: set_rate => update_data_rate_set
    !> Finalize the rate update data
    final :: update_data_rate_finalize
  end type rxn_update_data_wet_deposition_rate_t

  !> Interface to c reaction functions
  interface

    !> Allocate space for a rate update
    function rxn_wet_deposition_create_rate_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function rxn_wet_deposition_create_rate_update_data

    !> Set a new wet_deposition rate
    subroutine rxn_wet_deposition_set_rate_update_data(update_data, &
              rxn_id, base_rate) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Reaction id from pmc_rxn_wet_deposition::rxn_wet_deposition_t::set_rxn_id
      integer(kind=c_int), value :: rxn_id
      !> New pre-scaling base wet_deposition rate
      real(kind=c_double), value :: base_rate
    end subroutine rxn_wet_deposition_set_rate_update_data

    !> Free an update rate data object
    pure subroutine rxn_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine rxn_free_update_data

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Wet Deposition reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_wet_deposition_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_wet_deposition_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props
    type(string_t), allocatable :: unique_names(:)
    character(len=:), allocatable :: key_name, phase_name
    integer(kind=i_kind) :: i_rep, i_spec, i_rep_spec, num_spec

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the reaction property set
    call assert_msg(368664748, associated(this%property_set), &
            "Missing property set needed to initialize reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(133499444, &
            this%property_set%get_string(key_name, phase_name), &
            "Wet Deposition reaction is missing aerosol phase name")

    ! Check for aerosol representations
    call assert_msg(674938531, associated(aero_rep), &
            "Missing aerosol representation for wet deposition reaction")
    call assert_msg(731323851, size(aero_rep).gt.0, &
            "Missing aerosol representation for wet deposition reaction")

    ! Count the total number of species in the specified phase in each
    ! aerosol representation
    num_spec = 0
    do i_rep = 1, size(aero_rep)
      unique_names = aero_rep(i_rep)%val%unique_names( phase_name = &
                                                       phase_name )
      num_spec = num_spec + size( unique_names )
    end do
    call assert_msg(332795980, num_spec.gt.0, &
                    "No species found for wet deposition aerosol phase "// &
                    phase_name)

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_+3*num_spec))
    allocate(this%condensed_data_real(NUM_REAL_PROP_))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Save the number of species
    NUM_SPEC_ = num_spec

    ! Get reaction parameters
    key_name = "rate const"
    if (.not. this%property_set%get_real(key_name, BASE_RATE_)) then
      BASE_RATE_ = real(0.0, kind=dp)
    end if
    key_name = "scaling factor"
    if (.not. this%property_set%get_real(key_name, SCALING_)) then
      SCALING_ = real(1.0, kind=dp)
    end if

    ! Save the indices of each species undergoing wet deposition
    i_spec = 0
    do i_rep = 1, size(aero_rep)
      unique_names = aero_rep(i_rep)%val%unique_names( phase_name = &
                                                       phase_name )
      do i_rep_spec = 1, size(unique_names)
        i_spec = i_spec + 1
        REACT_(i_spec) = aero_rep(i_rep)%val%spec_state_id( &
                                         unique_names( i_rep_spec )%string )
        call assert( 702159475, REACT_(i_spec) .gt. 0 )
      end do
    end do
    call assert(312643342, i_spec .eq. NUM_SPEC_)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this reaction that can be used by an external
  !! module to update the base (unscaled) rate constant during the model run.
  subroutine set_rxn_id(this, rxn_id)

    !> Reaction data
    class(rxn_wet_deposition_t), intent(inout) :: this
    !> Reaction id
    integer(kind=i_kind), intent(in) :: rxn_id

    RXN_ID_ = rxn_id

  end subroutine set_rxn_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction properties. (For use by external modules.)
  function get_property_set(this) result(prop_set)

    !> Reaction properties
    type(property_t), pointer :: prop_set
    !> Reaction data
    class(rxn_wet_deposition_t), intent(in) :: this

    prop_set => this%property_set

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_wet_deposition_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for wet_deposition rate constants
  subroutine update_data_rate_set(this, rxn_id, base_rate)

    !> Update data
    class(rxn_update_data_wet_deposition_rate_t), intent(inout) :: this
    !> Reaction id from
    !! \c pmc_rxn_wet_deposition::rxn_wet_deposition_t::set_rxn_id
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Updated pre-scaling wet_deposition rate
    real(kind=dp), intent(in) :: base_rate

    call rxn_wet_deposition_set_rate_update_data(this%get_data(), rxn_id, &
            base_rate)

  end subroutine update_data_rate_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize update data
  subroutine update_data_rate_initialize(this, rxn_type)

    !> Update data object
    class(rxn_update_data_wet_deposition_rate_t) :: this
    !> Reaction type id
    integer(kind=i_kind), intent(in) :: rxn_type

    this%rxn_type = int(rxn_type, kind=c_int)
    this%update_data = rxn_wet_deposition_create_rate_update_data()
    this%is_malloced = .true.

  end subroutine update_data_rate_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize an update data object
  elemental subroutine update_data_rate_finalize(this)

    !> Update data object to free
    type(rxn_update_data_wet_deposition_rate_t), intent(inout) :: this

    if (this%is_malloced) call rxn_free_update_data(this%update_data)

  end subroutine update_data_rate_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef RXN_ID_
#undef NUM_SPEC_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef REACT_
#undef DERIV_ID_
#undef JAC_ID_
end module pmc_rxn_wet_deposition
