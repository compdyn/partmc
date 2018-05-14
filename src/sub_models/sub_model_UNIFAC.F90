! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_sub_model_UNIFAC module.

!> \page phlex_sub_model_UNIFAC Phlexible Module for Chemistry: UNIFAC Activity Coefficients
!!
!! The UNIFAC activity coefficient sub model calculates activity coefficients
!! for species in an aerosol phase based on the current aerosol phase
!! composition. The \c json object for this \ref phlex_sub_model "sub model"
!! is of the form:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "type" : "SUB_MODEL_UNIFAC",
!!      TODO finish
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be \b SUB_MODEL_UNIFAC.

!> The sub_model_UNIFAC_t type and assocatiated subroutines
module pmc_sub_model_UNIFAC

  use pmc_util,                                 only : dp, i_kind, &
                                                       string_t, assert_msg, &
                                                       die_msg, to_string
  use pmc_property
  use pmc_sub_model_data
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_phlex_state

  implicit none
  private

#define _NUM_SPEC_ this%condensed_data_int(1)
#define _NUM_INT_PROP_ 1
#define _NUM_REAL_PROP_ 0

  ! Update types (These must match values in sub_model_UNIFAC.c)
  ! (none for now)

  public :: sub_model_UNIFAC_t

  !> UNIFAC activity coefficient calculation
  !!
  !! Time-invariant data required by the UNIFAC activity coefficient sub model
  type, extends(sub_model_data_t) :: sub_model_UNIFAC_t
  contains
    !> Initialize the sub model data, validating input parameters and
    !! loading any required information form the \c
    !! sub_model_data_t::property_set. This routine should be called
    !! once for each sub model at the beginning of the model run after all
    !! the input files have been read in. It ensures all data required
    !! during the model run are included in the condensed data arrays.
    procedure :: initialize
  end type sub_model_UNIFAC_t

  ! Constructor for sub_model_UNIFAC_t
  interface sub_model_UNIFAC_t
    procedure :: constructor
  end interface sub_model_UNIFAC_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for sub_model_UNIFAC_t
  function constructor() result (new_obj)

    !> New sub model
    type(sub_model_UNIFAC_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !> Initialize the sub model data, validating input parameters and
  !! loading any required information form the \c
  !! sub_model_data_t::property_set. This routine should be called
  !! once for each sub model at the beginning of the model run after all
  !! the input files have been read in. It ensures all data required
  !! during the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_rep_set, chem_spec_data)

    !> Sub model data
    class(sub_model_UNIFAC_t), intent(inout) :: this
    !> The set of aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep_set(:)
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    ! TODO finish

  end subroutine initialize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_SPEC_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
end module pmc_sub_model_UNIFAC
