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
!!      "phases" : [
!!        "some phase",
!!        "some other phase"
!!      ],
!!      "functional groups" : {
!!        "CH3" : {
!!          "main group" : "CH2",
!!          "volume param" : 0.9011,
!!          "surface param" : 0.8480
!!        },
!!        "CH2" : {
!!          "main group" : "CH2",
!!          "volume param" : 0.6744,
!!          "surface param" : 0.5400
!!        },
!!        "CH=CH" : {
!!          "main group" : "C=C",
!!          "volume param" : 1.1167,
!!          "suface param" : 0.8670
!!        }
!!      },
!!      "main groups" : {
!!        "CH2" : {
!!          "interactions with" : {
!!            "C=C" : -35.36
!!          }
!!        },
!!        "C=C" : {
!!          "interactions with" : {
!!            "CH2" : 86.02
!!          }
!!        }
!!      }
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be \b SUB_MODEL_UNIFAC.
!! The key-value pair \b phases is also required, and its value must be an
!! array of strings that correspond to valid 
!! \ref phlex_aero_phase "aerosol phases". The key-value pair \b "functional
!! groups" is also required, and must contain a set of key-value pairs whose
!! keys are the names of UNIFAC functions groups, and whose values are a set
!! of key value pairs that contain, at minimum:
!!   - \b "main group" : a string that corresponds to a key in the \b
!!                       "main groups" set.
!!   - \b "volume param" : the floating-point volume parameter for this 
!!                         functional group.
!!   - \b "surface param" : this floating-point surface parameter for this
!!                          functional group.
!! The last required key-value pair is \b "main groups" whose value must
!! be a set of key-value pairs whose keys are the names of the UNIFAC main
!! groups and whose values are a set key-pairs that contain, at minimum,
!! \b "interaction with" whose value is a set of key-value pairs whose keys
!! are the names of the other \b "main groups" and whose values are the
!! floating-point interation parameters for that interaction. Each main group
!! may contain up to one interaction with each other main group, and may
!! not contain an interaction with itself. Missing interactions are assumed
!! to be 0.0.
!!
!! Species in the specified phase for whom acitivity coefficients will be
!! calculated must contain a key-value pair \b "UNIFAC groups" whose value
!! is a set of key value pairs that correspond with members of the
!! \b "functional groups" set and whose values are the integer number of
!! instances of a particular functional group in this species. For the
!! above example UNIFAC model, the following species would be valid and
!! included in activity coefficient calculations:
!! \code{.json}
!! { "pmc-data" : [
!!   {
!!     "name" : "my species",
!!     "type" : "CHEM_SPEC",
!!     "phase" : "AEROSOL",
!!     "UNIFAC groups" : {
!!       "CH3" : 4,
!!       "C=C" : 1
!!     }
!!   },
!!   {
!!     "name" : "my other species",
!!     "type" : "CHEM_SPEC",
!!     "phase" : "AEROSOL",
!!     "UNIFAC groups" : {
!!       "CH3" : 2,
!!       "CH2" : 4
!!     },
!!   },
!!   {
!!     "name" : "some phase",
!!     "type" : "AERO_PHASE",
!!     "species" : { "my species", "my other species" }
!!   }
!! ]}
!! \endcode

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
