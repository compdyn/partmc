! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_factory module.

!> The abstract aero_rep_factory_t structure and associated subroutines.
module pmc_aero_rep_factory

  use pmc_constants,                  only : i_kind, dp
  use pmc_util,                       only : die_msg, string_t, assert_msg
  use pmc_aero_rep_data
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: aero_rep_factory_t

  !> Factory type for aerosol representations
  !!
  !! Provides new instances of types extending aero_rep_data_t by name or
  !! from input file data
  type :: aero_rep_factory_t
  contains
    !> Create a new aerosol representation by type name
    procedure :: create => pmc_aero_rep_factory_create
    !> Create a new aerosol representation from input data
    procedure :: load => pmc_aero_rep_factory_load
  end type aero_rep_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new aerosol representation by type name
  function pmc_aero_rep_factory_create(this, type_name) result (new_obj)

    use pmc_aero_rep_single_particle

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Name of the aerosol representation
    character(len=:), allocatable :: type_name

    new_obj => null()

    select case (type_name)
      case ("AERO_REP_SINGLE_PARTICLE")
        new_obj => aero_rep_single_particle_t()
      case default
        call die_msg(792930166, "Unknown aerosol representation type: " &
                //type_name) 
    end select

  end function pmc_aero_rep_factory_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol represenation based on its type
#ifdef PMC_USE_JSON
  function pmc_aero_rep_factory_load(this, json, j_obj) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child
    character(kind=json_ck, len=:), allocatable :: key, unicode_type_name
    character(len=:), allocatable :: type_name
    logical(kind=json_lk) :: found

    new_obj => null()

    call json%get(j_obj, "type", unicode_type_name, found)
    call assert_msg(904082395, found, 'Missing aerosol representation type.')
    type_name = unicode_type_name
    new_obj => this%create(type_name)    
    call new_obj%load(json, j_obj)
#else
  function pmc_aero_rep_factory_load(this) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(723960750, "No support for input files.")
#endif
  end function pmc_aero_rep_factory_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_factory
