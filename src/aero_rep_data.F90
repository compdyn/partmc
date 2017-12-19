! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_data module.

!> The abstract aero_rep_data_t structure and associated subroutines.
module pmc_aero_rep_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
  use pmc_chem_spec_data
  use pmc_model_state
  use pmc_aero_phase_data
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: aero_rep_data_t, aero_rep_data_ptr

  !> Abstract aerosol representation data type
  !!
  !! Time-invariant data related to an aerosol representation. Derived types
  !! extending aero_rep_data_t should describe specific types of aerosol
  !! schemes (e.g., binned, modal, particle-resolved).
  type, abstract :: aero_rep_data_t
    private
    !> Name of the aerosol representation
    character(len=:), allocatable :: rep_name
    !> Aerosol phases associated with this aerosol scheme
    class(aero_phase_data_ptr), allocatable :: aero_phase
    !> Aerosol representation parameters. These will be available during 
    !! initialization, but not during integration. All information required
    !! by functions of the aerosol representation  must be saved by the
    !! exdending type in the condensed data arrays.
    type(property_t), pointer :: property_set => null()
    !> Condensed representaiton data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the model_state_t object. (floating-point)
    real(kind=dp), allocatable :: condensed_data_real(:)
    !> Condensed reaction data (integers)
    integer(kind=i_kind), allocatable ::  condensed_data_int(:)
  contains
    !> Aerosol representation initialization
    procedure(pmc_aero_rep_data_initialize), deferred :: initialize
    !> Get aerosol species state id
    procedure(pmc_aero_rep_data_species_state_id), deferred :: &
            species_state_id
    !> Get surface area concentration (m^2/m^3)
    procedure(pmc_aero_rep_data_surface_area_conc), deferred :: &
            surface_area_conc
    !> Get the surface area concentration for a specific aerosol species
    !! (m^2/m^3)
    procedure(pmc_aero_rep_data_species_surface_area_conc), deferred :: &
            species_surface_area_conc
    !> Get the vapor pressure scaling for a particular species (unitless)
    procedure(pmc_aero_rep_data_vapor_pressure_scaling), deferred :: &
            vapor_pressure_scaling
    !> Get the name of the aerosol representation
    procedure :: name => pmc_aero_rep_data_name
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size => pmc_aero_rep_data_pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack => pmc_aero_rep_data_bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_aero_rep_data_bin_unpack
    !> Load data from an input file
    procedure :: load => pmc_aero_rep_data_load
    !> Print the aerosol representation data
    procedure :: print => pmc_aero_rep_data_print
  end type aero_rep_data_t

  !> Pointer to aero_rep_data_t extending types
  type :: aero_rep_data_ptr
    !> Pointer to an aerosol representation
    class(aero_rep_data_t), pointer :: val
  end type aero_rep_data_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the  property_t object. This 
  !! routine should be called once for each aerosol representation
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  interface pmc_aero_rep_data_initialize_if
    subroutine pmc_aero_rep_data_initialize(this, aero_phase_set)
      import :: aero_rep_data_t, aero_phase_data_t

      !> Aerosol representation data
      class(aero_rep_data_t), intent(inout) :: this
      !> Aerosol phase data
      type(aero_phase_data_t), pointer, intent(in) :: aero_phase_set(:)

    end subroutine pmc_aero_rep_data_initialize
  end interface pmc_aero_rep_data_initialize_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an aerosol species state id(s). The returned array will have the
  !! index of the species in the specified phase in each aerosol group (e.g.,
  !! bin, mode). If the species is not present in a certain group-phase the
  !! index will be 0.
  interface pmc_aero_rep_data_species_state_id_if
    function pmc_aero_rep_data_species_state_id(this, phase, species_name) &
                    result(spec_index)
      use pmc_util,                                     only : i_kind
      import :: aero_rep_data_t

      !> Species index array
      integer(kind=i_kind), allocatable :: spec_index(:)
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase
      character(len=:), allocatable, intent(in) :: phase
      !> Species name
      character(len=:), allocatable, intent(in) :: species_name

    end function pmc_aero_rep_data_species_state_id
  end interface pmc_aero_rep_data_species_state_id_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get surface area concentration (m^2/m^3) between two phases. One phase
  !! may be 'gas' to indicate the gas-phase, or both phases may be aerosol
  !! phases.
  interface pmc_aero_rep_data_surface_area_conc_if 
    function pmc_aero_rep_data_surface_area_conc(this, phase1, phase2) &
                    result(surface_area_conc)
      use pmc_util,                                     only : dp
      import :: aero_rep_data_t

      !> Surface area concentration
      real(kind=dp) :: surface_area_conc
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase1
      character(len=:), allocatable, intent(in) :: phase1
      !> Aerosol phase2
      character(len=:), allocatable, intent(in) :: phase2

    end function pmc_aero_rep_data_surface_area_conc
  end interface pmc_aero_rep_data_surface_area_conc_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the surface area concentration for a specific aerosol species
  !! (m^2/m^3). It is assumed the surface is between the gas-phase and an
  !! aerosol phase.
  interface pmc_aero_rep_data_species_surface_area_conc_if
    function pmc_aero_rep_data_species_surface_area_conc(this, phase, species_name) &
                    result(surface_area_conc)
      use pmc_util,                                     only : dp
      import :: aero_rep_data_t

      !> Surface area concentration
      real(kind=dp) :: surface_area_conc
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase
      character(len=:), allocatable, intent(in) :: phase
      !> Species name
      character(len=:), allocatable, intent(in) :: species_name

    end function pmc_aero_rep_data_species_surface_area_conc
  end interface pmc_aero_rep_data_species_surface_area_conc_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
  !> Get the vapor pressure scaling for a particular species (unitless)
  interface pmc_aero_rep_data_vapor_pressure_scaling_if  
    function pmc_aero_rep_data_vapor_pressure_scaling(this, species_name) &
                    result(vapor_pressure_scaling)
      use pmc_util,                                     only : dp
      import :: aero_rep_data_t

      !> Vapor pressure scaling
      real(kind=dp) :: vapor_pressure_scaling
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Species name
      character(len=:), allocatable, intent(in) :: species_name

    end function pmc_aero_rep_data_vapor_pressure_scaling
  end interface pmc_aero_rep_data_vapor_pressure_scaling_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the aerosol representation
  function pmc_aero_rep_data_name(this) result(name)

    !> Aerosol representation name
    character(len=:), allocatable :: name
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    name = this%rep_name

  end function pmc_aero_rep_data_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol 
  !! representation data
  integer(kind=i_kind) function pmc_aero_rep_data_pack_size(this) &
                  result (pack_size)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int) 

  end function pmc_aero_rep_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_aero_rep_data_bin_pack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_pack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(257024095, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_aero_rep_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine pmc_aero_rep_data_bin_unpack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(954732699, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_aero_rep_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol representation from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to an
  !! aerosol representation. It should be of the form:
  !! 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my representation",
  !!     "type" : "AERO_REP_TYPE"
  !!     "some parameter" : 123.34,
  !!     "some other parameter" : true,
  !!     "nested parameters" : {
  !!        "sub param 1" : 12.43,
  !!        "sub param other" : "some text"
  !!   },
  !!   ...
  !! ]}
  !!
  !! Aerosol representations must have a unique name and have a "type" that
  !! corresponds to a valid aerosol representation type.
  !!
  !! All remaining data are optional and may include any valid JSON value, 
  !! including nested objects. However, extending types will have specific
  !! requirements for the remaining data. 
  subroutine pmc_aero_rep_data_load(this, json, j_obj)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(out) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next, species
    character(kind=json_ck, len=:), allocatable :: key

    this%property_set => property_t()

    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key)
      if (key.ne."type" .and. key.ne."name") &
              call this%property_set%load(json, child, .false.)
      call json%get_next(child, next)
      child => next
    end do
#else
  subroutine pmc_aero_rep_data_load(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this

    call warn_msg(433045149, "No support for input files")
#endif

  end subroutine pmc_aero_rep_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the aerosol representation data
  subroutine pmc_aero_rep_data_print(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    if (associated(this%property_set)) call this%property_set%print()

  end subroutine pmc_aero_rep_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_data
