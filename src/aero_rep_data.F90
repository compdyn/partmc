! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_data module.

!> \page phlex_aero_rep Phlexible Module for Chemistry: Aerosol Representation (general)
!!
!! An aerosol representation acts as an interface between the aerosol
!! micro-physics module and the chemistry module. A unique derived type that
!! extends the abstract \c aero_rep_data_t type should be developed for each
!! type of aerosol representation used by an external model (e.g., binned,
!! modal, single particle). This derived type can either calculate the aerosol
!! microphysics directly, or use the microphysics results of the host model to
!! update the \ref phlex_aero_phase "aerosol phase" states in the chemistry
!! module.
!!
!! The \c aero_rep_data_t type includes three deferred functions for
!! calculating physical aerosol conditions &mdash; surface area
!! concentration, species surface area concentration, and Kelvin Effect
!! vapor pressure scaling. Each type extending \c aero_rep_data_t should
!! provide these three functions. As more physical aerosol parameter functions
!! are required by specific reaction types, overridable type-bound procedures
!! should be added to the \c aero_rep_data_t type, whose default action is to
!! throw an error, so that only extending types that override the new
!! procedure can be used by reactions that require it.
!!
!! The available aerosol representations are:
!!  - \ref phlex_aero_rep_single_particle "Single Particle"
!!
!! The general input format for an aerosol representation can be found \ref
!! input_format_aero_rep "here".
!!
!! General instructions for adding a new aerosol representation can be found
!! \ref phlex_aero_rep_add "here".

!> The abstract aero_rep_data_t structure and associated subroutines.
module pmc_aero_rep_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
  use pmc_chem_spec_data
  use pmc_phlex_state
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
  !!
  !! See \ref phlex_aero_rep "Aerosol Representations" for details.
  type, abstract :: aero_rep_data_t
    private
    !> Name of the aerosol representation
    character(len=:), allocatable, public :: rep_name
    !> Aerosol phases associated with this aerosol scheme
    !!
    !! See \ref phlex_aero_phase "Aerosol Phases" for details.
    type(aero_phase_data_ptr), allocatable, public :: aero_phase(:)
    !> Aerosol representation parameters. These will be available during 
    !! initialization, but not during integration. All information required
    !! by functions of the aerosol representation  must be saved by the
    !! exdending type in the condensed data arrays.
    type(property_t), pointer, public :: property_set => null()
    !> Condensed representaiton data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the pmc_phlex_state::phlex_state_t object. (floating-point)
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed representaiton data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the pmc_phlex_state::phlex_state_t object. (integer)
    integer(kind=i_kind), allocatable, public ::  condensed_data_int(:)
  contains
    !> Initialize the aerosol representation data, validating component data and
    !! loading any required information from the \c
    !! aero_rep_data_t::property_set. This routine should be called once for
    !! each aerosol representation at the beginning of a model run after all
    !! the input files have been read in. It ensures all data required during
    !! the model run are included in the condensed data arrays.
    procedure(initialize), deferred :: initialize
    !> Get the size of the section of the
    !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
    !! aerosol representation
    procedure(get_size), deferred :: size
    !> Get a list of unique names for each element on the
    !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
    !! representation.
    procedure(unique_names), deferred :: unique_names
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
    procedure(state_id_by_unique_name), deferred :: &
            state_id_by_unique_name
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
    procedure(species_state_id), deferred :: &
            species_state_id
    !> Get the absolute integration tolerance for each variable on the
    !! pmc_phlex_state::state_var array from this aerosol representation
    procedure(get_abs_tol), deferred :: get_abs_tol
    !> Get an instance of the pmc_aero_rep_state::aero_rep_state_t variable
    !! for this aerosol representation.
    !!
    !! See phlex_aero_rep_state "Aerosol Representation States" for details.
    procedure(new_state), deferred :: new_state
    !> Get surface area concentration (m^2/m^3) between two phases. One phase
    !! may be set to 0 to indicate the gas phase, or both phases may be
    !! aerosol phases.
    procedure(surface_area_conc), deferred :: &
            surface_area_conc
    !> Get the surface area concentration for a specific aerosol species
    !! (m^2/m^3) between a specified aerosol phase and the gas phase.
    procedure(species_surface_area_conc), deferred :: &
            species_surface_area_conc
    !> Get the Kelvin effect vapor pressure adjustment for a particular 
    !! species (unitless)
    procedure(kelvin_effect), deferred :: &
            kelvin_effect
    !> Get the name of the aerosol representation
    procedure :: name => get_name
    !> Get a phase id in this aerosol representation. This id should be used
    !! for \c aero_rep_data_t functions requiring a phase id. Note that a
    !! particiular aerosol representation may include multiple instances of
    !! the same phase. (e.g., a binned aerosol representation with 8 bins may
    !! include 8 "aqueous" aerosol phases). It is assumed that the
    !! intra-phase chemistry and inter-phase mass transfer depend only on the
    !! type of phase(s) involed, and not the particular instance of a phase
    !! or interface. Thus, each instance \f$i\f$ of phase \f$p\f$ would be
    !! related to the same phase id \f$x_p\f$.
    !!
    !! Returns 0 if the phase is not present in the aerosol representation.
    !!
    !! This function should only be called during initialization
    procedure :: phase_id
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Load data from an input file
    procedure :: load
    !> Print the aerosol representation data
    procedure :: print => do_print
  end type aero_rep_data_t

  !> Pointer to aero_rep_data_t extending types
  type :: aero_rep_data_ptr
    !> Pointer to an aerosol representation
    class(aero_rep_data_t), pointer :: val
  end type aero_rep_data_ptr

interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the \c
  !! aero_rep_data_t::property_set. This routine should be called once for
  !! each aerosol representation at the beginning of a model run after all
  !! the input files have been read in. It ensures all data required during
  !! the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_phase_set, &
                spec_state_id, aero_state_id, chem_spec_data)
    use pmc_util,                                     only : i_kind
    use pmc_chem_spec_data
    use pmc_aero_phase_data
    import :: aero_rep_data_t

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this
    !> The set of aerosol phases. Note that an aerosol representation may
    !! implement any number of instances of each phase.
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representation in the
    !! \c pmc_phlex_state::phlex_state_t::state_var array
    integer(kind=i_kind), intent(in) :: spec_state_id
    !> Index for this representation in then 
    !! \c pmc_phlex_state::phlex_state_t::aero_rep_state array
    integer(kind=i_kind), intent(in) :: aero_state_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the section of the
  !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
  !! aerosol representation
  function get_size(this) result (state_size)
    use pmc_util,                                     only : i_kind
    import :: aero_rep_data_t

    !> Size of the state array section
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element on the
  !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
  !! representation.
  function unique_names(this) 
    use pmc_util,                                     only : string_t
    import :: aero_rep_data_t

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

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
    use pmc_util,                                     only : i_kind
    import :: aero_rep_data_t

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Unique name
    character(len=:), allocatable, intent(in) :: unique_name

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
  function species_state_id(this, phase_name, &
                    species_name) result(spec_index)
    use pmc_util,                                     only : i_kind
    import :: aero_rep_data_t

    !> Species index array
    integer(kind=i_kind), allocatable :: spec_index(:)
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase
    character(len=:), allocatable, intent(in) :: phase_name
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

  end function species_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the absolute integration tolerance for each variable on the
  !! pmc_phlex_state::state_var array from this aerosol representation
  subroutine get_abs_tol(this, chem_spec_data, abs_tol)
    use pmc_util,                                       only : dp
    use pmc_chem_spec_data
    import :: aero_rep_data_t

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Absolute integration tolerance
    real(kind=dp), allocatable :: abs_tol(:)

  end subroutine get_abs_tol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an instance of the pmc_aero_rep_state::aero_rep_state_t variable
  !! for this aerosol representation.
  !!
  !! See phlex_aero_rep_state "Aerosol Representation States" for details.
  function new_state(this) result (aero_rep_state)
    use pmc_aero_rep_state
    import :: aero_rep_data_t

    !> Aerosol representation state
    class(aero_rep_state_t), pointer :: aero_rep_state
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

  end function new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get surface area concentration (m^2/m^3) between two phases. One phase
  !! may be set to 0 to indicate the gas phase, or both phases may be
  !! aerosol phases.
  function surface_area_conc(this, i_phase1, i_phase2, &
                    phlex_state, jac_contrib)
    use pmc_util,                                     only : dp, i_kind
    use pmc_phlex_state
    import :: aero_rep_data_t

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> First aerosol phase index. Use \c aero_rep_data_t%phase_id() to get this
    !! index
    integer(kind=i_kind), intent(in) :: i_phase1
    !> Second aerosol phase index. Use \c aero_rep_data_t%phase_id() to get this
    !! index
    integer(kind=i_kind), intent(in) :: i_phase2
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! \c pmc_phlex_state::phlex_state_t::state_var array that, when present,
    !! will be filled with the partial derivatives of the result of this
    !! calculation by each state variable.
    real(kind=dp), allocatable, intent(inout), optional :: jac_contrib(:)

  end function surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the surface area concentration for a specific aerosol species
  !! (m^2/m^3) between a specified aerosol phase and the gas phase.
  function species_surface_area_conc(this, i_phase, &
                    i_spec, phlex_state, jac_contrib) result(surface_area_conc)
    use pmc_util,                                     only : dp, i_kind
    use pmc_phlex_state
    import :: aero_rep_data_t

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase id. Use the \c phase_id() function of \c aero_rep_data_t
    !! to get this index.
    integer(kind=i_kind), intent(in) :: i_phase
    !> Species id
    integer(kind=i_kind), intent(in) :: i_spec
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! \c pmc_phlex_state::phlex_state_t::state_var array that, when present,
    !! will be filled with the partial derivatives of the result of this
    !! calculation by each state variable.
    real(kind=dp), allocatable, intent(inout), optional :: jac_contrib(:)

  end function species_surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
  !> Get the Kelvin effect vapor pressure adjustment for a particular species 
  !! (unitless)
  function kelvin_effect(this, i_spec, phlex_state, jac_contrib)
    use pmc_util,                                     only : dp, i_kind
    use pmc_phlex_state
    import :: aero_rep_data_t

    !> Vapor pressure scaling
    real(kind=dp) :: kelvin_effect
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Species index. This index is the index of species \f$s\f$ in phase
    !! \f$p\f$ and should be found using the \c spec_id() function of
    !! \c aero_phase_t.
    integer(kind=i_kind), intent(in) :: i_spec
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! \c pmc_phlex_state::phlex_state_t::state_var array that, when present,
    !! will be filled with the partial derivatives of the result of this
    !! calculation by each state variable.
    real(kind=dp), allocatable, intent(inout), optional :: jac_contrib(:)

  end function kelvin_effect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the aerosol representation
  function get_name(this)

    !> Aerosol representation name
    character(len=:), allocatable :: get_name
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    get_name = this%rep_name

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a phase id in this aerosol representation. This id should be used
  !! for \c aero_rep_data_t functions requiring a phase id. Note that a
  !! particiular aerosol representation may include multiple instances of
  !! the same phase. (e.g., a binned aerosol representation with 8 bins may
  !! include 8 "aqueous" aerosol phases). It is assumed that the
  !! intra-phase chemistry and inter-phase mass transfer depend only on the
  !! type of phase(s) involed, and not the particular instance of a phase
  !! or interface. Thus, each instance \f$i\f$ of phase \f$p\f$ would be
  !! related to the same phase id \f$x_p\f$.
  !!
  !! Returns 0 if the phase is not present in the aerosol representation.
  !!
  !! This function should only be called during initialization
  integer(kind=i_kind) function phase_id(this, phase_name)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase to find
    character(len=:), allocatable, intent(in) :: phase_name

    integer(kind=i_kind) :: i_phase

    phase_id = 0
    do i_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        phase_id = i_phase
        return
      end if
    end do

  end function phase_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol 
  !! representation data
  integer(kind=i_kind) function pack_size(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int) 

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

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

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

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

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_aero_rep Input JSON Object Format: Aerosol Representation (general)
  !!
  !! A \c json object containing information about an \ref phlex_aero_rep
  !! "aerosol representation" of the form:
  !! \code{.json}
  !! { "pmc-data" : [
  !!   {
  !!     "type" : "AERO_REP_TYPE",
  !!     "some parameter" : 123.34,
  !!     "some other parameter" : true,
  !!     "nested parameters" : {
  !!        "sub param 1" : 12.43,
  !!        "sub param other" : "some text"
  !!   },
  !!   ...
  !! ]}
  !! \endcode
  !! Aerosol representations must have a unique \b type that corresponds to a
  !! valid aerosol representation type. These include:
  !!
  !!   - \subpage phlex_aero_rep_single_particle "Single Particle"
  !!
  !! All remaining data are optional and may include any valid \c json value, 
  !! including nested objects. However, extending types will have specific
  !! requirements for the remaining data. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol representation from an input file
#ifdef PMC_USE_JSON
  subroutine load(this, json, j_obj)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this
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
  subroutine load(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this

    call warn_msg(433045149, "No support for input files")
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the aerosol representation data
  subroutine do_print(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    if (associated(this%property_set)) call this%property_set%print()

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_data
