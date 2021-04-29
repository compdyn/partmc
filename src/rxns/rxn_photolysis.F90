! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_photolysis module.

!> \page camp_rxn_photolysis CAMP: Photolysis
!!
!! Photolysis reactions take the form:
!!
!! \f[\ce{
!!   X + h $\nu$ -> Y_1 ( + Y_2 \dots )
!! }\f]
!!
!! where \f$\ce{X}\f$ is the species being photolyzed, and
!! \f$\ce{Y_n}\f$ are the photolysis products.
!!
!! Photolysis rate constants (including the \f$\ce{h $\nu$}\f$ term)
!! can be constant or set from an external photolysis module using the
!! \c pmc_rxn_photolysis::rxn_update_data_photolysis_t object.
!! External modules can use the
!! \c pmc_rxn_photolysis::rxn_photolysis_t::get_property_set() function during
!! initilialization to access any needed reaction parameters to identify
!! certain photolysis reactions.
!! An \c pmc_rxn_photolysis::update_data_photolysis_t object should be
!! initialized for each photolysis reaction. These objects can then be used
!! during solving to update the photolysis rate from an external module.
!!
!! Input data for photolysis reactions have the following format :
!! \code{.json}
!!   {
!!     "type" : "PHOTOLYSIS",
!!     "reactants" : {
!!       "spec1" : {}
!!     },
!!     "products" : {
!!       "spec2" : {},
!!       "spec3" : { "yield" : 0.65 },
!!       ...
!!     },
!!     "scaling factor" : 1.2,
!!     ...
!!   }
!! \endcode
!! The key-value pairs \b reactants, and \b products are required. There must
!! be exactly one key-value pair in the \b reactants object whose name is the
!! species being photolyzed and whose value is an empty \c json object. Any
!! number of products may be present. Products without a specified \b yield
!! are assumed to have a \b yield of 1.0. The \b scaling \b factor is
!! optional, and can be used to set a constant scaling factor for the rate
!! constant. When the \b scaling \b factor is not provided, it is assumed to
!! be 1.0. All other data is optional and will be available to external
!! photolysis modules during initialization. Rate constants are in units of
!! \f$s^{-1}\f$, and must be set using a \c rxn_photolysis_update_data_t
!! object.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_photolysis_t type and associated functions.
module pmc_rxn_photolysis

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_camp_state
  use pmc_mpi
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg

  use iso_c_binding

  implicit none
  private

#define NUM_REACT_ this%condensed_data_int(1)
#define NUM_PROD_ this%condensed_data_int(2)
#define RXN_ID_ this%condensed_data_int(3)
#define SCALING_ this%condensed_data_real(1)
#define NUM_INT_PROP_ 3
#define NUM_REAL_PROP_ 1
#define NUM_ENV_PARAM_ 2
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_ + x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_ + x)

public :: rxn_photolysis_t, rxn_update_data_photolysis_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_photolysis_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Get the reaction property set
    procedure :: get_property_set
    !> Initialize update data
    procedure :: update_data_initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_photolysis_t

  !> Constructor for rxn_photolysis_t
  interface rxn_photolysis_t
    procedure :: constructor
  end interface rxn_photolysis_t

  !> Photolysis rate update object
  type, extends(rxn_update_data_t) :: rxn_update_data_photolysis_t
  private
    !> Flag indicating whether the update data as been allocated
    logical :: is_malloced = .false.
    !> Unique id for finding reactions during model initialization
    integer(kind=i_kind) :: rxn_unique_id = 0
  contains
    !> Update the rate data
    procedure :: set_rate => update_data_rate_set
    !> Determine the pack size of the local update data
    procedure :: internal_pack_size
    !> Pack the local update data to a binary
    procedure :: internal_bin_pack
    !> Unpack the local update data from a binary
    procedure :: internal_bin_unpack
    !> Finalize the rate update data
    final :: update_data_finalize
  end type rxn_update_data_photolysis_t

  !> Interface to c reaction functions
  interface

    !> Allocate space for a rate update
    function rxn_photolysis_create_rate_update_data() result (update_data) &
              bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function rxn_photolysis_create_rate_update_data

    !> Set a new photolysis rate
    subroutine rxn_photolysis_set_rate_update_data(update_data, photo_id, &
              base_rate) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Photo id
      integer(kind=c_int), value :: photo_id
      !> New pre-scaling base photolysis rate
      real(kind=c_double), value :: base_rate
    end subroutine rxn_photolysis_set_rate_update_data

    !> Free an update rate data object
    pure subroutine rxn_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine rxn_free_update_data

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Photolysis reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_photolysis_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_photolysis_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(408416753, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(415586594, &
            this%property_set%get_property_t(key_name, reactants), &
            "Photolysis reaction is missing reactants")
    key_name = "products"
    call assert_msg(180421290, &
            this%property_set%get_property_t(key_name, products), &
            "Photolysis reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(240165383, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the
    ! reactants and products, one molecular property for each reactant,
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(NUM_INT_PROP_ + &
            (i_spec + 2) * (i_spec + products%size())))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + products%size()))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Save the size of the reactant and product arrays (for reactions where
    ! these can vary)
    NUM_REACT_ = i_spec
    NUM_PROD_ = products%size()

    ! Get reaction parameters (it might be easiest to keep these at the
    ! beginning of the condensed data array, so they can be accessed using
    ! compliler flags)
    key_name = "scaling factor"
    if (.not. this%property_set%get_real(key_name, SCALING_)) then
      SCALING_ = real(1.0, kind=dp)
    end if

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(747277322, REACT_(i_spec).gt.0, &
              "Missing photolysis reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(231542303, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_qty = 1, temp_int - 1
          REACT_(i_spec + i_qty) = REACT_(i_spec)
        end do
        i_spec = i_spec + temp_int - 1
      end if

      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Make sure exactly one reactant is present
    call assert_msg(908486656, i_spec.eq.2, "Incorrect number of reactants"//&
            " for Photolysis reaction: "//to_string(i_spec-1))

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(360988742, PROD_(i_spec).gt.0, &
              "Missing photolysis product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(173703744, products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        YIELD_(i_spec) = temp_real
      else
        YIELD_(i_spec) = 1.0
      end if

      call products%iter_next()
      i_spec = i_spec + 1
    end do

    ! Initialize the reaction id
    RXN_ID_ = -1

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction properties. (For use by external photolysis modules.)
  function get_property_set(this) result(prop_set)

    !> Reaction properties
    type(property_t), pointer :: prop_set
    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this

    prop_set => this%property_set

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_photolysis_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for photolysis rate constants
  subroutine update_data_rate_set(this, base_rate)

    !> Update data
    class(rxn_update_data_photolysis_t), intent(inout) :: this
    !> Updated pre-scaling photolysis rate
    real(kind=dp), intent(in) :: base_rate

    call rxn_photolysis_set_rate_update_data(this%get_data(), &
            this%rxn_unique_id, base_rate)

  end subroutine update_data_rate_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize update data
  subroutine update_data_initialize(this, update_data, rxn_type)

    use pmc_rand,                                only : generate_int_id

    !> The reaction to update
    class(rxn_photolysis_t), intent(inout) :: this
    !> Update data object
    class(rxn_update_data_photolysis_t), intent(out) :: update_data
    !> Reaction type id
    integer(kind=i_kind), intent(in) :: rxn_type

    ! If a reaction id has not yet been generated, do it now
    if (RXN_ID_.eq.-1) then
      RXN_ID_ = generate_int_id()
    endif

    update_data%rxn_unique_id = RXN_ID_
    update_data%rxn_type = int(rxn_type, kind=c_int)
    update_data%update_data = rxn_photolysis_create_rate_update_data()
    update_data%is_malloced = .true.

  end subroutine update_data_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the reaction data
  integer(kind=i_kind) function internal_pack_size(this, comm) &
      result(pack_size)

    !> Reaction update data
    class(rxn_update_data_photolysis_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size = &
      pmc_mpi_pack_size_logical(this%is_malloced, comm) + &
      pmc_mpi_pack_size_integer(this%rxn_unique_id, comm)

  end function internal_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine internal_bin_pack(this, buffer, pos, comm)

    !> Reaction update data
    class(rxn_update_data_photolysis_t), intent(in) :: this
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
    call pmc_mpi_pack_integer(buffer, pos, this%rxn_unique_id, comm)
    call assert(649543400, &
         pos - prev_position <= this%pack_size(comm))
#endif

  end subroutine internal_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine internal_bin_unpack(this, buffer, pos, comm)

    !> Reaction update data
    class(rxn_update_data_photolysis_t), intent(inout) :: this
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
    call pmc_mpi_unpack_integer(buffer, pos, this%rxn_unique_id, comm)
    call assert(254749806, &
         pos - prev_position <= this%pack_size(comm))
    this%update_data = rxn_photolysis_create_rate_update_data()
#endif

  end subroutine internal_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize an update data object
  elemental subroutine update_data_finalize(this)

    !> Update data object to free
    type(rxn_update_data_photolysis_t), intent(inout) :: this

    if (this%is_malloced) call rxn_free_update_data(this%update_data)

  end subroutine update_data_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_photolysis
