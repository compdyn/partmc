! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_phase_transfer module.

!> \page phlex_rxn_phase_transfer Phlexible Mechanism for Chemistry: Phase-Transfer Reaction
!!
!! Phase transfer reactions are based on Henry's Law equilibrium constants
!! whose equations take the form:
!!
!! \f[
!!   Ae^{C({1/T-1/298})}
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor (\f$s^{-1}\f$), \f$C\f$ is a
!! constant and \f$T\f$ is the temperature (K). Uptake kinetics are based on
!! the particle size, the gas-phase species diffusion coefficient and
!! molecular weight, and \f$N^{*}\f$, which is used to calculate the mass
!! accomodation coefficient. Details of the calculations can be found in:
!!
!! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
!! and condensed tropospheric aqueous mechanism and its application."
!! J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
!!
!! Input data for Phase transfer equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "PHASE_TRANSFER",
!!     "A" : 123.45,
!!     "C" : 123.45,
!!     "gas-phase species" : "my gas spec",
!!     "aerosol-phase species" : "my aero spec"
!!       ...
!!   }
!! \endcode
!! The key-value pairs \b gas-phase species, and \b aerosol-phase species are
!! required. Only one gas- and one aerosol-phase species are allowed per
!! phase-transfer reaction. Additionally, gas-phase species must include
!! parameters named "diffusion coeff", which specifies the diffusion
!! coefficient in (\f$m^2s^{-1}\f$), and "molecular weight", which specifies the molecular
!! weight of the species in (kg/mol). They may optionally include the
!! parameter "N star", which will be used to calculate the mass accomodation
!! coefficient. When this parameter is not included, the mass accomodation
!! coefficient is assumed to be 1.0.
!!
!! When \b A is not included, it is assumed to be 1.0, when \b C is not
!! included, it is assumed to be 0.0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_phase_transfer_t type and associated functions. 
module pmc_rxn_phase_transfer

  use pmc_constants,                        only: const
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg, &
                                                  string_t
  use pmc_rxn_data
  use pmc_chem_spec_data
  use pmc_property
  use pmc_phlex_state
  use pmc_aero_rep_data
  use pmc_aero_phase_data

  implicit none
  private

#define _del_H_ this%condensed_data_real(1)
#define _del_S_ this%condensed_data_real(2)
#define _Dg_ this%condensed_data_real(3)
#define _pre_c_rms_ this%condensed_data_real(4)
#define _A_ this%condensed_data_real(5)
#define _C_ this%condensed_data_real(6)
#define _c_rms_alpha_ this%condensed_data_real(7)
#define _equil_const_ this%condensed_data_real(8)
#define _CONV_ this%condensed_data_real(9)
#define _MW_ this%condensed_data_real(10)
#define _ug_m3_TO_ppm_ this%condensed_data_real(11)
#define _NUM_AERO_PHASE_ this%condensed_data_int(1)
#define _GAS_SPEC_ this%condensed_data_int(2)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 11
#define _AERO_SPEC_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _AERO_SPEC_ACT_COEFF_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_AERO_PHASE_+x)
#define _AERO_WATER_(x) this%condensed_data_int(_NUM_INT_PROP_+2*_NUM_AERO_PHASE_+x)
#define _AERO_WATER_ACT_COEFF_(x) this%condensed_data_int(_NUM_INT_PROP_+3*_NUM_AERO_PHASE_+x)
#define _AERO_PHASE_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+4*_NUM_AERO_PHASE_+x)
#define _AERO_REP_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+5*_NUM_AERO_PHASE_+x)
#define _DERIV_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+6*_NUM_AERO_PHASE_+x)
#define _JAC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+1+7*_NUM_AERO_PHASE_+x)

  public :: rxn_phase_transfer_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_phase_transfer_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_phase_transfer_t

  !> Constructor for rxn_phase_transfer_t
  interface rxn_phase_transfer_t
    procedure :: constructor
  end interface rxn_phase_transfer_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Phase transfer reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_phase_transfer_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information from reactant, product and reaction 
  !! property_t objects. This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_phase_transfer_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props
    character(len=:), allocatable :: key_name, spec_name, water_name, phase_name, &
            string_val
    integer(kind=i_kind) :: i_spec, i_qty, i_aero_rep, i_aero_phase, n_aero_ids
    integer(kind=i_kind) :: i_aero_id
    class(string_t), allocatable :: unique_spec_names(:), unique_water_names(:)
    integer(kind=i_kind), allocatable :: aero_spec_ids(:)
    integer(kind=i_kind), allocatable :: water_spec_ids(:)

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real, N_star

    ! Get the property set
    if (.not. associated(this%property_set)) call die_msg(318525776, &
            "Missing property set needed to initialize reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(448087197, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in phase-transfer reaction")

    ! Get the aerosol-phase species name
    key_name = "aerosol-phase species"
    call assert_msg(797902545, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase species in phase-transfer reaction")

    ! Get the aerosol-phase water name
    key_name = "aerosol-phase water"
    call assert_msg(386889865, &
            this%property_set%get_string(key_name, water_name), &
            "Missing aerosol-phase water in phase-transfer reaction")

    ! Check for aerosol representations
    call assert_msg(234155350, associated(aero_rep), &
            "Missing aerosol representation for phase transfer reaction")
    call assert_msg(207961800, size(aero_rep).gt.0, &
            "Missing aerosol representation for phase transfer reaction")
    
    ! Count the instances of this phase/species pair
    n_aero_ids = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Get the unique names in this aerosol representation for the 
      ! partitioning species and aerosol-phase water
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      unique_water_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)

      ! Skip aerosol representations that do not contain this phase
      if (.not.allocated(unique_spec_names)) cycle

      ! Check the size of the unique name lists
      call assert_msg(598091463, size(unique_spec_names).eq. &
              size(unique_water_names), "Missing species "// &
              spec_name//" or "//water_name//" in phase "//phase_name// &
              " or improper implementation of aerosol phase in aerosol "// &
              "representation")
 
      ! Add these instances to the list     
      n_aero_ids = n_aero_ids + size(unique_spec_names)
    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(_NUM_INT_PROP_ + 1 + n_aero_ids * 12))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set the number of aerosol-species instances
    _NUM_AERO_PHASE_ = n_aero_ids

    ! Get the properties required of the aerosol species
    call assert_msg(669162256, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "aerosol-phase species "//trim(spec_name))

    ! Get the aerosol species molecular weight
    key_name = "molecular weight"
    call assert_msg(209812557, spec_props%get_real(key_name, _MW_), &
            "Missing property 'MW' for aerosol species "//trim(spec_name)// &
            " required for phase-transfer reaction")

    ! Set the ug/m3 -> ppm conversion prefactor (multiply by T/P to get conversion)
    _CONV_ = const%univ_gas_const / _MW_

    ! Get the aerosol-phase water species
    key_name = "aerosol-phase water"
    call assert_msg(374667967, &
            this%property_set%get_string(key_name, water_name), &
            "Missing aerosol-phase water in phase-transfer reaction")

    ! Set the ids of each aerosol-phase species instance
    i_aero_id = 1
    do i_aero_rep = 1, size(aero_rep)
        
      ! Get the unique names in this aerosol representation for the 
      ! partitioning species and aerosol-phase water
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      unique_water_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)
     
      ! Add the species concentration and activity coefficient ids to
      ! the condensed data 
      do i_spec = 1, size(unique_spec_names)
        _AERO_SPEC_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_spec_names(i_spec)%string)
        _AERO_SPEC_ACT_COEFF_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%activity_coeff_state_id( &
              unique_spec_names(i_spec)%string)
        _AERO_WATER_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_water_names(i_spec)%string)
        _AERO_WATER_ACT_COEFF_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%activity_coeff_state_id( &
              unique_water_names(i_spec)%string)
        _AERO_PHASE_ID_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%phase_id(phase_name)
        _AERO_REP_ID_(i_aero_id) = i_aero_rep
        i_aero_id = i_aero_id + 1
      end do
    end do

    ! Get reaction parameters 
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, _A_)) then
      _A_ = 1.0
    end if
    key_name = "C"
    if (.not. this%property_set%get_real(key_name, _C_)) then
      _C_ = 0.0
    end if

    ! Get the gas-phase species and find the required species properties and index
    key_name = "gas-phase species"
    call assert_msg(847983010, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase species in phase-transfer reaction")

    ! Save the index of this species in the state variable array
    _GAS_SPEC_ = chem_spec_data%gas_state_id(spec_name)

    ! Make sure the species exists
    call assert_msg(751684145, _GAS_SPEC_.gt.0, &
            "Missing phase-transfer gas-phase species: "//spec_name)

    ! Get the required properties for the gas-phase species
    call assert_msg(757296139, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "gas-phase species "//trim(spec_name))

    ! Get N* to calculate the mass accomodation coefficient. If it is not
    ! present, set _del_H_ and _del_S_ to zero to indicate a mass accomodation
    ! coefficient of 1.0
    ! Mass accomodation equation is based on equations in:
    ! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
    ! and condensed tropospheric aqueous mechanism and its application."
    ! J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
    key_name = "N star"
    if (spec_props%get_real(key_name, N_star)) then     
      ! enthalpy change (kcal mol-1)
      _del_H_ = real(- 10.0d0*(N_star-1.0d0) + &
              7.53d0*(N_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0, kind=dp)
      ! entropy change (cal mol-1)
      _del_S_ = real(- 32.0d0*(N_star-1.0d0) + &
              9.21d0*(N_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0, kind=dp)
      ! Convert dH and dS to (J mol-1)
      _del_H_ = real(_del_H_ * 4184.0d0, kind=dp)
      _del_S_ = real(_del_S_ * 4.184d0, kind=dp)
    else
      _del_H_ = real(0.0, kind=dp)
      _del_S_ = real(0.0, kind=dp)
    end if

    ! Get the diffusion coefficient (m^2/s)
    key_name = "diffusion coeff"
    call assert_msg(100205531, spec_props%get_real(key_name, _Dg_), &
            "Missing diffusion coefficient for species "//spec_name)
    
    ! Calculate the constant portion of c_rms [m/(K^2*s)]
    key_name = "molecular weight"
    call assert_msg(469582180, spec_props%get_real(key_name, temp_real), &
            "Missing molecular weight for species "//spec_name)
    _pre_c_rms_ = sqrt(8.0*const%univ_gas_const/(const%pi*temp_real))

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_phase_transfer_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id

    expr = ""

  end function build_rate_const_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build time derivative expression
  function build_deriv_expr(this, rxn_id, spec_id, chem_spec_data) &
                  result (expr)

    !> Contribution to time derivative expression for species spec_id
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_phase_transfer_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _del_H_
#undef _del_S_
#undef _Dg_
#undef _pre_c_rms_
#undef _A_
#undef _C_
#undef _c_rms_alpha_
#undef _equil_const_
#undef _CONV_
#undef _MW_
#undef _ug_m3_TO_ppm_
#undef _NUM_AERO_PHASE_
#undef _GAS_SPEC_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _AERO_SPEC_
#undef _AERO_SPEC_ACT_COEFF_
#undef _AERO_WATER_
#undef _AERO_WATER_ACT_COEFF_
#undef _AERO_PHASE_ID_
#undef _AERO_REP_ID_
#undef _DERIV_ID_
#undef _JAC_ID_
end module pmc_rxn_phase_transfer
