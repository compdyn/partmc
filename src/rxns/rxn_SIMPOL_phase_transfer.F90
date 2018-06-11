! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_SIMPOL_phase_transfer module.

!> \page phlex_rxn_SIMPOL_phase_transfer Phlexible Module for Chemistry: SIMPOL.1 Phase-Transfer Reaction
!!
!! SIMPOL phase transfer reactions are based on the SIMPOL.1 model 
!! calculations of vapor pressure described by Pankow and Asher (2008)
!! \cite Pankow2008. 
!!
!! Mass accomodation coefficients and condensation rate constants are
!! calculated using the method of Ervans et al. (2003) \cite Ervens2003 and
!! references therein. Mass accomodation coefficients (\f$\alpha\f$) are
!! calculated as:
!!
!! \f[
!!   \Delta H_{obs} = -10 \times (N^*-1) + 7.53 \times (N^{*2/3}-1) - 0.1 \times 10 (\mbox{\si{kcal.M^{-1}}})
!! \f]
!! \f[
!!   \Delta S_{obs} = -13 \times (N^*-1) - 19 \times (N^*-1) + 9.21 \times (N^{*2/3}-1) - 0.1 \times 13 (\mbox{\si{cal.M^{-1}.K^{-1}}})
!! \f]
!! \f[
!!   \frac{\alpha}{1-\alpha} = e^{\frac{-\Delta G^{\*}}{RT}}
!! \f]
!!
!! Condensation rate constants are calculated as:
!! \f[
!!   k_{f} = (\frac{r^2}{3D_g} + \frac{4r}{3 \langle c \rangle \alpha})
!! \f]
!! where \f$r\f$ is the particle radius (\f$\mbox{m}\f$) and 
!! \f$\langle c \rangle \f$ is the mean speed of the gas-phase molecules:
!! \f[
!!   \langle c \rangle = \sqrt{\frac{8RT}{\pi MW}}
!! \f]
!! where \f$R\f$ is the ideal gas constant
!! (\f$\mbox{\si{\joule\per\kelvin\per\mole}}\f$). The particle radius used
!! to calculate \f$k_{f}\f$ is the effective radius (\f$r_{eff}\f$), which is
!! taken as the "least-wrong" choice for condensation rates, as it is weighted
!! to surface area \cite Zender2002 .
!!
!! Input data for SIMPOL phase transfer reactions have the following format :
!! \code{.json}
!!   {
!!     "type" : "SIMPOL_PHASE_TRANSFER",
!!     "gas-phase species" : "my gas spec",
!!     "aerosol phase" : "my aero phase",
!!     "aerosol-phase species" : "my aero spec",
!!     "B" : [ 123.2e3, -41.24, 2951.2, -1.245e-4 ]
!!       ...
!!   }
!! \endcode
!! The key-value pairs \b gas-phase \b species, \b aerosol \b phase and 
!! \b aerosol-phase \b species are required. Only one gas- and one 
!! aerosol-phase species are allowed per phase-transfer reaction. 
!! Additionally, gas-phase species must include parameters named 
!! \b diffusion \b coeff, which specifies the diffusion coefficient in
!! \f$\mbox{\si{\square\metre\per\second}}\f$, and \b molecular \b weight,
!! which specifies the molecular weight of the species in 
!! \f$\mbox{\si{\kilo\gram\per\mole}}\f$. They may optionally include the
!! parameter \b N \b star, which will be used to calculate the mass
!! accomodation coefficient. When this parameter is not included, the mass
!! accomodation coefficient is assumed to be 1.0.
!!
!! The key-value pair \b B is also required and must have a value of an array
!! of exactly four members that specifies the SIMPOL parameters for the
!! partitioning species. The \b B parameters can be obtained by summing the
!! contributions of each functional group present in the partitioning species
!! to the overall \f$B_{n,i}\f$ for species \f$i\f$, such that:
!! \f[
!!   B_{n,i} = \sum_{k} \nu_{k,i} B_{n,k} \forall n \in [1...4]
!! \f]
!! where \f$\nu_{k,i}\f$ is the number of functional groups \f$k\f$ in species
!! \f$i\f$ and the parameters \f$B_{n,k}\f$ for each functional group \f$k\f$
!! can be found in table 5 of Pankow and Asher (2008) \cite Pankow2008.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_SIMPOL_phase_transfer_t type and associated functions. 
module pmc_rxn_SIMPOL_phase_transfer

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, &
                                                  die_msg, string_t

  implicit none
  private

#define DELTA_H_ this%condensed_data_real(1)
#define DELTA_S_ this%condensed_data_real(2)
#define DIFF_COEFF_ this%condensed_data_real(3)
#define PRE_C_AVG_ this%condensed_data_real(4)
#define B1_ this%condensed_data_real(5)
#define B2_ this%condensed_data_real(6)
#define B3_ this%condensed_data_real(7)
#define B4_ this%condensed_data_real(8)
#define C_AVG_ALHPA_ this%condensed_data_real(9)
#define EQUIL_CONST_ this%condensed_data_real(10)
#define CONV_ this%condensed_data_real(11)
#define MW_ this%condensed_data_real(12)
#define UGM3_TO_PPM_ this%condensed_data_real(13)
#define NUM_AERO_PHASE_ this%condensed_data_int(1)
#define GAS_SPEC_ this%condensed_data_int(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 13
#define AERO_SPEC_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define AERO_ACT_ID_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_AERO_PHASE_+x)
#define AERO_PHASE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+2*NUM_AERO_PHASE_+x)
#define AERO_REP_ID_(x) this%condensed_data_int(NUM_INT_PROP_+3*NUM_AERO_PHASE_+x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+4*NUM_AERO_PHASE_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+5*NUM_AERO_PHASE_+x)

  public :: rxn_SIMPOL_phase_transfer_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_SIMPOL_phase_transfer_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_SIMPOL_phase_transfer_t

  !> Constructor for rxn_SIMPOL_phase_transfer_t
  interface rxn_SIMPOL_phase_transfer_t
    procedure :: constructor
  end interface rxn_SIMPOL_phase_transfer_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Phase transfer reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_SIMPOL_phase_transfer_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_SIMPOL_phase_transfer_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, b_params
    character(len=:), allocatable :: key_name, spec_name, phase_name, &
            string_val
    integer(kind=i_kind) :: i_spec, i_qty, i_aero_rep, i_aero_phase, &
            n_aero_ids, i_aero_id, temp_int
    type(string_t), allocatable :: unique_spec_names(:)
    integer(kind=i_kind), allocatable :: aero_spec_ids(:), phase_ids(:)
    real(kind=dp) :: temp_real, N_star

    ! Get the property set
    if (.not. associated(this%property_set)) call die_msg(382913491, &
            "Missing property set needed to initialize reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(325074932, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in phase-transfer reaction")

    ! Get the aerosol-phase species name
    key_name = "aerosol-phase species"
    call assert_msg(988456388, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase species in phase-transfer reaction")

    ! Check for aerosol representations
    call assert_msg(260518827, associated(aero_rep), &
            "Missing aerosol representation for phase transfer reaction")
    call assert_msg(590304021, size(aero_rep).gt.0, &
            "Missing aerosol representation for phase transfer reaction")
    
    ! Count the instances of this phase/species pair
    n_aero_ids = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Get the unique names in this aerosol representation for the 
      ! partitioning species 
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)

      ! Skip aerosol representations that do not contain this phase
      if (.not.allocated(unique_spec_names)) cycle

      ! Add these instances to the list     
      n_aero_ids = n_aero_ids + size(unique_spec_names)

      deallocate(unique_spec_names)

    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + 2 + n_aero_ids * 8))
    allocate(this%condensed_data_real(NUM_REAL_PROP_))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set the number of aerosol-species instances
    NUM_AERO_PHASE_ = n_aero_ids

    ! Get the properties required of the aerosol species
    call assert_msg(162662115, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "aerosol-phase species "//trim(spec_name))

    ! Get the aerosol species molecular weight
    key_name = "molecular weight"
    call assert_msg(839930958, spec_props%get_real(key_name, MW_), &
            "Missing property 'MW' for aerosol species "//trim(spec_name)// &
            " required for phase-transfer reaction")

    ! Set the ug/m3 -> ppm conversion prefactor (multiply by T/P to get 
    ! conversion)
    ! (ppm_x*Pa_air*m^3/K/ug_x) = Pa_air*m^3/mol_air/K * mol_x/kg_x * 
    !                   1.0e-9kg_x/ug_x * 1.0e6ppm_x*mol_air/mol_x
    CONV_ = const%univ_gas_const / MW_ / 1.0e3

    ! Set the ids of each aerosol-phase species instance
    i_aero_id = 1
    do i_aero_rep = 1, size(aero_rep)
        
      ! Get the unique names in this aerosol representation for the 
      ! partitioning species
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
     
      ! Get the phase ids for this aerosol phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(phase_name)

      ! Add the species concentration and activity coefficient ids to
      ! the condensed data 
      do i_spec = 1, size(unique_spec_names)
        AERO_SPEC_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_spec_names(i_spec)%string)
        AERO_PHASE_ID_(i_aero_id) = phase_ids(i_spec)
        AERO_REP_ID_(i_aero_id) = i_aero_rep
        i_aero_id = i_aero_id + 1
      end do

      deallocate(unique_spec_names)

    end do

    ! Get the SIMPOL.1 parameters 
    key_name = "B"
    call assert_msg(882881186, &
            this%property_set%get_property_t(key_name, b_params), &
            "Missing 'B' parameters for SIMPOL.1 vapor pressure "// &
            "calculation in phase transfer reaction of "//spec_name)
    call assert_msg(654885723, b_params%size().eq.4, &
            "Incorrect number of 'B' parameters for SIMPOL.1 vapor "// &
            "pressure calculation in phase transfer reactions of "// &
            spec_name)
    call b_params%iter_reset()
    call assert_msg(694024883, b_params%get_real(val = B1_), &
            "Got non-real 'B' parameter in phase transfer reaction of "// &
            spec_name)
    call b_params%iter_next()
    call assert_msg(231316411, b_params%get_real(val = B2_), &
            "Got non-real 'B' parameter in phase transfer reaction of "// &
            spec_name)
    call b_params%iter_next()
    call assert_msg(126167907, b_params%get_real(val = B3_), &
            "Got non-real 'B' parameter in phase transfer reaction of "// &
            spec_name)
    call b_params%iter_next()
    call assert_msg(573535753, b_params%get_real(val = B4_), &
            "Got non-real 'B' parameter in phase transfer reaction of "// &
            spec_name)

    ! Get the gas-phase species and find the required species properties and
    ! index
    key_name = "gas-phase species"
    call assert_msg(740333884, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase species in phase-transfer reaction")

    ! Save the index of this species in the state variable array
    GAS_SPEC_ = chem_spec_data%gas_state_id(spec_name)

    ! Make sure the species exists
    call assert_msg(551477581, GAS_SPEC_.gt.0, &
            "Missing phase-transfer gas-phase species: "//spec_name)

    ! Get the required properties for the gas-phase species
    call assert_msg(611221674, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "gas-phase species "//trim(spec_name))

    ! Get N* to calculate the mass accomodation coefficient. If it is not
    ! present, set DELTA_H_ and DELTA_S_ to zero to indicate a mass
    ! accomodation coefficient of 1.0
    ! Mass accomodation equation is based on equations in:
    ! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
    ! and condensed tropospheric aqueous mechanism and its application."
    ! J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
    key_name = "N star"
    if (spec_props%get_real(key_name, N_star)) then     
      ! enthalpy change (kcal mol-1)
      DELTA_H_ = real(- 10.0d0*(N_star-1.0d0) + &
              7.53d0*(N_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0, kind=dp)
      ! entropy change (cal mol-1)
      DELTA_S_ = real(- 32.0d0*(N_star-1.0d0) + &
              9.21d0*(N_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0, kind=dp)
      ! Convert dH and dS to (J mol-1)
      DELTA_H_ = real(DELTA_H_ * 4184.0d0, kind=dp)
      DELTA_S_ = real(DELTA_S_ * 4.184d0, kind=dp)
    else
      DELTA_H_ = real(0.0, kind=dp)
      DELTA_S_ = real(0.0, kind=dp)
    end if

    ! Get the diffusion coefficient (m^2/s)
    key_name = "diffusion coeff"
    call assert_msg(948176709, spec_props%get_real(key_name, DIFF_COEFF_), &
            "Missing diffusion coefficient for species "//spec_name)
    
    ! Calculate the constant portion of c_rms [m/(K^2*s)]
    key_name = "molecular weight"
    call assert_msg(272813400, spec_props%get_real(key_name, temp_real), &
            "Missing molecular weight for species "//spec_name)
    PRE_C_AVG_ = sqrt(8.0*const%univ_gas_const/(const%pi*temp_real*1.0e3))

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_SIMPOL_phase_transfer_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef DELTA_H_
#undef DELTA_S_
#undef DIFF_COEFF_
#undef PRE_C_AVG_
#undef B1_
#undef B2_
#undef B3_
#undef B4_
#undef C_AVG_ALHPA_
#undef EQUIL_CONST_
#undef CONV_
#undef MW_
#undef UGM3_TO_PPM_
#undef NUM_AERO_PHASE_
#undef GAS_SPEC_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef AERO_SPEC_
#undef AERO_ACT_ID_
#undef AERO_PHASE_ID_
#undef AERO_REP_ID_
#undef DERIV_ID_
#undef JAC_ID_
end module pmc_rxn_SIMPOL_phase_transfer
