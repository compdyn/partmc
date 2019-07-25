! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_HL_phase_transfer module.

!> \page phlex_rxn_HL_phase_transfer Phlexible Module for Chemistry: Henry's Law Phase-Transfer Reaction
!!
!! Henry's Law phase-trasfer reactions use equilibrium rate constants that
!! are calculated as:
!!
!! \f[
!!   HLC(298K) * e^{C({1/T-1/298})}
!! \f]
!!
!! where \f$HLC(298K)\f$ is the Henry's Law constant at 298 K
!! [\f$\mbox{\si{M Pa^{-1}}}\f$], \f$C\f$ is a constant [\f$\mbox{K}\f$] and
!! \f$T\f$ is the temperature [\f$\mbox{K}\f$]. Uptake kinetics are based on
!! the particle effective radius, \f$r_{eff}\f$ [\f$\mbox{m}\f$], the
!! condensing species gas-phase diffusion coefficient, \f$D_g\f$
!! [\f$\mbox{\si{\square\metre\per\second}}\f$], its molecular weight \f$MW\f$
!! [\f$\mbox{\si{\kilo\gram\per\mole}}\f$], and \f$N^{*}\f$, which is
!! used to calculate the mass accomodation coefficient.
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
!! [\f$\mbox{\si{\joule\per\kelvin\per\mole}}\f$]. The particle radius used
!! to calculate \f$k_{f}\f$ is the effective radius [\f$r_{eff}\f$], which is
!! taken as the "least-wrong" choice for condensation rates, as it is weighted
!! to surface area \cite Zender2002 .
!!
!! Input data for Phase transfer equations have the following format :
!! \code{.json}
!!   {
!!     "type" : "HL_PHASE_TRANSFER",
!!     "gas-phase species" : "my gas spec",
!!     "aerosol-phase species" : "my aero spec",
!!     "areosol phase" : "my aqueous phase",
!!     "aerosol-phase water" : "H2O_aq",
!!       ...
!!   }
!! \endcode
!! The key-value pairs \b gas-phase \b species, and \b aerosol-phase
!! \b species are required. Only one gas- and one aerosol-phase species are
!! allowed per phase-transfer reaction. Additionally, gas-phase species must
!! include parameters named \b HLC(298K) \b [M \b Pa-1], which is the Henry's
!! Law constant at 298 K, \b HLC \b exp \b factor \b [K], which is the
!! Henry's Law constant exponential factor "C", \b diffusion \b coeff \b [m2
!! \b s-1], which specifies the diffusion coefficient in
!! \f$\mbox{\si{\square\metre\per\second}}\f$, and \b molecular \b weight
!! \b [kg \b mol-1], which specifies the molecular weight of the species in
!! \f$\mbox{\si{\kilo\gram\per\mole}}\f$. They may optionally include the
!! parameter \b N \b star, which will be used to calculate the mass
!! accomodation coefficient. When this parameter is not included, the mass
!! accomodation coefficient is assumed to be 1.0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_HL_phase_transfer_t type and associated functions.
module pmc_rxn_HL_phase_transfer

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
#define A_ this%condensed_data_real(5)
#define C_ this%condensed_data_real(6)
#define C_AVG_ALPHA_ this%condensed_data_real(7)
#define EQUIL_CONST_ this%condensed_data_real(8)
#define CONV_ this%condensed_data_real(9)
#define MW_ this%condensed_data_real(10)
#define UGM3_TO_PPM_ this%condensed_data_real(11)
#define SMALL_NUMBER_ this%condensed_data_real(12)
#define NUM_AERO_PHASE_ this%condensed_data_int(1)
#define GAS_SPEC_ this%condensed_data_int(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 12
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+NUM_AERO_PHASE_+x)
#define PHASE_INT_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+6*NUM_AERO_PHASE_+x)
#define PHASE_REAL_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+7*NUM_AERO_PHASE_+x)
#define AERO_SPEC_(x) this%condensed_data_int(PHASE_INT_LOC_(x))
#define AERO_WATER_(x) this%condensed_data_int(PHASE_INT_LOC_(x)+1)
#define AERO_PHASE_ID_(x) this%condensed_data_int(PHASE_INT_LOC_(x)+2)
#define AERO_REP_ID_(x) this%condensed_data_int(PHASE_INT_LOC_(x)+3)
#define NUM_AERO_PHASE_JAC_ELEM_(x) this%condensed_data_int(PHASE_INT_LOC_(x)+4)
#define PHASE_JAC_ID_(x,s,e) this%condensed_data_int(PHASE_INT_LOC_(x)+5+s*NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define SMALL_WATER_CONC_(x) this%condensed_data_real(PHASE_REAL_LOC_(x))
#define FAST_FLUX_(x) this%condensed_data_real(PHASE_REAL_LOC_(x)+1)
#define AERO_ADJ_(x) this%condensed_data_real(PHASE_REAL_LOC_(x)+2)
#define EFF_RAD_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)+2+e)
#define NUM_CONC_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)+2+NUM_AERO_PHASE_JAC_ELEM_(x,e)+e)

  public :: rxn_HL_phase_transfer_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_HL_phase_transfer_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_HL_phase_transfer_t

  !> Constructor for rxn_HL_phase_transfer_t
  interface rxn_HL_phase_transfer_t
    procedure :: constructor
  end interface rxn_HL_phase_transfer_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Phase transfer reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_HL_phase_transfer_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)

    !> Reaction data
    class(rxn_HL_phase_transfer_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props
    character(len=:), allocatable :: key_name, spec_name, water_name, &
            phase_name
    integer(kind=i_kind) :: i_spec, i_aero_rep, n_aero_ids, i_aero_id, &
            n_aero_jac_elem, i_phase
    type(string_t), allocatable :: unique_spec_names(:), &
            unique_water_names(:)
    integer(kind=i_kind), allocatable :: phase_ids(:)
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

    ! Count the instances of this phase/species pair, and the number of
    ! Jacobian elements needed in calculations of mass, volume, etc.
    n_aero_ids = 0
    n_aero_jac_elem = 0
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

      ! Get the number of Jacobian elements for calculations of mass, volume,
      ! number, etc. for this partitioning into this phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(phase_name)
      do i_phase = 1, size(phase_ids)
        n_aero_jac_elem = n_aero_jac_elem + &
                2*aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_phase))
      end do

      deallocate(unique_spec_names)
      deallocate(unique_water_names)

    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + 2 + n_aero_ids * 13 + &
                                      n_aero_jac_elem * 2))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + n_aero_ids * 3 + &
                                      n_aero_jac_elem * 2))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set the number of aerosol-species instances
    NUM_AERO_PHASE_ = n_aero_ids

    ! Get the properties required of the aerosol species
    call assert_msg(669162256, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "aerosol-phase species "//trim(spec_name))

    ! Get the aerosol species molecular weight
    key_name = "molecular weight [kg mol-1]"
    call assert_msg(209812557, spec_props%get_real(key_name, MW_), &
            "Missing property 'MW' for aerosol species "//trim(spec_name)// &
            " required for phase-transfer reaction")

    ! Set the ug/m3 -> ppm conversion prefactor (multiply by T/P to get
    ! conversion)
    CONV_ = const%univ_gas_const / MW_ / 1.0d3

    ! Get the aerosol-phase water species
    key_name = "aerosol-phase water"
    call assert_msg(374667967, &
            this%property_set%get_string(key_name, water_name), &
            "Missing aerosol-phase water in phase-transfer reaction")

    ! Set the ids of each aerosol-phase species instance
    i_aero_id = 1
    PHASE_INT_LOC_(i_aero_id)  = NUM_INT_PROP_+8*NUM_AERO_PHASE_+3
    PHASE_REAL_LOC_(i_aero_id) = NUM_REAL_PROP_+1
    do i_aero_rep = 1, size(aero_rep)

      ! Get the unique names in this aerosol representation for the
      ! partitioning species and aerosol-phase water
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      unique_water_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)

      ! Get the phase ids for this aerosol phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(phase_name)

      ! Add the species concentration and activity coefficient ids to
      ! the condensed data, and set the number of jacobian elements for
      ! the aerosol representations and the locations of the real data
      do i_spec = 1, size(unique_spec_names)
        NUM_AERO_PHASE_JAC_ELEM_(i_aero_id) = &
              2*aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_spec))
        AERO_SPEC_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_spec_names(i_spec)%string)
        AERO_WATER_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_water_names(i_spec)%string)
        AERO_PHASE_ID_(i_aero_id) = phase_ids(i_spec)
        AERO_REP_ID_(i_aero_id) = i_aero_rep
        i_aero_id = i_aero_id + 1
        if (i_aero_id .le. NUM_AERO_PHASE_) then
          PHASE_INT_LOC_(i_aero_id)  = PHASE_INT_LOC_(i_aero_id - 1) + 5 + &
                                     2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
          PHASE_REAL_LOC_(i_aero_id) = PHASE_REAL_LOC_(i_aero_id - 1) + 3 + &
                                     2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
        end if
      end do

      deallocate(unique_spec_names)
      deallocate(unique_water_names)

    end do

    ! Get the gas-phase species and find the required species properties and
    ! index
    key_name = "gas-phase species"
    call assert_msg(847983010, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase species in phase-transfer reaction")

    ! Save the index of this species in the state variable array
    GAS_SPEC_ = chem_spec_data%gas_state_id(spec_name)

    ! Make sure the species exists
    call assert_msg(955306778, GAS_SPEC_.gt.0, &
            "Missing phase-transfer gas-phase species: "//spec_name)

    ! Get the required properties for the gas-phase species
    call assert_msg(757296139, &
            chem_spec_data%get_property_set(spec_name, spec_props), &
            "Missing properties required for phase-transfer of "// &
            "gas-phase species "//trim(spec_name))

    ! Get Henry's Law constant parameters
    key_name = "HLC(298K) [M Pa-1]"
    call assert_msg(637925661, spec_props%get_real(key_name, A_), &
                    "Missing Henry's Law constant at 298 K for "// &
                    "partitioning species "//trim(spec_name))
    A_ = A_ * 101325.0 * 1.0d-6; ! save A in (M/ppm)

    key_name = "HLC exp factor [K]"
    call assert_msg(801365019, spec_props%get_real(key_name, C_), &
                    "Missing Henry's Law constant exponential factor "// &
                    "for partitioning species "//trim(spec_name))

    ! Get N* to calculate the mass accomodation coefficient. If it is not
    ! present, set DELTA_H_ and DELTA_S_ to zero to indicate a mass accomodation
    ! coefficient of 1.0
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
    key_name = "diffusion coeff [m2 s-1]"
    call assert_msg(100205531, spec_props%get_real(key_name, DIFF_COEFF_), &
            "Missing diffusion coefficient for species "//spec_name)

    ! Calculate the constant portion of c_rms [m/(K^2*s)]
    key_name = "molecular weight [kg mol-1]"
    call assert_msg(469582180, spec_props%get_real(key_name, temp_real), &
            "Missing molecular weight for species "//spec_name)
    PRE_C_AVG_ = sqrt(8.0*const%univ_gas_const/(const%pi*temp_real))

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_HL_phase_transfer_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_HL_phase_transfer
