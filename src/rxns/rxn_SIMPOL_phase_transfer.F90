! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_SIMPOL_phase_transfer module.

!> \page camp_rxn_SIMPOL_phase_transfer CAMP: SIMPOL.1 Phase-Transfer Reaction
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
!!   \frac{\alpha}{1-\alpha} = e^{\frac{-\Delta G^{*}}{RT}}
!! \f]
!! If \f$\Delta H\f$ and \f$\Delta S\f$ are not provided, the mass accomodation
!! coefficient is assumed to be 0.1 (\cite Zaveri2008).
!!
!! Condensation rate constants are calculated as:
!! \f[
!!   k_{f} = (\frac{r^2}{3D_g} + \frac{4r}{3 \langle c \rangle \alpha})^{-1}
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
!!     "aerosol-phase activity coefficient" : "my aero act coeff"
!!     "B" : [ 123.2e3, -41.24, 2951.2, -1.245e-4 ]
!!       ...
!!   }
!! \endcode
!! The key-value pairs \b gas-phase \b species, \b aerosol \b phase and
!! \b aerosol-phase \b species are required. Only one gas- and one
!! aerosol-phase species are allowed per phase-transfer reaction. The
!! key-value pair \b aerosol-phase \b activity \b coefficient is optional.
!! When it is included its value must be the name of an species of type
!! \c ACTIVITY_COEFF that is present in the specified aerosol phase. When
!! it is not included, activity coefficients are assume to be 1.0.
!!
!! Gas-phase species must include parameters named
!! \b diffusion \b coeff \b [m2 \b s-1], which specifies the diffusion
!! coefficient in \f$\mbox{\si{\square\metre\per\second}}\f$, and \b molecular
!! \b weight \b [kg \b mol-1], which specifies the molecular weight of the
!! species in \f$\mbox{\si{\kilo\gram\per\mole}}\f$. They may optionally
!! include the parameter \b N \b star, which will be used to calculate th
!! mass accomodation coefficient. When this parameter is not included, the
!! mass accomodation coefficient is assumed to be 1.0.
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
  use pmc_camp_state
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
#define CONV_ this%condensed_data_real(9)
#define MW_ this%condensed_data_real(10)
#define NUM_AERO_PHASE_ this%condensed_data_int(1)
#define GAS_SPEC_ this%condensed_data_int(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 10
#define NUM_ENV_PARAM_ 4
#define AERO_SPEC_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define AERO_ACT_ID_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_AERO_PHASE_+x)
#define AERO_PHASE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+2*NUM_AERO_PHASE_+x)
#define AERO_REP_ID_(x) this%condensed_data_int(NUM_INT_PROP_+3*NUM_AERO_PHASE_+x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+4*NUM_AERO_PHASE_+x)
#define GAS_ACT_JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+5*NUM_AERO_PHASE_+x)
#define AERO_ACT_JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+6*NUM_AERO_PHASE_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+7*NUM_AERO_PHASE_+x)
#define PHASE_INT_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+10*NUM_AERO_PHASE_+x)
#define PHASE_REAL_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+11*NUM_AERO_PHASE_+x)
#define NUM_AERO_PHASE_JAC_ELEM_(x) this%condensed_data_int(PHASE_INT_LOC_(x))
#define PHASE_JAC_ID_(x,s,e) this%condensed_data_int(PHASE_INT_LOC_(x)+(s-1)*NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define EFF_RAD_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+e)
#define NUM_CONC_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define MASS_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+2*NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define MW_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+3*NUM_AERO_PHASE_JAC_ELEM_(x)+e)

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
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_SIMPOL_phase_transfer_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props, b_params
    character(len=:), allocatable :: key_name, gas_spec_name, aero_spec_name
    character(len=:), allocatable :: phase_name, act_name, error_msg
    integer(kind=i_kind) :: i_spec, i_aero_rep, n_aero_ids, i_aero_id
    integer(kind=i_kind) :: i_phase, n_aero_jac_elem, tmp_size
    type(string_t), allocatable :: unique_spec_names(:), unique_act_names(:)
    integer(kind=i_kind), allocatable :: phase_ids(:)
    real(kind=dp) :: temp_real, N_star
    logical :: has_act_coeff

    ! Get the property set
    if (.not. associated(this%property_set)) call die_msg(382913491, &
            "Missing property set needed to initialize reaction")

    ! Get the gas-phase species name
    key_name = "gas-phase species"
    call assert_msg(740333884, &
            this%property_set%get_string(key_name, gas_spec_name), &
            "Missing gas-phase species in SIMPOL.1 phase transfer reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(325074932, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in SIMPOL.1 phase-transfer reaction")

    ! Get the aerosol-phase species name
    key_name = "aerosol-phase species"
    call assert_msg(988456388, &
            this%property_set%get_string(key_name, aero_spec_name), &
            "Missing aerosol-phase species in SIMPOL.1 phase-transfer "// &
            "reaction")

    ! Set up a general error message
    error_msg = " for SIMPOL.1 phase transfer of gas species '"// &
                gas_spec_name//"' to aerosol-phase species '"// &
                aero_spec_name//"' in phase '"//phase_name//"'"

    ! Get the aerosol-phase activity coeffcient name
    key_name = "aerosol-phase activity coefficient"
    has_act_coeff = this%property_set%get_string(key_name, act_name)

    ! Check for aerosol representations
    call assert_msg(260518827, associated(aero_rep), &
            "Missing aerosol representation"//error_msg)
    call assert_msg(590304021, size(aero_rep).gt.0, &
            "Missing aerosol representation"//error_msg)

    ! Count the instances of this phase/species pair
    n_aero_ids = 0
    n_aero_jac_elem = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Get the unique names in this aerosol representation for the
      ! partitioning species
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = aero_spec_name)

      ! Skip aerosol representations that do not contain this phase
      if (.not.allocated(unique_spec_names)) cycle

      ! Add these instances to the list
      n_aero_ids = n_aero_ids + size(unique_spec_names)

      ! Get the number of Jacobian elements for calculations of mass, volume,
      ! number, etc. for this partitioning into this phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(phase_name)
      do i_phase = 1, size(phase_ids)
        n_aero_jac_elem = n_aero_jac_elem + &
                aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_phase))
      end do

      deallocate(unique_spec_names)

    end do

    call assert_msg(314000134, n_aero_ids.gt.0, &
                    "Aerosol species not found"//error_msg)

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + 2 + n_aero_ids * 13 + &
                                     n_aero_jac_elem * 2))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + n_aero_jac_elem * 4))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Set the number of aerosol-species instances
    NUM_AERO_PHASE_ = n_aero_ids

    ! Get the properties required of the aerosol species
    call assert_msg(162662115, &
            chem_spec_data%get_property_set(aero_spec_name, spec_props), &
            "Missing properties"//error_msg)

    ! Get the aerosol species molecular weight
    key_name = "molecular weight [kg mol-1]"
    call assert_msg(839930958, spec_props%get_real(key_name, MW_), &
            "Missing property 'MW'"//error_msg)

    ! Set the kg/m3 -> ppm conversion prefactor (multiply by T/P to get
    ! conversion)
    ! (ppm_x*Pa_air*m^3/K/kg_x) = Pa_air*m^3/mol_air/K * mol_x/kg_x *
    !                   1.0e6ppm_x*mol_air/mol_x
    CONV_ = const%univ_gas_const / MW_ * 1.0e6

    ! Set the ids of each aerosol-phase species instance
    i_aero_id = 1
    PHASE_INT_LOC_(i_aero_id)  = NUM_INT_PROP_+12*NUM_AERO_PHASE_+3
    PHASE_REAL_LOC_(i_aero_id) = NUM_REAL_PROP_+1
    do i_aero_rep = 1, size(aero_rep)

      ! Get the unique names in this aerosol representation for the
      ! partitioning species
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = aero_spec_name)

      ! Find the corresponding activity coefficients, if specified
      if (has_act_coeff) then
        unique_act_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = act_name)
        call assert_msg(236251734, size(unique_act_names).eq. &
                        size(unique_spec_names), &
                        "Mismatch of SIMPOL species and activity coeffs"// &
                        error_msg)
      end if

      ! Get the phase ids for this aerosol phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(phase_name)

      ! Add the species concentration and activity coefficient ids to
      ! the condensed data, and set the number of Jacobian elements for
      ! the aerosol representations and the locations of the real data
      do i_spec = 1, size(unique_spec_names)
        NUM_AERO_PHASE_JAC_ELEM_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_spec))
        AERO_SPEC_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_spec_names(i_spec)%string)
        if (has_act_coeff) then
          AERO_ACT_ID_(i_aero_id) = &
              aero_rep(i_aero_rep)%val%spec_state_id( &
              unique_act_names(i_spec)%string)
        else
          AERO_ACT_ID_(i_aero_id) = -1
        end if
        AERO_PHASE_ID_(i_aero_id) = phase_ids(i_spec)
        AERO_REP_ID_(i_aero_id) = i_aero_rep
        i_aero_id = i_aero_id + 1
        if (i_aero_id .le. NUM_AERO_PHASE_) then
          PHASE_INT_LOC_(i_aero_id)  = PHASE_INT_LOC_(i_aero_id - 1) + 1 + &
                                     2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
          PHASE_REAL_LOC_(i_aero_id) = PHASE_REAL_LOC_(i_aero_id - 1) + &
                                     4*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
        end if
      end do

      deallocate(unique_spec_names)

    end do

    ! Get the SIMPOL.1 parameters
    key_name = "B"
    call assert_msg(882881186, &
            this%property_set%get_property_t(key_name, b_params), &
            "Missing 'B' parameters"//error_msg)
    call assert_msg(654885723, b_params%size().eq.4, &
            "Incorrect number of 'B' parameters"//error_msg)
    call b_params%iter_reset()
    call assert_msg(694024883, b_params%get_real(val = B1_), &
            "Got non-real 'B1' parameter"//error_msg)
    call b_params%iter_next()
    call assert_msg(231316411, b_params%get_real(val = B2_), &
            "Got non-real 'B2' parameter"//error_msg)
    call b_params%iter_next()
    call assert_msg(126167907, b_params%get_real(val = B3_), &
            "Got non-real 'B3' parameter"//error_msg)
    call b_params%iter_next()
    call assert_msg(573535753, b_params%get_real(val = B4_), &
            "Got non-real 'B4' parameter"//error_msg)

    ! Save the index of the gas-phase species in the state variable array
    GAS_SPEC_ = chem_spec_data%gas_state_id(gas_spec_name)

    ! Make sure the species exists
    call assert_msg(551477581, GAS_SPEC_.gt.0, &
            "Missing gas-phase species"//error_msg)

    ! Get the required properties for the gas-phase species
    call assert_msg(611221674, &
            chem_spec_data%get_property_set(gas_spec_name, spec_props), &
            "Missing properties for gas-phase species"//error_msg)

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
    key_name = "diffusion coeff [m2 s-1]"
    call assert_msg(948176709, spec_props%get_real(key_name, DIFF_COEFF_), &
            "Missing diffusion coefficient"//error_msg)

    ! Calculate the constant portion of c_rms [m /( K^(1/2) * s )]
    key_name = "molecular weight [kg mol-1]"
    call assert_msg(272813400, spec_props%get_real(key_name, temp_real), &
            "Missing molecular weight"//error_msg)
    PRE_C_AVG_ = sqrt(8.0*const%univ_gas_const/(const%pi*temp_real))

    ! Check the sizes of the data arrays
    tmp_size = PHASE_INT_LOC_(i_aero_id - 1) + 1 + &
               2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1) - 1
    call assert_msg(625802519, size(this%condensed_data_int) .eq. tmp_size, &
                    "int array size mismatch"//error_msg)
    tmp_size = PHASE_REAL_LOC_(i_aero_id - 1) + &
               4*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1) - 1
    call assert_msg(391089510, size(this%condensed_data_real) .eq. tmp_size, &
                    "real array size mismatch"//error_msg)

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

end module pmc_rxn_SIMPOL_phase_transfer
