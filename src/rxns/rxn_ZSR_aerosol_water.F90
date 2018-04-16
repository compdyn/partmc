! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_ZSR_aerosol_water module.

!> \page phlex_rxn_ZSR_aerosol_water Phlexible Mechanism for Chemistry: ZSR Aerosol Water Reaction
!!
!! ZSR aerosol water reactions calculate equilibrium aerosol water content
!! based on the Zdanovski-Stokes-Robinson mixing rule \cite{Stokes1966,
!! Jacobson1996} in the following generalized format:
!!
!! \f[
!!   W = \sum\limits_{i=0}^{n}\frac{1000 M_i}{MW_i m_{i}(a_w)}
!! \f]
!!
!! where \f$M\f$ is the concentration of binary electrolyte \f$i\f$ 
!! (\f$\mu g m^{-3}\f$) with molecular weight \f$MW_i\f$ (g/mol) and 
!! molality \f$m_{i}\f$ at a given water activity \f$a_w\f$ (RH; 0-1)
!! contributing to the total aerosol water content \f$W\f$ 
!! (\f$\mu g m^{-3}\f$).
!!
!! Input data for ZSR aerosol water equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "ZSR_AEROSOL_WATER",
!!     "aerosol phase" : "my aero phase",
!!     "gas-phase water" : "H2O",
!!     "aerosol-phase water" : "H2O_aq",
!!     "ion pairs" : {
!!       "Na2SO4" : {
!!         "type" : "JACOBSON",
!!         "ions" : {
!!           "Nap" : { "qty" : 2 },
!!           "SO4mm" : {}
!!         },
!!         "Y_j" : [-3.295311e3, 3.188349e4, -1.305168e5, 2.935608e5],
!!         "low RH" : 51.0
!!       },
!!       "H2SO4" : {
!!         "type" : "EQSAM",
!!         "ions" : {
!!           "SO4mm" : {}
!!         },
!!         "NW" : 4.5,
!!         "ZW" : 0.5,
!!         "MW" : 98.0
!!       }
!!       ...
!!     }
!!   }
!! \endcode
!! The key-value pair \b "aerosol phase" is required to specify the aerosol
!! phase for which to calculate water content. Key-value pairs 
!! \b "gas-phase water" and \b "aerosol-phase water" must also be present and
!! specify the names for the water species in each phase. The final required
!! key-value pair is \b "ion pairs" which should contain a set of key-value
!! pairs where the key of each member of the set is the name of a binary
!! electrolyte and the contents contain parameters required to estimate
!! the contribution of the this electrolyte to total aerosol water. The
!! name of the electrolyte may or may not refer to an actual aerosol-phase
!! species.
!!
!! Each binary electrolyte must include a \b "type" that refers to a method
!! of calculating ion-pair contributions to aerosol water. Valid values for
!! \b "type" are "JACOBSON" and "EQSAM". These are described next.
!!
!! Aerosol water from ion pairs with type "JACOBSON" use equations (28) and
!! (29) in Jacobson et al. \cite{Jacobson1996} where experimentally determined
!! binary solution molalities are fit to a polynomial as:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = Y_0 + Y_1 a_w + Y_2 a_w^2 + Y_3 a_w^3 + ...,
!! \f]
!!
!! where \f$Y_j\f$ are the fitting parameters. Thus, \f$m_i(a_w)\f$ is
!! calculated at each time step, assuming constant \f$a_w\f$. These values
!! must be included in a key-value pair \b "Y_j" whose value is an array
!! with the \f$Y_j\f$ parameters. The size of the array corresponds to the
!! order of the polynomial equation, which must be greater than 1. The
!! key-value pair \b "low RH" is required to specify the lowest RH for which
!! this fit is valid. This value for RH will be used for all lower RH in
!! calculations of \f$m_i(a_w)\f$ as per Jacobson et al. \cite{1996}.
!!
!! The key-value pair "ions" must contain the set of ions this binary
!! electrolyte includes. Each species must correspond to a species present in
!! \b "aerosol phase" and  have a \b "charge" parameter that specifies their
!! charge (uncharged species are not permitted in this set) and a 
!! \b "molecular weight" (g/mol) property. Ions without a \b "qty" specified
!! are assumed to appear once in the binary electrolyte. The total
!! molecular weight for the binary electrolye \f$MW_i\f$ is calculated as a 
!! sum of its ionic components, and the ion species concentrations are used
!! to determine the \f$M_i\f$ during integration.
!!
!! For the above example, the following input data should be present:
!! \code{.json}
!! {
!!   "name" : "H2O",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "GAS",
!! },  
!! {
!!   "name" : "Nap",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : 1,
!!   "molecular weight" : 22.9898
!! },
!! {
!!   "name" : "SO4mm",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -2
!!   "molecular weight" : 96.06
!! },
!! {
!!   "name" : "my aero phase",
!!   "type" : "AERO_PHASE",
!!   "species" : ["Nap", "SO4mm", H2O_aq"]
!! }
!!
!! Aerosol water from ion pairs with type "EQSAM" use the parameterization of
!! Metzger et al. \cite{Metzget2002} for aerosol water content:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = (NW_i MW_{H2O}/MW_i 1/(a_w-1))^{ZW_i}
!! \f]
!!
!! where \f$NW_i\f$ and \f$ZW_i\f$ are fitting parameters \cite{Metger2002},
!! and must be provided in key-value pairs \b "NW" and \b "ZW", along with the
!! binary electrolyte molecular weight \b "MW" (g/mol). The key-value pair
!! \b "ions" must contain a set of ions that can be summed to calculate
!! \f$M_i\f$ at runtime.
!!
!! TODO Find a way to incorporate the "regimes" in EQSAM
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_ZSR_aerosol_water_t type and associated functions. 
module pmc_rxn_ZSR_aerosol_water

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

#define _NUM_PHASE_ this%condensed_data_int(1)
#define _GAS_WATER_ID_ this%condensed_data_int(2)
#define _AERO_WATER_ID_ this%condensed_data_int(3)
#define _NUM_ION_PAIR_ this%condensed_data_int(4)
#define _TOTAL_INT_PARAM_ this%condensed_data_int(5)
#define _TOTAL_FLOAT_PARAM_ this%condensed_data_int(6)
#define _NUM_INT_PROP_ 6
#define _NUM_REAL_PROP_ 0
#define _PHASE_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _PAIR_INT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+x)
#define _PAIR_FLOAT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+_NUM_ION_PAIR_+x)
#define _TYPE_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x))
#define _JACOB_NUM_CATION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+1)
#define _JACOB_NUM_ANION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+2)
#define _JACOB_CATION_ID_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+3)
#define _JACOB_ANION_ID_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+4)
#define _JACOB_NUM_Y_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+5)
#define _EQSAM_NUM_ION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x))
#define _EQSAM_ION_ID_(x,y) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+1+y)
#define _JACOB_low_RH_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x))
#define _JACOB_Y_(x,y) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+1+y)
#define _EQSAM_NW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x))
#define _EQSAM_ZW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+1)
#define _EQSAM_MW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+2)

  public :: rxn_ZSR_aerosol_water_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_ZSR_aerosol_water_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_ZSR_aerosol_water_t

  !> Constructor for rxn_ZSR_aerosol_water_t
  interface rxn_ZSR_aerosol_water_t
    procedure :: constructor
  end interface rxn_ZSR_aerosol_water_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for ZSR aerosol water reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_ZSR_aerosol_water_t), pointer :: new_obj

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
    class(rxn_ZSR_aerosol_water_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)


  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_ZSR_aerosol_water_t), intent(in) :: this
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
    class(rxn_ZSR_aerosol_water_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_PHASE_
#undef _GAS_WATER_ID_
#undef _AERO_WATER_ID_
#undef _NUM_ION_PAIR_
#undef _TOTAL_INT_PARAM_
#undef _TOTAL_FLOAT_PARAM_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_ID_
#undef _PAIR_INT_PARAM_LOC_
#undef _PAIR_FLOAT_PARAM_LOC_
#undef _TYPE_
#undef _JACOB_NUM_CATION_
#undef _JACOB_NUM_ANION_
#undef _JACOB_CATION_ID_
#undef _JACOB_ANION_ID_
#undef _JACOB_NUM_Y_
#undef _EQSAM_NUM_ION_
#undef _EQSAM_ION_ID_
#undef _JACOB_low_RH_
#undef _JACOB_Y_
#undef _EQSAM_NW_
#undef _EQSAM_ZW_
#undef _EQSAM_MW_

end module pmc_rxn_ZSR_aerosol_water
