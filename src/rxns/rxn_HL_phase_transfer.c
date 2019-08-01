/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Phase Transfer reaction solver functions
 *
*/
/** \file
 * \brief Phase Transfer reaction solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aero_rep_solver.h"
#include "../rxns.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472
// Small number for ignoring low concentrations
#define VERY_SMALL_NUMBER_ 1.0e-30
// Factor used to calculate minimum aerosol water concentrations for
// HL phase transfer
#define MIN_WATER_ 1.0e-4

// Jacobian set indices
#define JAC_GAS 0
#define JAC_AERO 1

#define DELTA_H_ float_data[0]
#define DELTA_S_ float_data[1]
#define DIFF_COEFF_ float_data[2]
#define PRE_C_AVG_ float_data[3]
#define A_ float_data[4]
#define C_ float_data[5]
#define C_AVG_ALPHA_ float_data[6]
#define EQUIL_CONST_ float_data[7]
#define CONV_ float_data[8]
#define MW_ float_data[9]
#define UGM3_TO_PPM_ float_data[10]
#define SMALL_NUMBER_ float_data[11]
#define NUM_AERO_PHASE_ int_data[0]
#define GAS_SPEC_ (int_data[1]-1)
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 12
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 1 + NUM_AERO_PHASE_ + x]
#define PHASE_INT_LOC_(x) (int_data[NUM_INT_PROP_ + 2 + 6*NUM_AERO_PHASE_ + x]-1)
#define PHASE_REAL_LOC_(x) (int_data[NUM_INT_PROP_ + 2 + 7*NUM_AERO_PHASE_ + x]-1)
#define AERO_SPEC_(x) (int_data[PHASE_INT_LOC_(x)]-1)
#define AERO_WATER_(x) (int_data[PHASE_INT_LOC_(x) + 1]-1)
#define AERO_PHASE_ID_(x) (int_data[PHASE_INT_LOC_(x) + 2]-1)
#define AERO_REP_ID_(x) (int_data[PHASE_INT_LOC_(x) + 3]-1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[PHASE_INT_LOC_(x) + 4])
#define PHASE_JAC_ID_(x,s,e) int_data[PHASE_INT_LOC_(x) + 5 + s*NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define SMALL_WATER_CONC_(x) (float_data[PHASE_REAL_LOC_(x)])
#define EFF_RAD_JAC_ELEM_(x,e) float_data[PHASE_REAL_LOC_(x) + 1 + e]
#define NUM_CONC_JAC_ELEM_(x,e) float_data[PHASE_REAL_LOC_(x) + 1 + NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define INT_DATA_SIZE_ (PHASE_INT_LOC_(NUM_AERO_PHASE_-1)+5+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))
#define FLOAT_DATA_SIZE_ (PHASE_REAL_LOC_(NUM_AERO_PHASE_-1)+1+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param model_data Pointer to the model data
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_get_used_jac_elem(ModelData *model_data,
          void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  bool *aero_jac_elem = (bool*) malloc(sizeof(bool) * model_data->n_state_var);
  if (aero_jac_elem==NULL) {
    printf("\n\nERROR allocating space for 1D jacobian structure array for HL "
           "partitioning reaction\n\n");
    exit(1);
  }

  jac_struct[GAS_SPEC_][GAS_SPEC_] = true;
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    jac_struct[AERO_SPEC_(i_aero_phase)][GAS_SPEC_] = true;
    jac_struct[GAS_SPEC_][AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[AERO_SPEC_(i_aero_phase)][AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[GAS_SPEC_][AERO_WATER_(i_aero_phase)] = true;
    jac_struct[AERO_SPEC_(i_aero_phase)][AERO_WATER_(i_aero_phase)] = true;

    for (int i_elem = 0; i_elem < model_data->n_state_var; i_elem++)
      aero_jac_elem[i_elem] = false;

    int n_jac_elem = aero_rep_get_used_jac_elem( model_data,
                                                 AERO_REP_ID_(i_aero_phase),
                                                 AERO_PHASE_ID_(i_aero_phase),
                                                 aero_jac_elem );
    int i_used_elem = 0;
    for (int i_elem = 0; i_elem < model_data->n_state_var; i_elem++) {
      if (aero_jac_elem[i_elem] == true) {
        jac_struct[GAS_SPEC_][i_elem] = true;
        jac_struct[AERO_SPEC_(i_aero_phase)][i_elem] = true;
        PHASE_JAC_ID_(i_aero_phase,JAC_GAS,i_used_elem) = i_elem;
        PHASE_JAC_ID_(i_aero_phase,JAC_AERO,i_used_elem) = i_elem;
        i_used_elem++;
      }
    }
    for (; i_used_elem<NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         ++i_used_elem) {
      PHASE_JAC_ID_(i_aero_phase,JAC_GAS,i_used_elem) = -1;
      PHASE_JAC_ID_(i_aero_phase,JAC_AERO,i_used_elem) = -1;
    }
    if (i_used_elem != n_jac_elem) {
      printf("\n\nERROR Error setting used Jacobian elements in HL "
             "partitioning reaction %d %d\n\n", i_used_elem, n_jac_elem);
      rxn_HL_phase_transfer_print(rxn_data);
      exit(1);
    }

  }

  free(aero_jac_elem);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Update the time derivative ids
  DERIV_ID_(0) = deriv_ids[GAS_SPEC_];
  for (int i=0; i < NUM_AERO_PHASE_; i++)
	  DERIV_ID_(i + 1) = deriv_ids[AERO_SPEC_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][GAS_SPEC_];
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    JAC_ID_(i_jac++) = jac_ids[AERO_SPEC_(i_aero_phase)][GAS_SPEC_];
    JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][AERO_SPEC_(i_aero_phase)];
    JAC_ID_(i_jac++) =
            jac_ids[AERO_SPEC_(i_aero_phase)][AERO_SPEC_(i_aero_phase)];
    JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][AERO_WATER_(i_aero_phase)];
    JAC_ID_(i_jac++) =
            jac_ids[AERO_SPEC_(i_aero_phase)][AERO_WATER_(i_aero_phase)];
    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         i_elem++) {
      if (PHASE_JAC_ID_(i_aero_phase,JAC_GAS,i_elem) > 0) {
        PHASE_JAC_ID_(i_aero_phase,JAC_GAS,i_elem) =
              jac_ids[GAS_SPEC_][PHASE_JAC_ID_(i_aero_phase,JAC_GAS,i_elem)];
      }
      if (PHASE_JAC_ID_(i_aero_phase,JAC_AERO,i_elem) > 0) {
        PHASE_JAC_ID_(i_aero_phase,JAC_AERO,i_elem) =
              jac_ids[AERO_SPEC_(i_aero_phase)]
                     [PHASE_JAC_ID_(i_aero_phase,JAC_AERO,i_elem)];
      }
    }
  }

  // Calculate a small number based on the integration tolerances to use
  // during solving. TODO find a better place to do this
  realtype *abs_tol = model_data->abs_tol;
  SMALL_NUMBER_ = ( abs_tol[GAS_SPEC_] > abs_tol[AERO_SPEC_(0)] ?
                    abs_tol[AERO_SPEC_(0)] / 10.0 : abs_tol[GAS_SPEC_] / 10.0 );

  // Calculate a small concentration for aerosol-phase water based on the
  // integration tolerances to use during solving. TODO find a better place
  // to do this
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    SMALL_WATER_CONC_(i_aero_phase) = abs_tol[AERO_WATER_(i_aero_phase)] / 10.0;
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Phase Transfer reaction this only involves recalculating the rate
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_update_env_state(double *env_data,
          void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the mass accomodation coefficient if the N* parameter
  // was provided, otherwise set it to 1.0
  double mass_acc = 1.0;
  if (DELTA_H_!=0.0 || DELTA_S_!=0.0) {
    double del_G = DELTA_H_ - TEMPERATURE_K_ * DELTA_S_;
    mass_acc = exp(-del_G/(UNIV_GAS_CONST_ * TEMPERATURE_K_));
    mass_acc = mass_acc / (1.0 + mass_acc);
  }

  // Save c_rms * mass_acc for use in mass transfer rate calc
  C_AVG_ALPHA_ = PRE_C_AVG_ * sqrt(TEMPERATURE_K_) * mass_acc;

  // Calculate the Henry's Law equilibrium rate constant in units of
  // (ug_x/ug_H2O/ppm) where x is the aerosol-phase species. (A was saved in
  // units of M/ppm.)
  if (C_==0.0) {
    EQUIL_CONST_ = A_ * MW_;
  } else {
    EQUIL_CONST_ = A_ * exp(C_ * (1.0/TEMPERATURE_K_ - 1.0/298.0)) * MW_;
  }

  // Calculate the conversion from ug/m^3 -> ppm
  UGM3_TO_PPM_ = CONV_ * TEMPERATURE_K_ / PRESSURE_PA_;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}

/** \brief Calculate the overall per-particle reaction rate ([ppm]/s)
 *
 * This function is called by the deriv and Jac functions to get the overall
 * reaction rate on a per-particle basis, trying to avoid floating-point
 * subtraction errors.
 *
 * \param rxn_data Pointer to the reaction data
 * \param state State array
 * \param cond_rc Condensation rate constant ([ppm]^(-n_react+1)s^-1)
 * \param evap_rc Evaporation rate constant ([ppm]^(-n_prod+1)s^-1)
 * \param i_phase Index for the aerosol phase being calculated
 * \return The calculated overall rate ([ppm]/s)
 */
#ifdef PMC_USE_SUNDIALS
realtype calculate_overall_rate(void *rxn_data, realtype *state,
          realtype cond_rc, realtype evap_rc, int i_phase)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the overall rate.
  // These equations are set up to try to avoid loss of accuracy from
  // subtracting two almost-equal numbers when rate_cond ~ rate_evap.
  // When modifying these calculations, be sure to use the Jacobian checker
  // during unit testing.
  long double rate       = ZERO;
  long double aero_conc  = state[AERO_SPEC_(i_phase)];
  long double aero_water = state[AERO_WATER_(i_phase)];
  long double gas_conc   = state[GAS_SPEC_];
  long double l_cond_rc  = cond_rc;
  long double l_evap_rc  = evap_rc;
  if (l_evap_rc== ZERO || l_cond_rc == ZERO) {
    rate = l_evap_rc * aero_conc / aero_water - l_cond_rc * gas_conc;
  } else if (l_evap_rc * aero_conc / aero_water <
             l_cond_rc * gas_conc) {
    realtype gas_eq = aero_conc * ( l_evap_rc / l_cond_rc );
    rate = ( gas_eq - gas_conc * aero_water) * ( l_cond_rc / aero_water );
  } else {
    realtype aero_eq = gas_conc * aero_water * ( l_cond_rc / l_evap_rc );
    rate = ( aero_conc - aero_eq ) * ( l_evap_rc / aero_water );
  }
  return (long double) rate;
}
#endif

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being computed (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void * rxn_HL_phase_transfer_calc_deriv_contrib(ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius,                      // particle effective radius (m)
                  NULL);                        // partial derivative

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the particle number concentration (#/cc) for per-particle mass
    // concentrations; otherwise set to 1
    realtype number_conc = ONE;
    if (aero_conc_type == 0) {
      aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc, 		// particle number concentration
                                                // (#/cc)
                  NULL);                        // partial derivative
    }

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] * number_conc <
        MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) continue;

    // Slow down rates as water approaches the minimum value
#if 0
    realtype water_adj = state[AERO_WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    realtype water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    water_scaling *= water_scaling;
#endif
    realtype water_scaling = ONE;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (EQUIL_CONST_);

    // Slow down condensation rate as gas-phase concentrations become small
#if 0
    realtype gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    realtype cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    cond_scaling *= cond_scaling;
#endif
    realtype cond_scaling = ONE;

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= cond_scaling * water_scaling;

    // Slow down evaporation as aerosol-phase concentrations become small
#if 0
    realtype aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    evap_scaling *= evap_scaling;
#endif
    realtype evap_scaling = ONE;

    // Calculate aerosol-phase evaporation rate (ug/m^3/s)
    evap_rate *= evap_scaling * water_scaling;

    // Get the overall rate
    realtype rate = calculate_overall_rate(rxn_data, state, cond_rate,
                                           evap_rate, i_phase);

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (DERIV_ID_(0)>=0) {
      deriv[DERIV_ID_(0)] += number_conc * rate;
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (DERIV_ID_(1+i_phase)>=0)
      deriv[DERIV_ID_(1+i_phase)] -= rate / UGM3_TO_PPM_;

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void * rxn_HL_phase_transfer_calc_jac_contrib(ModelData *model_data,
          realtype *J, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			     // model data
		  AERO_REP_ID_(i_phase),	     // aerosol representation index
		  AERO_PHASE_ID_(i_phase),	     // aerosol phase index
		  &radius,                           // particle effective radius (m)
                  &(EFF_RAD_JAC_ELEM_(i_phase,0)));  // partial derivative

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			     // model data
		  AERO_REP_ID_(i_phase),	     // aerosol representation index
		  AERO_PHASE_ID_(i_phase));	     // aerosol phase index

    // Get the particle number concentration (#/cc) for per-particle
    // concentrations
    realtype number_conc = ONE;
    if (aero_conc_type == 0) {
    aero_rep_get_number_conc(
		  model_data,			     // model data
		  AERO_REP_ID_(i_phase),	     // aerosol representation index
		  AERO_PHASE_ID_(i_phase),	     // aerosol phase index
		  &number_conc, 		     // particle number concentration
                                                     // (#/cc)
                  &(NUM_CONC_JAC_ELEM_(i_phase,0))); // partial derivative
    } else {
      for (int i_elem=0; i_elem<NUM_AERO_PHASE_JAC_ELEM_(i_phase); ++i_elem)
        NUM_CONC_JAC_ELEM_(i_phase, i_elem) = ZERO;
    }

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] * number_conc <
        MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) continue;

    // Slow down rates as water approaches the minimum value
#if 0
    realtype water_adj = state[AERO_WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    realtype water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    realtype water_scaling_deriv =
      2.0 / ( SMALL_WATER_CONC_(i_phase) *
              ( exp(  water_adj / SMALL_WATER_CONC_(i_phase) ) + 2.0 +
                exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) );
    water_scaling_deriv *= 2.0 * water_scaling;
    water_scaling *= water_scaling;
#endif
    realtype water_scaling = ONE;
    realtype water_scaling_deriv = ZERO;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (EQUIL_CONST_);

    // Slow down condensation rate as gas-phase concentrations become small
#if 0
    realtype gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    realtype cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    realtype cond_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  gas_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -gas_adj / SMALL_NUMBER_ ) ) );
    cond_scaling_deriv *= 2.0 * cond_scaling;
    cond_scaling *= cond_scaling;
#endif
    realtype cond_scaling = ONE;
    realtype cond_scaling_deriv = ZERO;

    // Slow down evaporation as aerosol-phase concentrations become small
#if 0
    realtype aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    realtype evap_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  aero_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -aero_adj / SMALL_NUMBER_ ) ) );
    evap_scaling_deriv *= 2.0 * evap_scaling;
    evap_scaling *= evap_scaling;
#endif
    realtype evap_scaling = ONE;
    realtype evap_scaling_deriv = ZERO;

    // Get the overall rate for certain Jac elements
    realtype rate = calculate_overall_rate(rxn_data, state,
        cond_rate * cond_scaling * water_scaling,
        evap_rate * evap_scaling * water_scaling, i_phase);

    // Update evap rate to be for aerosol species concentrations
    evap_rate /= (UGM3_TO_PPM_ * state[AERO_WATER_(i_phase)]);

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (JAC_ID_(1+i_phase*5+1)>=0)
              J[JAC_ID_(1+i_phase*5+1)] += number_conc * evap_rate *
                    water_scaling * UGM3_TO_PPM_ *
                    ( evap_scaling + state[AERO_SPEC_(i_phase)] *
                                     evap_scaling_deriv );
    if (JAC_ID_(1+i_phase*5+3)>=0)
	      J[JAC_ID_(1+i_phase*5+3)] += number_conc * evap_rate *
                    evap_scaling * UGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)] *
                    ( - water_scaling / state[AERO_WATER_(i_phase)] +
                      water_scaling_deriv );
    if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= number_conc * cond_rate *
                      water_scaling *
                      ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (JAC_ID_(1+i_phase*5)>=0) J[JAC_ID_(1+i_phase*5)] +=
            cond_rate * water_scaling / UGM3_TO_PPM_ *
            ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
    if (JAC_ID_(1+i_phase*5+2)>=0) J[JAC_ID_(1+i_phase*5+2)] -= evap_rate *
            water_scaling *
            ( evap_scaling + state[AERO_SPEC_(i_phase)] * evap_scaling_deriv );
    if (JAC_ID_(1+i_phase*5+4)>=0) J[JAC_ID_(1+i_phase*5+4)] -= evap_rate *
	    evap_scaling * state[AERO_SPEC_(i_phase)] *
            ( - water_scaling / state[AERO_WATER_(i_phase)] +
              water_scaling_deriv );

    // Add contributions from species used in aerosol property calculations

    // Calculate d_rate/d_effecive_radius and d_rate/d_number_concentration
    realtype d_rate_d_radius = rate *
        -(2.0*radius/(3.0*DIFF_COEFF_) + 4.0/(3.0*C_AVG_ALPHA_)) /
        (2.0*radius*radius/(3.0*DIFF_COEFF_) + 4.0*radius/(3.0*C_AVG_ALPHA_));
    realtype d_rate_d_number = rate / number_conc;

    // Loop through Jac elements and update
    for (int i_elem=0; i_elem<NUM_AERO_PHASE_JAC_ELEM_(i_phase); ++i_elem) {

      // Gas-phase species dependencies
      if (PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem) > 0) {

        // species involved in effective radius calculation
        J[PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem)] += number_conc *
            d_rate_d_radius * EFF_RAD_JAC_ELEM_(i_phase, i_elem);

        // species involved in numer concentration
        J[PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem)] += number_conc *
            d_rate_d_number * NUM_CONC_JAC_ELEM_(i_phase, i_elem);

      }

      // Aerosol-phase species dependencies
      if (PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem) > 0) {

        // species involved in effective radius calculation
        J[PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem)] -= d_rate_d_radius /
            UGM3_TO_PPM_ * EFF_RAD_JAC_ELEM_(i_phase, i_elem);

        // species involved in numer concentration
        J[PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem)] -= d_rate_d_number /
            UGM3_TO_PPM_ * NUM_CONC_JAC_ELEM_(i_phase, i_elem);

      }

    }

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nPhase Transfer reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}
