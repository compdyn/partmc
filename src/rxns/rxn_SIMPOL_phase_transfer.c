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
#include "../sub_model_solver.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472
// Small number for ignoring low concentrations
#define VERY_SMALL_NUMBER_ 1.0e-30

#define DELTA_H_ float_data[0]
#define DELTA_S_ float_data[1]
#define DIFF_COEFF_ float_data[2]
#define PRE_C_AVG_ float_data[3]
#define B1_ float_data[4]
#define B2_ float_data[5]
#define B3_ float_data[6]
#define B4_ float_data[7]
#define C_AVG_ALHPA_ float_data[8]
#define EQUIL_CONST_ float_data[9]
#define CONV_ float_data[10]
#define MW_ float_data[11]
#define UGM3_TO_PPM_ float_data[12]
#define SMALL_NUMBER_ float_data[13]
#define NUM_AERO_PHASE_ int_data[0]
#define GAS_SPEC_ (int_data[1]-1)
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 14
#define AERO_SPEC_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define AERO_ACT_ID_(x) (int_data[NUM_INT_PROP_ + NUM_AERO_PHASE_ + x])
#define AERO_PHASE_ID_(x) (int_data[NUM_INT_PROP_ + 2*(NUM_AERO_PHASE_) + x]-1)
#define AERO_REP_ID_(x) (int_data[NUM_INT_PROP_ + 3*(NUM_AERO_PHASE_) + x]-1)
#define DERIV_ID_(x) (int_data[NUM_INT_PROP_ + 4*(NUM_AERO_PHASE_) + x])
#define JAC_ID_(x) (int_data[NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x])
#define FAST_FLUX_(x) (float_data[NUM_FLOAT_PROP_ + x])
#define AERO_ADJ_(x) (float_data[NUM_FLOAT_PROP_ + NUM_AERO_PHASE_ + x])
#define INT_DATA_SIZE_ (NUM_INT_PROP_+2+8*(NUM_AERO_PHASE_))
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+2*(NUM_AERO_PHASE_))

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_get_used_jac_elem(void *rxn_data,
          bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  jac_struct[GAS_SPEC_][GAS_SPEC_] = true;
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    jac_struct[AERO_SPEC_(i_aero_phase)][GAS_SPEC_] = true;
    jac_struct[GAS_SPEC_][AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[AERO_SPEC_(i_aero_phase)][AERO_SPEC_(i_aero_phase)] = true;
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data for finding sub model ids
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_update_ids(ModelData *model_data,
          int *deriv_ids, int **jac_ids, void *rxn_data)
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
    }

  // Find activity coefficient ids, if they exist
  // FIXME Don't hard-code sub model ids
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    int aero_state_id = AERO_SPEC_(i_aero_phase);
    AERO_ACT_ID_(i_aero_phase) =
      sub_model_get_parameter_id(model_data, 1, (void*) (&aero_state_id));
  }

  // Calculate a small number based on the integration tolerances to use
  // during solving. TODO find a better place to do this
  realtype *abs_tol = model_data->abs_tol;
  SMALL_NUMBER_ = ( abs_tol[GAS_SPEC_] > abs_tol[AERO_SPEC_(0)] ?
                    abs_tol[AERO_SPEC_(0)] / 10.0 : abs_tol[GAS_SPEC_] / 10.0 );

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
void * rxn_SIMPOL_phase_transfer_update_env_state(double *env_data,
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
  C_AVG_ALHPA_ = PRE_C_AVG_ * sqrt(TEMPERATURE_K_) * mass_acc;

  // SIMPOL.1 vapor pressure (Pa)
  double vp = B1_ / TEMPERATURE_K_
                + B2_ + B3_ * TEMPERATURE_K_
                + B4_ * log(TEMPERATURE_K_);
  vp = 101325.0 * pow(10, vp);

  // Calculate the conversion from ug_x/m^3 -> ppm_x
  UGM3_TO_PPM_ = CONV_ * TEMPERATURE_K_ / PRESSURE_PA_;

  // Calculate the partitioning coefficient K_eq (ppm_x/ug_x*ug_tot/kg_tot)
  // such that for partitioning species X at equilibrium:
  //   [X]_gas = [X]_aero * activity_coeff_X * K_eq * MW_tot_aero / [tot]_aero
  // where 'tot' indicates all species within an aerosol phase combined
  // with []_gas in (ppm) and []_aero in (ug/m^3)
  EQUIL_CONST_ = vp                    // (Pa_x*mol_tot/mol_x)
                  / PRESSURE_PA_       // (1/Pa_air)
                  / MW_                // (mol_x/kg_x)
                  * 1.0e6;             // 1.0e6ppm_x*Pa_air/Pa_x *
                                       //  1.0e-9kg_x/ug_x * 1.0e9ug_tot/kg_tot

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Check whether phase transfer is fast enough at this timestep to treat
 * it as essentially going to completion
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current solver time step (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_pre_calc(ModelData *model_data, void *rxn_data,
          double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Get the total mass available for transfer (ppm)
  realtype total_mass = state[GAS_SPEC_];

  // Calculate the net equilibrium constant across all phases
  realtype net_eq_const = ZERO;

  // Calculate fast phase-transfer rates for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number conc (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    realtype aero_phase_mass;
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),        // aerosol representation index
                  AERO_PHASE_ID_(i_phase),      // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  &aero_phase_avg_MW);          // avg MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_mass <= ZERO) {
      FAST_FLUX_(i_phase) = ZERO;
      AERO_ADJ_(i_phase)  = ZERO;
      continue;
    }

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = number_conc/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALHPA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    realtype evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW /
              aero_phase_mass);

    // Determine if transfer is fast enough to treat as being at equilibrium
    if ( ONE / cond_rate > time_step / 1000.0 &&
         ONE / evap_rate > time_step / 1000.0 ) {
      FAST_FLUX_(i_phase) = ZERO;
      AERO_ADJ_(i_phase)  = ZERO;
      continue;
    }

    // Add the mass in this phase to the total mass available for transfer
    total_mass += state[AERO_SPEC_(i_phase)] * UGM3_TO_PPM_;

    // Temporarily store the equilibrium constant for this phase
    AERO_ADJ_(i_phase) = cond_rate / evap_rate;

    // Add the equilibrium constant for this phase to the net value
    net_eq_const += AERO_ADJ_(i_phase);

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[GAS_SPEC_];

    // Get the activity coefficient (if one exists)
    realtype act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_get_parameter_value(model_data,
                AERO_ACT_ID_(i_phase));
    }

    // Calculate aerosol-phase evaporation rate (ppm/s)
    evap_rate *= state[AERO_SPEC_(i_phase)] * act_coeff;

    // Calculate the flux for this aerosol phase instance
    FAST_FLUX_(i_phase) = (cond_rate / UGM3_TO_PPM_ - evap_rate);
    model_data->rel_flux[AERO_SPEC_(i_phase)] += FAST_FLUX_(i_phase);
    model_data->rel_flux[GAS_SPEC_] -= FAST_FLUX_(i_phase) * UGM3_TO_PPM_;

  }

  // Return if no rates are fast enough for equilibrium treatment
  if (total_mass == ZERO) return (void*) &(float_data[FLOAT_DATA_SIZE_]);

  // Calculate the new gas-phase concentration
  realtype gas_conc = total_mass / ( net_eq_const + ONE );
  model_data->state_adj[GAS_SPEC_] = gas_conc - state[GAS_SPEC_];

  // Calculate the new aerosol-phase concentrations
  realtype total_aero_conc = (total_mass - gas_conc) / UGM3_TO_PPM_;
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Skip slow phase transfer
    if (AERO_ADJ_(i_phase)==ZERO) continue;

    // Calculate the new aerosol conc for this phase
    AERO_ADJ_(i_phase) = AERO_ADJ_(i_phase) / net_eq_const * total_aero_conc -
                         state[AERO_SPEC_(i_phase)];
    model_data->state_adj[AERO_SPEC_(i_phase)] += AERO_ADJ_(i_phase);

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Scale state adjustments from this reaction
 *
 * \param model_data Pointer to the model data
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction
 */
void * rxn_SIMPOL_phase_transfer_scale_adj(ModelData *model_data, void *rxn_data)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // FIXME Figure out a better way to scale the adjustments

  // Gas-phase over-depletion (ppm)
  realtype gas_over_dep = -(state[GAS_SPEC_] +
                            model_data->state_adj[GAS_SPEC_]);
  if (gas_over_dep > ZERO) {

    // Scale the adjustment for each aerosol phase
    for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {
      realtype scale_factor = FAST_FLUX_(i_phase) /
                              model_data->rel_flux[GAS_SPEC_];
      model_data->state_adj[GAS_SPEC_] += gas_over_dep * scale_factor;
      model_data->state_adj[AERO_SPEC_(i_phase)] -= gas_over_dep *
                                                    scale_factor /
                                                    UGM3_TO_PPM_;
    }

  } else {

    // Check and scale the adjustment for each aerosol phase
    for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

      realtype aero_over_dep = -(state[AERO_SPEC_(i_phase)] +
                                 model_data->state_adj[AERO_SPEC_(i_phase)]);

      if (aero_over_dep < ZERO) continue;

      realtype scale_factor = FAST_FLUX_(i_phase) /
                              model_data->rel_flux[AERO_SPEC_(i_phase)];
      model_data->state_adj[GAS_SPEC_] -= aero_over_dep * scale_factor *
                                          UGM3_TO_PPM_;
      model_data->state_adj[AERO_SPEC_(i_phase)] += aero_over_dep *
                                                    scale_factor;
    }

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

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
void * rxn_SIMPOL_phase_transfer_calc_deriv_contrib(ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Skip reactions that are being treated as instantaneous
    if (FAST_FLUX_(i_phase) != 0.0) continue;

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number conc (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    realtype aero_phase_mass;
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),        // aerosol representation index
                  AERO_PHASE_ID_(i_phase),      // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  &aero_phase_avg_MW);          // avg MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_mass <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = number_conc/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALHPA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    realtype evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW /
              aero_phase_mass);

    // Slow down condensation rate as gas-phase concentrations become small
    realtype gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    realtype cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    cond_scaling *= cond_scaling;

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[GAS_SPEC_] * cond_scaling;

    // Get the activity coefficient (if one exists)
    realtype act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_get_parameter_value(model_data,
                AERO_ACT_ID_(i_phase));
    }

    // Slow down evaporation as aerosol-phase activity becomes small
    realtype aero_adj = state[AERO_SPEC_(i_phase)] * act_coeff -
                        VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    evap_scaling *= evap_scaling;

    // Calculate aerosol-phase evaporation rate (ppm/s)
    // (Slow down evaporation as aerosol-phase concentrations approach zero
    //  to help out the solver.)
    evap_rate *= state[AERO_SPEC_(i_phase)] * act_coeff * evap_scaling;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (DERIV_ID_(0)>=0) deriv[DERIV_ID_(0)] += evap_rate - cond_rate;

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (DERIV_ID_(1+i_phase)>=0) {
      if (aero_conc_type==0) {
        // Per-particle condensation
        deriv[DERIV_ID_(1+i_phase)] += (cond_rate - evap_rate) /
                UGM3_TO_PPM_ / number_conc;
      } else {
        // Total aerosol mass condensation
        deriv[DERIV_ID_(1+i_phase)] += (cond_rate - evap_rate) /
                UGM3_TO_PPM_;
      }
    }
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
void * rxn_SIMPOL_phase_transfer_calc_jac_contrib(ModelData *model_data,
          realtype *J, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Skip reactions that are being treated as instantaneous
    if (FAST_FLUX_(i_phase) != 0.0) continue;

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number conc (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    realtype aero_phase_mass;
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),       // aerosol representation index
                  AERO_PHASE_ID_(i_phase),     // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  &aero_phase_avg_MW);          // avg MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_mass <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = number_conc/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALHPA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    realtype evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW /
              aero_phase_mass);

    // Slow down condensation rate as gas-phase concentrations become small
    realtype gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    realtype cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    realtype cond_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  gas_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -gas_adj / SMALL_NUMBER_ ) ) );
    cond_scaling_deriv *= 2.0 * cond_scaling;
    cond_scaling *= cond_scaling;

    // Get the activity coefficient (if one exists)
    realtype act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_get_parameter_value(model_data,
                AERO_ACT_ID_(i_phase));
    }

    // Slow down evaporation as aerosol-phase activity becomes small
    realtype aero_adj = state[AERO_SPEC_(i_phase)] * act_coeff -
                        VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    realtype evap_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  aero_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -aero_adj / SMALL_NUMBER_ ) ) );
    evap_scaling_deriv *= 2.0 * evap_scaling;
    evap_scaling *= evap_scaling;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
      if (JAC_ID_(1+i_phase*3+1)>=0)
          J[JAC_ID_(1+i_phase*3+1)] += evap_rate * act_coeff *
                                       ( evap_scaling +
                                         state[AERO_SPEC_(i_phase)] *
                                         evap_scaling_deriv );
      if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= cond_rate *
                                          ( cond_scaling +
                                            state[GAS_SPEC_] *
                                            cond_scaling_deriv );

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (aero_conc_type==0) {
      // Per-particle condensation
      if (JAC_ID_(1+i_phase*3)>=0) J[JAC_ID_(1+i_phase*3)] +=
          cond_rate / number_conc / UGM3_TO_PPM_ *
          ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
      if (JAC_ID_(1+i_phase*3+2)>=0) J[JAC_ID_(1+i_phase*3+2)] -=
          evap_rate * act_coeff / number_conc / UGM3_TO_PPM_ *
          ( evap_scaling + state[AERO_SPEC_(i_phase)] * evap_scaling_deriv );
    } else {
      // Total aerosol mass condensation
      if (JAC_ID_(1+i_phase*3)>=0) J[JAC_ID_(1+i_phase*3)] +=
          cond_rate / UGM3_TO_PPM_ *
          ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
      if (JAC_ID_(1+i_phase*3+2)>=0) J[JAC_ID_(1+i_phase*3+2)] -=
          evap_rate * act_coeff / UGM3_TO_PPM_ *
          ( evap_scaling + state[AERO_SPEC_(i_phase)] * evap_scaling_deriv );
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
void * rxn_SIMPOL_phase_transfer_skip(void *rxn_data)
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
void * rxn_SIMPOL_phase_transfer_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nSIMPOL.1 Phase Transfer reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef UNIV_GAS_CONST_
#undef VERY_SMALL_NUMBER_
#undef RATE_SCALING_

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
#undef SMALL_NUMBER_
#undef NUM_AERO_PHASE_
#undef GAS_SPEC_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef AERO_SPEC_
#undef AERO_ACT_ID_
#undef AERO_PHASE_ID_
#undef AERO_REP_ID_
#undef DERIV_ID_
#undef JAC_ID_
#undef FAST_FLUX_
#undef AERO_ADJ_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
