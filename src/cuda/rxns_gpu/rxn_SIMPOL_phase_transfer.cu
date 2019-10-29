/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Phase Transfer reaction solver functions
 *
*/
/** \file
 * \brief Phase Transfer reaction solver functions
*/
extern "C"{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_rep_solver_gpu.h"
#include "../aeros_gpu/sub_model_solver_gpu.h"

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]
 //todo fix
// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472
// Small number for ignoring low concentrations
#define VERY_SMALL_NUMBER_ 1.0e-30

// Jacobian set indices
#define JAC_GAS 0
#define JAC_AERO 1

#define DELTA_H_ float_data[0*n_rxn]
#define DELTA_S_ float_data[1*n_rxn]
#define DIFF_COEFF_ float_data[2*n_rxn]
#define PRE_C_AVG_ float_data[3*n_rxn]
#define B1_ float_data[4*n_rxn]
#define B2_ float_data[5*n_rxn]
#define B3_ float_data[6*n_rxn]
#define B4_ float_data[7*n_rxn]
#define CONV_ float_data[8*n_rxn]
#define MW_ float_data[9*n_rxn]
#define SMALL_NUMBER_ float_data[10*n_rxn]
#define NUM_AERO_PHASE_ int_data[0*n_rxn]
#define GAS_SPEC_ (int_data[1*n_rxn]-1)
#define C_AVG_ALPHA_ rate_constants[0]
#define EQUIL_CONST_ rate_constants[1]
#define UGM3_TO_PPM_ rate_constants[2]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 11
#define NUM_ENV_PARAM_ 3
#define AERO_SPEC_(x) (int_data[(NUM_INT_PROP_ + x)*n_rxn]-1)
#define AERO_ACT_ID_(x) (int_data[(NUM_INT_PROP_ + NUM_AERO_PHASE_ + x)*n_rxn])
#define AERO_PHASE_ID_(x) (int_data[(NUM_INT_PROP_ + 2*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define AERO_REP_ID_(x) (int_data[(NUM_INT_PROP_ + 3*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define DERIV_ID_(x) (int_data[(NUM_INT_PROP_ + 4*(NUM_AERO_PHASE_) + x)*n_rxn])
#define GAS_ACT_JAC_ID_(x) int_data[(NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x)*n_rxn]
#define AERO_ACT_JAC_ID_(x) int_data[(NUM_INT_PROP_ + 1 + 6*(NUM_AERO_PHASE_) + x)*n_rxn]
#define JAC_ID_(x) (int_data[(NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x)*n_rxn])
#define PHASE_INT_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 10*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define PHASE_FLOAT_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 11*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[PHASE_INT_LOC_(x)*n_rxn])
#define PHASE_JAC_ID_(x,s,e) int_data[(PHASE_INT_LOC_(x)+1+s*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define EFF_RAD_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+e]
#define NUM_CONC_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define MASS_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+2*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define MW_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+3*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]


#ifdef PMC_USE_GPU
__device__
#endif
double rxn_gpu_SIMPOL_phase_transfer_calc_overall_rate(int *rxn_data,
    double *double_pointer_gpu, double *rate_constants, double *state,
    double cond_rc, double evap_rc, int i_phase, int n_rxn2)
{
#ifndef FORCE_CPU
  int n_rxn=n_rxn2;
#else
  int n_rxn=1;
#endif

  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the overall rate.
  // These equations are set up to try to avoid loss of accuracy from
  // subtracting two almost-equal numbers when rate_cond ~ rate_evap.
  // When modifying these calculations, be sure to use the Jacobian checker
  // during unit testing.
  long double rate      = ZERO;
  long double aero_conc = state[AERO_SPEC_(i_phase)];
  long double gas_conc  = state[GAS_SPEC_];
  long double l_cond_rc  = cond_rc;
  long double l_evap_rc  = evap_rc;
  if (l_evap_rc == ZERO || l_cond_rc == ZERO) {
    rate = l_evap_rc * aero_conc - l_cond_rc * gas_conc;
  } else if (l_evap_rc * aero_conc < l_cond_rc * gas_conc) {
    long double gas_eq = aero_conc * ( l_evap_rc / l_cond_rc );
    rate = ( gas_eq - gas_conc ) * l_cond_rc;
  } else {
    long double aero_eq = gas_conc * (l_cond_rc / l_evap_rc);
    rate = (aero_conc - aero_eq) * l_evap_rc;
  }

  return (double) rate;
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
#ifdef PMC_USE_GPU
__device__
#endif
void rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(ModelData *model_data, double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
#ifndef FORCE_CPU
  int n_rxn=n_rxn2;
#else
  int n_rxn=1;
#endif

  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  /*

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    double radius;
    aero_rep_gpu_get_effective_radius(
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
    double number_conc = ONE;
    if (aero_conc_type == 0) {
      aero_rep_gpu_get_number_conc(
  		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc, 		// particle number conc (#/cc)
                  NULL);                        // partial derivative
    }

    // Get the total mass of the aerosol phase
    double aero_phase_mass;
    aero_rep_gpu_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),        // aerosol representation index
                  AERO_PHASE_ID_(i_phase),      // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  NULL);                        // partial derivatives

    // Get the total mass of the aerosol phase
    double aero_phase_avg_MW;
    aero_rep_gpu_get_aero_phase_avg_MW(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),        // aerosol representation index
                  AERO_PHASE_ID_(i_phase),      // aerosol phase index
                  &aero_phase_avg_MW,           // avg MW in the aerosol phase
                  NULL);                        // partial derivatives

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_mass <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    double cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    double evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW /
              aero_phase_mass);

    // Get the activity coefficient (if one exists)
    double act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = state[AERO_ACT_ID_(i_phase)];
    }

    // Calculate aerosol-phase evaporation rate (ppm/s)
    evap_rate *= act_coeff;

    // Calculate the overall rate
    double rate = rxn_gpu_SIMPOL_phase_transfer_calc_overall_rate(
                        (int*) rxn_data, double_pointer_gpu, rate_constants,
                        state, cond_rate, evap_rate, i_phase, n_rxn2);

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    #ifndef FORCE_CPU
      if (DERIV_ID_(0)>=0) deriv[DERIV_ID_(0)] += number_conc * rate;
    #else
      atomicAdd((double*)&(deriv[DERIV_ID_(0)]), number_conc * rate);
    #endif

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (DERIV_ID_(1+i_phase)>=0) {
      #ifndef FORCE_CPU
        deriv[DERIV_ID_(1+i_phase)] -= rate / UGM3_TO_PPM_;
      #else
      atomicAdd((double*)&(deriv[DERIV_ID_(1+i_phase)]),-rate / UGM3_TO_PPM_);
      #endif
    }
  }

*/

  /*
  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Skip reactions that are being treated as instantaneous
    if (FAST_FLUX_(i_phase) != 0.0) continue;

    // Get the particle effective radius (m)
    double radius;
    aero_rep_gpu_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    double number_conc;
    aero_rep_gpu_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number conc (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_gpu_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    double aero_phase_gpu_mass;
    double aero_phase_gpu_avg_MW;
    aero_rep_gpu_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),        // aerosol representation index
                  AERO_PHASE_ID_(i_phase),      // aerosol phase index
                  &aero_phase_gpu_mass,             // total aerosol-phase mass
                  &aero_phase_gpu_avg_MW);          // avg MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_gpu_mass <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    double cond_rate = number_conc/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALHPA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    double evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_gpu_avg_MW /
              aero_phase_gpu_mass);

    // Slow down condensation rate as gas-phase concentrations become small
    double gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    double cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    cond_scaling *= cond_scaling;

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[GAS_SPEC_] * cond_scaling;

    // Get the activity coefficient (if one exists)
    double act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_gpu_get_parameter_value(model_data,
                AERO_ACT_ID_(i_phase));
    }

    // Slow down evaporation as aerosol-phase activity becomes small
    double aero_adj = state[AERO_SPEC_(i_phase)] * act_coeff -
                        VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    double evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    evap_scaling *= evap_scaling;

    // Calculate aerosol-phase evaporation rate (ppm/s)
    // (Slow down evaporation as aerosol-phase concentrations approach zero
    //  to help out the solver.)
    evap_rate *= state[AERO_SPEC_(i_phase)] * act_coeff * evap_scaling;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    //if (DERIV_ID_(0)>=0) deriv[DERIV_ID_(0)] += evap_rate - cond_rate;
    if (DERIV_ID_(0)>=0) atomicAdd((double*)&(deriv[DERIV_ID_(0)]),evap_rate - cond_rate);

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (DERIV_ID_(1+i_phase)>=0) {
      if (aero_conc_type==0) {
        // Per-particle condensation
        //deriv[DERIV_ID_(1+i_phase)] += (cond_rate - evap_rate) /
        //        UGM3_TO_PPM_ / number_conc;
        atomicAdd((double*)&(deriv[DERIV_ID_(1+i_phase)]),(cond_rate - evap_rate) /
          UGM3_TO_PPM_ / number_conc);

      } else {
        // Total aerosol mass condensation
        //deriv[DERIV_ID_(1+i_phase)] += (cond_rate - evap_rate) /
        //        UGM3_TO_PPM_;
        atomicAdd((double*)&(deriv[DERIV_ID_(1+i_phase)]),(cond_rate - evap_rate) /
          UGM3_TO_PPM_);
      }
    }
  }
*/

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

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */

#ifdef PMC_USE_SUNDIALS
#ifdef PMC_USE_GPU
__device__
#endif
void rxn_gpu_SIMPOL_phase_transfer_calc_jac_contrib(ModelData *model_data, double *rate_constants, double *state,
          double *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  /*
  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Skip reactions that are being treated as instantaneous
    if (FAST_FLUX_(i_phase) != 0.0) continue;

    // Get the particle effective radius (m)
    double radius;
    aero_rep_gpu_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    double number_conc;
    aero_rep_gpu_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number conc (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_gpu_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    double aero_phase_gpu_mass;
    double aero_phase_gpu_avg_MW;
    aero_rep_gpu_get_aero_phase_mass(
                  model_data,                   // model data
                  AERO_REP_ID_(i_phase),       // aerosol representation index
                  AERO_PHASE_ID_(i_phase),     // aerosol phase index
                  &aero_phase_gpu_mass,             // total aerosol-phase mass
                  &aero_phase_gpu_avg_MW);          // avg MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO || aero_phase_gpu_mass <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    double cond_rate = number_conc/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALHPA_));

    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    double evap_rate = cond_rate * (EQUIL_CONST_ * aero_phase_gpu_avg_MW /
              aero_phase_gpu_mass);

    // Slow down condensation rate as gas-phase concentrations become small
    double gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    double cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    double cond_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  gas_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -gas_adj / SMALL_NUMBER_ ) ) );
    cond_scaling_deriv *= 2.0 * cond_scaling;
    cond_scaling *= cond_scaling;

    // Get the activity coefficient (if one exists)
    double act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_gpu_get_parameter_value(model_data,
                AERO_ACT_ID_(i_phase));
    }

    // Slow down evaporation as aerosol-phase activity becomes small
    double aero_adj = state[AERO_SPEC_(i_phase)] * act_coeff -
                        VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    double evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    double evap_scaling_deriv =
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
*/

}
#endif

}