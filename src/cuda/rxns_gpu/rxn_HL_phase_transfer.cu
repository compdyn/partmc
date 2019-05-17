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
extern "C" {
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_rep_solver_gpu.h"
#include "../rxns_gpu.h"

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
#define PHASE_JAC_ID_(x, s, e) int_data[PHASE_INT_LOC_(x) + 5 + s*NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define SMALL_WATER_CONC_(x) (float_data[PHASE_REAL_LOC_(x)])
#define FAST_FLUX_(x) (float_data[PHASE_REAL_LOC_(x) + 1])
#define AERO_ADJ_(x) (float_data[PHASE_REAL_LOC_(x) + 2])
#define EFF_RAD_JAC_ELEM_(x, e) float_data[PHASE_REAL_LOC_(x) + 3 + e]
#define NUM_CONC_JAC_ELEM_(x, e) float_data[PHASE_REAL_LOC_(x) + 3 + NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define INT_DATA_SIZE_ (PHASE_INT_LOC_(NUM_AERO_PHASE_-1)+5+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))
#define FLOAT_DATA_SIZE_ (PHASE_REAL_LOC_(NUM_AERO_PHASE_-1)+3+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))

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
__device__ void rxn_gpu_HL_phase_transfer_calc_deriv_contrib(ModelDatagpu *model_data,
          realtype *deriv, void *rxn_data, double time_step)
  //void rxn_gpu_HL_phase_transfer_calc_deriv_contrib(ModelDatagpu *model_data,
  //        realtype *deriv, void *rxn_data, double time_step)
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
    aero_rep_gpu_get_effective_radius(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase),	// aerosol phase index
          &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_gpu_get_number_conc(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase),	// aerosol phase index
          &number_conc);		// particle number concentration
                                                // (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_gpu_get_aero_conc_type(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] *
        (aero_conc_type==0?number_conc:1.0) <
        MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) continue;

    // Slow down rates as water approaches the minimum value
    realtype water_adj = state[AERO_WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    realtype water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    water_scaling *= water_scaling;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (UGM3_TO_PPM_ *
        EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

    // Slow down condensation rate as gas-phase concentrations become small
    realtype gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    realtype cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    cond_scaling *= cond_scaling;

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[GAS_SPEC_] * cond_scaling * water_scaling;

    // Slow down evaporation as aerosol-phase concentrations become small
    realtype aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    evap_scaling *= evap_scaling;

    // Calculate aerosol-phase evaporation rate (ug/m^3/s)
    evap_rate *= state[AERO_SPEC_(i_phase)] * evap_scaling * water_scaling;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (DERIV_ID_(0)>=0) {
      if (aero_conc_type==0) {
        // Scale the changes to the gas-phase by the number of particles for
        // per-particle aerosol concentrations
        //deriv[DERIV_ID_(0)] += number_conc *
           //     (evap_rate * UGM3_TO_PPM_ - cond_rate);
           atomicAdd((double*)&(deriv[DERIV_ID_(0)]),number_conc *
                (evap_rate * UGM3_TO_PPM_ - cond_rate));
      } else {
        // No scaling for aerosol concentrations with total mass per aerosol
        // phase
        //deriv[DERIV_ID_(0)] += evap_rate * UGM3_TO_PPM_ - cond_rate;
        atomicAdd((double*)&(deriv[DERIV_ID_(0)]),evap_rate * UGM3_TO_PPM_ - cond_rate);
      }
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    //if (DERIV_ID_(1+i_phase)>=0) deriv[DERIV_ID_(1+i_phase)] +=
      //      cond_rate / UGM3_TO_PPM_ - evap_rate;
      if (DERIV_ID_(1+i_phase)>=0)
            atomicAdd((double*)&(deriv[DERIV_ID_(1+i_phase)]),cond_rate / UGM3_TO_PPM_ - evap_rate);

  }

  //return (void*) &(float_data[FLOAT_DATA_SIZE_]);

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
__device__ void rxn_gpu_HL_phase_transfer_calc_jac_contrib(ModelDatagpu *model_data,
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
    aero_rep_gpu_get_effective_radius(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase),	// aerosol phase index
          &radius);			// particle effective radius (m)

    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_gpu_get_number_conc(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase),	// aerosol phase index
          &number_conc);		// particle number concentration
                                                // (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_gpu_get_aero_conc_type(
          model_data,			// model data
          AERO_REP_ID_(i_phase),	// aerosol representation index
          AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] *
        (aero_conc_type==0?number_conc:1.0) <
        MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) continue;

    // Slow down rates as water approaches the minimum value
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

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (UGM3_TO_PPM_ *
        EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

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

    // Slow down evaporation as aerosol-phase concentrations become small
    realtype aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    realtype evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    realtype evap_scaling_deriv =
      2.0 / ( SMALL_NUMBER_ * ( exp(  aero_adj / SMALL_NUMBER_ ) + 2.0 +
                                exp( -aero_adj / SMALL_NUMBER_ ) ) );
    evap_scaling_deriv *= 2.0 * evap_scaling;
    evap_scaling *= evap_scaling;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (aero_conc_type==0) {
      // Scale the changes to the gas-phase by the number of particles for
      // per-particle aerosol concentrations
      if (JAC_ID_(1+i_phase*5+1)>=0)
          J[JAC_ID_(1+i_phase*5+1)] += number_conc * evap_rate *
                      water_scaling * UGM3_TO_PPM_ *
                      ( evap_scaling + state[AERO_SPEC_(i_phase)] *
                                       evap_scaling_deriv );
      if (JAC_ID_(1+i_phase*5+3)>=0)
          J[JAC_ID_(1+i_phase*5+3)] -= number_conc * evap_rate *
                      UGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)] *
                      ( water_scaling / state[AERO_WATER_(i_phase)] +
                        water_scaling_deriv );
      if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= number_conc * cond_rate *
                      water_scaling *
                      ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
    } else {
      // No scaling for aerosol concentrations with total mass per aerosol phase
      if (JAC_ID_(1+i_phase*5+1)>=0)
              J[JAC_ID_(1+i_phase*5+1)] += evap_rate * water_scaling *
                      UGM3_TO_PPM_ *
                      ( evap_scaling + state[AERO_SPEC_(i_phase)] *
                                       evap_scaling_deriv );
      if (JAC_ID_(1+i_phase*5+3)>=0)
              J[JAC_ID_(1+i_phase*5+3)] -= evap_rate *
                      UGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)] *
                      ( water_scaling / state[AERO_WATER_(i_phase)] +
                        water_scaling_deriv );
      if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= cond_rate * water_scaling *
                      ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (JAC_ID_(1+i_phase*5)>=0) J[JAC_ID_(1+i_phase*5)] +=
            cond_rate * water_scaling / UGM3_TO_PPM_ *
            ( cond_scaling + state[GAS_SPEC_] * cond_scaling_deriv );
    if (JAC_ID_(1+i_phase*5+2)>=0) J[JAC_ID_(1+i_phase*5+2)] -= evap_rate *
            water_scaling *
            ( evap_scaling + state[AERO_SPEC_(i_phase)] * evap_scaling_deriv );
    if (JAC_ID_(1+i_phase*5+4)>=0) J[JAC_ID_(1+i_phase*5+4)] += evap_rate *
        state[AERO_SPEC_(i_phase)] *
            ( water_scaling / state[AERO_WATER_(i_phase)] +
              water_scaling_deriv );
  }

  //return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void *rxn_gpu_HL_phase_transfer_skip(void *rxn_data) {
  int *int_data = (int *) rxn_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void *rxn_gpu_HL_phase_transfer_print(void *rxn_data) {
  int *int_data = (int *) rxn_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nPhase Transfer reaction\n");
  for (int i = 0; i < INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i = 0; i < FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);

}

}