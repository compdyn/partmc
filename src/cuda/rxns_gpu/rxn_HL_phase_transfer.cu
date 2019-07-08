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
#define TEMPERATURE_K_ env_data[0*n_rxn]
#define PRESSURE_PA_ env_data[1*n_rxn]

// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472
// Small number for ignoring low concentrations
#define VERY_SMALL_NUMBER_ 1.0e-30
// Factor used to calculate minimum aerosol water concentrations for
// HL phase transfer
#define MIN_WATER_ 1.0e-4

#define DELTA_H_ float_data[0*n_rxn]
#define DELTA_S_ float_data[1*n_rxn]
#define DIFF_COEFF_ float_data[2*n_rxn]
#define PRE_C_AVG_ float_data[3*n_rxn]
#define A_ float_data[4*n_rxn]
#define C_ float_data[5*n_rxn]
#define C_AVG_ALPHA_ float_data[6*n_rxn]
#define EQUIL_CONST_ float_data[7*n_rxn]
#define CONV_ float_data[8*n_rxn]
#define MW_ float_data[9*n_rxn]
#define UGM3_TO_PPM_ float_data[10*n_rxn]
#define SMALL_NUMBER_ float_data[11*n_rxn]
#define NUM_AERO_PHASE_ int_data[0*n_rxn]
#define GAS_SPEC_ (int_data[1*n_rxn]-1)
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 12
#define DERIV_ID_(x) int_data[(NUM_INT_PROP_ + x)*n_rxn]
#define JAC_ID_(x) int_data[(NUM_INT_PROP_ + 1 + NUM_AERO_PHASE_ + x)*n_rxn]
#define PHASE_INT_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 6*NUM_AERO_PHASE_ + x)*n_rxn]-1)
#define PHASE_REAL_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 7*NUM_AERO_PHASE_ + x)*n_rxn]-1)
#define AERO_SPEC_(x) (int_data[(PHASE_INT_LOC_(x))*n_rxn]-1)
#define AERO_WATER_(x) (int_data[(PHASE_INT_LOC_(x) + 1)*n_rxn]-1)
#define AERO_PHASE_ID_(x) (int_data[(PHASE_INT_LOC_(x) + 2)*n_rxn]-1)
#define AERO_REP_ID_(x) (int_data[(PHASE_INT_LOC_(x) + 3)*n_rxn]-1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[(PHASE_INT_LOC_(x) + 4)*n_rxn])
#define PHASE_JAC_ID_(x, s, e) int_data[(PHASE_INT_LOC_(x) + 5 + s*NUM_AERO_PHASE_JAC_ELEM_(x) + e)*n_rxn]
#define SMALL_WATER_CONC_(x) (float_data[(PHASE_REAL_LOC_(x))*n_rxn])
#define FAST_FLUX_(x) (float_data[(PHASE_REAL_LOC_(x) + 1)*n_rxn])
#define AERO_ADJ_(x) (float_data[(PHASE_REAL_LOC_(x) + 2)*n_rxn])
#define EFF_RAD_JAC_ELEM_(x, e) float_data[(PHASE_REAL_LOC_(x) + 3 + e)*n_rxn]
#define NUM_CONC_JAC_ELEM_(x, e) float_data[(PHASE_REAL_LOC_(x) + 3 + NUM_AERO_PHASE_JAC_ELEM_(x) + e)*n_rxn]
#define INT_DATA_SIZE_ (PHASE_INT_LOC_(NUM_AERO_PHASE_-1)+5+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))
#define FLOAT_DATA_SIZE_ (PHASE_REAL_LOC_(NUM_AERO_PHASE_-1)+3+2*NUM_AERO_PHASE_JAC_ELEM_(NUM_AERO_PHASE_-1))

/** \brief Update reaction data for new environmental conditions
 *
 * For Phase Transfer reaction this only involves recalculating the rate
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
__device__ void rxn_gpu_HL_phase_transfer_update_env_state(int n_rxn2, double *double_pointer_gpu, double *env_data,
                                              void *rxn_data)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

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
__device__ void rxn_gpu_HL_phase_transfer_calc_deriv_contrib(ModelDatagpu *model_data, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int deriv_length, int n_rxn2)
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
    double water_adj = state[AERO_WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    double water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    water_scaling *= water_scaling;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    double cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    double evap_rate = cond_rate / (UGM3_TO_PPM_ *
        EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

    // Slow down condensation rate as gas-phase concentrations become small
    double gas_adj = state[GAS_SPEC_] - VERY_SMALL_NUMBER_;
    gas_adj = ( gas_adj > ZERO ) ? gas_adj : ZERO;
    double cond_scaling =
      2.0 / ( 1.0 + exp( -gas_adj / SMALL_NUMBER_ ) ) - 1.0;
    cond_scaling *= cond_scaling;

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[GAS_SPEC_] * cond_scaling * water_scaling;

    // Slow down evaporation as aerosol-phase concentrations become small
    double aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    double evap_scaling =
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
*/
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
__device__ void rxn_gpu_HL_phase_transfer_calc_jac_contrib(ModelDatagpu *model_data, double *state,
          double *J, void *rxn_data, double * double_pointer_gpu, double time_step, int deriv_length, int n_rxn2)
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
    double water_adj = state[AERO_WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    double water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    double water_scaling_deriv =
      2.0 / ( SMALL_WATER_CONC_(i_phase) *
              ( exp(  water_adj / SMALL_WATER_CONC_(i_phase) ) + 2.0 +
                exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) );
    water_scaling_deriv *= 2.0 * water_scaling;
    water_scaling *= water_scaling;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    double cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) +
              4.0*radius/(3.0*C_AVG_ALPHA_));

    // Calculate the evaporation rate constant (1/s)
    double evap_rate = cond_rate / (UGM3_TO_PPM_ *
        EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

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

    // Slow down evaporation as aerosol-phase concentrations become small
    double aero_adj = state[AERO_SPEC_(i_phase)] - VERY_SMALL_NUMBER_;
    aero_adj = ( aero_adj > ZERO ) ? aero_adj : ZERO;
    double evap_scaling =
      2.0 / ( 1.0 + exp( -aero_adj / SMALL_NUMBER_ ) ) - 1.0;
    double evap_scaling_deriv =
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

*/

}
#endif


/** \brief Retrieve Int data size
 *
 * \param rxn_data Pointer to the reaction data
 * \return The data size of int array
 */
void * rxn_gpu_HL_phase_transfer_int_size(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) float_data;
}

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void *rxn_gpu_HL_phase_transfer_skip(void *rxn_data) {

  int n_rxn=1;
  int *int_data = (int *) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void *rxn_gpu_HL_phase_transfer_print(void *rxn_data) {
  int n_rxn=1;
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