/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aqueous Equilibrium reaction solver functions
 *
*/
/** \file
 * \brief Aqueous Equilibrium reaction solver functions
*/

extern "C"{
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns_gpu.h"

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Small number
#define SMALL_NUMBER_ 1.0e-30//1.0e-30

// Factor used to calculate minimum water concentration for aqueous
// phase equilibrium reactions
#define MIN_WATER_ 1.0e-4

#define NUM_REACT_ (int_data[0*n_rxn])
#define NUM_PROD_ (int_data[1*n_rxn])
#define NUM_AERO_PHASE_ (int_data[2*n_rxn])
#define A_ (float_data[0*n_rxn])
#define C_ (float_data[1*n_rxn])
#define RATE_CONST_REVERSE_ (float_data[2*n_rxn])
#define RATE_CONST_FORWARD_ (rxn_env_data[0*n_rxn])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 3
#define REACT_(x) (int_data[(NUM_INT_PROP_+x)*n_rxn]-1)
#define PROD_(x) (int_data[(NUM_INT_PROP_+NUM_REACT_*NUM_AERO_PHASE_+x)*n_rxn]-1)
#define WATER_(x) (int_data[(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_)*NUM_AERO_PHASE_+x)*n_rxn]-1)
#define ACTIVITY_COEFF_(x) (int_data[(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+1)*NUM_AERO_PHASE_+x)*n_rxn]-1)
#define DERIV_ID_(x) (int_data[(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+2)*NUM_AERO_PHASE_+x)*n_rxn])
#define JAC_ID_(x) (int_data[(NUM_INT_PROP_+(2*(NUM_REACT_+NUM_PROD_)+2)*NUM_AERO_PHASE_+x)*n_rxn])
#define MASS_FRAC_TO_M_(x) (float_data[(NUM_FLOAT_PROP_+x)*n_rxn])
#define SMALL_WATER_CONC_(x) (float_data[(NUM_FLOAT_PROP_+NUM_REACT_+NUM_PROD_+x)*n_rxn])
#define SMALL_CONC_(x) (float_data[(NUM_FLOAT_PROP_+NUM_REACT_+NUM_PROD_+NUM_AERO_PHASE_+x)*n_rxn])

#ifdef __CUDA_ARCH__
__host__ __device__
#endif
realtype rxn_gpu_aqueous_equilibrium_calc_overall_rate(int *rxn_int_data,
     double *rxn_float_data, double *rxn_env_data, realtype *state,
     realtype react_fact, realtype prod_fact, realtype water, int i_phase, int n_rxn2)
{
#ifdef __CUDA_ARCH__
  int n_rxn=n_rxn2;
#else
  int n_rxn=1;;
#endif
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  realtype rate = ONE;
  realtype rc_forward = RATE_CONST_FORWARD_;
  realtype rc_reverse = RATE_CONST_REVERSE_;
  realtype react_fact_l = react_fact;
  realtype prod_fact_l = prod_fact;
  realtype water_l = water;

  /// \todo explore higher precision variables to reduce Jac errors

  // Calculate the overall rate
  // These equations are set up to try to avoid loss of accuracy from
  // subtracting two almost-equal numbers when rate_forward ~ rate_backward.
  // When modifying these calculations, be sure to use the Jacobian checker
  // during unit testing.
  if (react_fact_l == ZERO || prod_fact_l == ZERO) {
      rate = rc_forward * react_fact_l - rc_reverse * prod_fact_l;
  } else if (rc_forward * react_fact_l > rc_reverse * prod_fact_l) {
      realtype mod_react = ONE;
      for (int i_react = 1; i_react < NUM_REACT_; ++i_react) {
          mod_react *= (realtype)state[REACT_(i_phase * NUM_REACT_ + i_react)] *
                       (realtype)MASS_FRAC_TO_M_(i_react) / water_l;
      }
      realtype r1_eq = (prod_fact_l / mod_react) * (rc_reverse / rc_forward);
      rate = ((realtype)state[REACT_(i_phase * NUM_REACT_)] *
              (realtype)MASS_FRAC_TO_M_(0) / water_l -
              r1_eq) *
             rc_forward * mod_react;
  } else {
      realtype mod_prod = ONE;
      for (int i_prod = 1; i_prod < NUM_PROD_; ++i_prod) {
          mod_prod *= (realtype)state[PROD_(i_phase * NUM_PROD_ + i_prod)] *
                      (realtype)MASS_FRAC_TO_M_(NUM_REACT_ + i_prod) / water_l;
      }
      if (ACTIVITY_COEFF_(i_phase) >= 0)
          mod_prod *= (realtype)state[ACTIVITY_COEFF_(i_phase)];
      realtype p1_eq = (react_fact_l / mod_prod) * (rc_forward / rc_reverse);
      rate = (p1_eq - (realtype)state[PROD_(i_phase * NUM_PROD_)] *
                      (realtype)MASS_FRAC_TO_M_(NUM_REACT_) / water_l) *
             rc_reverse * mod_prod;
  }

return (realtype)rate;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
#ifdef __CUDA_ARCH__
__host__ __device__
#endif
void rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(ModelData *model_data, realtype *deriv, int *rxn_int_data,
          double *rxn_float_data, double *rxn_env_data, double time_step)
{
#ifdef __CUDA_ARCH__
  int n_rxn=model_data->n_rxn;
#else
  int n_rxn=1;
#endif
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If no aerosol water is present, no reaction occurs
    double water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Get the product of all reactants
    double react_fact = ONE;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      react_fact *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / water;
    }

    // Get the product of all products
    double prod_fact = ONE;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      prod_fact *= state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / water;
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) prod_fact *=
            state[ACTIVITY_COEFF_(i_phase)];

    double rate = rxn_gpu_aqueous_equilibrium_calc_overall_rate(
                        rxn_int_data, rxn_float_data, rxn_env_data,
                        state, react_fact, prod_fact, water, i_phase, n_rxn);
    if (rate == ZERO) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(deriv[DERIV_ID_(i_deriv++)]), -(rate /
        MASS_FRAC_TO_M_(i_react)) * water);
#else
      deriv[DERIV_ID_(i_deriv++)] -= (rate / MASS_FRAC_TO_M_(i_react)) * water;
#endif
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(deriv[DERIV_ID_(i_deriv++)]), (rate /
        MASS_FRAC_TO_M_(NUM_REACT_+i_prod)) * water);
#else
      deriv[DERIV_ID_(i_deriv++)] +=
          (rate / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod)) * water;;
#endif

    }

  }

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */

#ifdef PMC_USE_SUNDIALS
#ifdef __CUDA_ARCH__
__host__ __device__
#endif
void rxn_gpu_aqueous_equilibrium_calc_jac_contrib(ModelData *model_data, realtype *J, int *rxn_int_data,
          double *rxn_float_data, double *rxn_env_data, double time_step)
{
#ifdef __CUDA_ARCH__
  int n_rxn=model_data->n_rxn;
#else
  int n_rxn=1;;
#endif
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  double water;
  double forward_rate;
  double reverse_rate;

  //todo: fix gpu by following the methodology of calc_deriv

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase=0, i_jac = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If not aerosol water is present, no reaction occurs
    water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + NUM_PROD_ + 2);
      continue;
    }

    // Calculate the forward rate (M/s)
    forward_rate = RATE_CONST_FORWARD_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      forward_rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / water;
    }

    // Calculate the reverse rate (M/s)
    reverse_rate = RATE_CONST_REVERSE_;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      reverse_rate *= state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / water;
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) reverse_rate *=
            state[ACTIVITY_COEFF_(i_phase)];

    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = 0; i_react_ind < NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (-forward_rate) /
        state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water);
#else
      J[JAC_ID_(i_jac++)] +=
        (-forward_rate) /
        state[REACT_(i_phase * NUM_REACT_ + i_react_ind)] /
        MASS_FRAC_TO_M_(i_react_dep) * water;
#endif
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (forward_rate) /
        state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
#else
      J[JAC_ID_(i_jac++)] +=
          (forward_rate) / state[REACT_(i_phase * NUM_REACT_ + i_react_ind)] /
          MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
#endif
      }
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = 0; i_prod_ind < NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (reverse_rate) /
        state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water);
#else
      J[JAC_ID_(i_jac++)] += (reverse_rate) /
        state[PROD_(i_phase * NUM_PROD_ + i_prod_ind)] /
        MASS_FRAC_TO_M_(i_react_dep) * water;
#endif
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
        (-reverse_rate) /
          state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
#else
      J[JAC_ID_(i_jac++)] += (-reverse_rate) /
       state[PROD_(i_phase * NUM_PROD_ + i_prod_ind)] /
       MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
#endif
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (reverse_rate) /
         ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(i_react_dep));
#else
      J[JAC_ID_(i_jac++)] +=
          (forward_rate * (NUM_REACT_ - 1) - reverse_rate * (NUM_PROD_ - 1)) /
          MASS_FRAC_TO_M_(i_react_dep);
#endif
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
         -( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep));
#else
      J[JAC_ID_(i_jac++)] -=
          (forward_rate * (NUM_REACT_ - 1) - reverse_rate * (NUM_PROD_ - 1)) /
          MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep);
#endif
    }

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase)<0) {
      i_jac += NUM_REACT_ + NUM_PROD_;
      continue;
    }
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
          reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(i_react_dep) * water);
#else
      J[JAC_ID_(i_jac++)] += reverse_rate / state[ACTIVITY_COEFF_(i_phase)] /
                             MASS_FRAC_TO_M_(i_react_dep) * water;
#endif
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
#ifdef __CUDA_ARCH__
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
          -reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
#else
      J[JAC_ID_(i_jac++)] -= reverse_rate / state[ACTIVITY_COEFF_(i_phase)] /
                       MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
#endif
    }

  }

}
#endif


}