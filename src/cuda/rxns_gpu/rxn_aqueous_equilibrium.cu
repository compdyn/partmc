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
#define RATE_CONST_FORWARD_ (rate_constants[0*n_rxn])
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

/** \brief Update reaction data for new environmental conditions
 *
 * For Aqueous Equilibrium reaction this only involves recalculating the
 * forward rate constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
__device__ void rxn_gpu_aqueous_equilibrium_update_env_state(double *rate_constants,
   int n_rxn2,double *double_pointer_gpu, double *env_data,
          void *rxn_data)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the equilibrium constant
  // (assumes reactant and product concentrations in M)
  double equil_const;
  if (C_==0.0) {
    equil_const = A_;
  } else {
    equil_const = A_ * exp(C_ * (1.0/TEMPERATURE_K_ - 1.0/298.0));
  }

  // Set the forward rate constant
  RATE_CONST_FORWARD_ = equil_const * RATE_CONST_REVERSE_;

  rate_constants[0] = RATE_CONST_FORWARD_;
}


__device__ double rxn_gpu_aqueous_equilibrium_calc_overall_rate(int *rxn_data,
     double *double_pointer_gpu, double *rate_constants, double *state,
     double react_fact, double prod_fact, double water, int i_phase, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  double rate         = ONE;
  double rc_forward   = RATE_CONST_FORWARD_;
  double rc_reverse   = RATE_CONST_REVERSE_;
  double react_fact_l = react_fact;
  double prod_fact_l  = prod_fact;
  double water_l      = water;

  /// \todo explore higher precision variables to reduce Jac errors

  // Calculate the overall rate
  // These equations are set up to try to avoid loss of accuracy from
  // subtracting two almost-equal numbers when rate_forward ~ rate_backward.
  // When modifying these calculations, be sure to use the Jacobian checker
  // during unit testing.
  if (react_fact_l == ZERO || prod_fact_l == ZERO) {
    rate = rc_forward * react_fact_l - rc_reverse * prod_fact_l;
  } else if (rc_forward * react_fact_l >
             rc_reverse * prod_fact_l) {
    double mod_react = ONE;
    for (int i_react = 1; i_react < NUM_REACT_; ++i_react) {
      mod_react *= (double) state[REACT_(i_phase*NUM_REACT_+i_react)] *
                   (double) MASS_FRAC_TO_M_(i_react) / water_l;
    }
    double r1_eq = (prod_fact_l / mod_react) *
                     (rc_reverse / rc_forward);
    rate = ((double) state[REACT_(i_phase*NUM_REACT_)] *
            (double) MASS_FRAC_TO_M_(0) / water_l - r1_eq) *
           rc_forward * mod_react;
  } else {
    double mod_prod = ONE;
    for (int i_prod = 1; i_prod < NUM_PROD_; ++i_prod) {
      mod_prod *= (double) state[PROD_(i_phase*NUM_PROD_+i_prod)] *
                  (double) MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / water_l;
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) mod_prod *=
                                             (double) state[ACTIVITY_COEFF_(i_phase)];
    double p1_eq = (react_fact_l / mod_prod) *
                     (rc_forward / rc_reverse);
    rate = (p1_eq - (double) state[PROD_(i_phase*NUM_PROD_)] *
                    (double) MASS_FRAC_TO_M_(NUM_REACT_) / water_l) *
           rc_reverse * mod_prod;
  }

  return (double) rate;
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

#ifdef PMC_USE_GPU
__device__
#endif
void rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

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
                        (int*) rxn_data, double_pointer_gpu, rate_constants,
                        state, react_fact, prod_fact, water, i_phase, n_rxn2);
    if (rate == ZERO) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      atomicAdd(&(deriv[DERIV_ID_(i_deriv++)]), -(rate /
        MASS_FRAC_TO_M_(i_react)) * water);
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      atomicAdd(&(deriv[DERIV_ID_(i_deriv++)]), (rate /
        MASS_FRAC_TO_M_(NUM_REACT_+i_prod)) * water);

    }

  }

}
#endif


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
void rxn_cpu_aqueous_equilibrium_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If no aerosol water is present, no reaction occurs
    if (state[WATER_(i_phase)] < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Slow down rates as water approaches the minimum value
    double water_adj = state[WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    double water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;

    // Get the lowest concentration to use in slowing rates
    double min_react_conc = 1.0e100;
    double min_prod_conc  = 1.0e100;

    // Calculate the forward rate (M/s)
    double forward_rate = RATE_CONST_FORWARD_ * water_scaling;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      forward_rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / state[WATER_(i_phase)];
      if (min_react_conc > state[REACT_(i_phase*NUM_REACT_+i_react)])
        min_react_conc = state[REACT_(i_phase*NUM_REACT_+i_react)];
    }

    // Calculate the reverse rate (M/s)
    double reverse_rate = RATE_CONST_REVERSE_ * water_scaling;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      reverse_rate *= state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / state[WATER_(i_phase)];

      if (min_prod_conc > state[PROD_(i_phase*NUM_PROD_+i_prod)])
        min_prod_conc = state[PROD_(i_phase*NUM_PROD_+i_prod)];
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) reverse_rate *=
            state[ACTIVITY_COEFF_(i_phase)];

    // Slow rates as concentrations become low

    double min_conc =
      (forward_rate > reverse_rate) ? min_react_conc : min_prod_conc;
    min_conc -= SMALL_NUMBER_;

    if (min_conc <= ZERO) continue;

    //continue is causing fail
    double spec_scaling =
      2.0 / ( 1.0 + exp( -min_conc / SMALL_CONC_(i_phase) ) ) - 1.0;
    forward_rate *= spec_scaling;
    reverse_rate *= spec_scaling;

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] += (reverse_rate - forward_rate) /
	      MASS_FRAC_TO_M_(i_react) * state[WATER_(i_phase)];
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] += (forward_rate - reverse_rate) /
	      MASS_FRAC_TO_M_(NUM_REACT_+i_prod) * state[WATER_(i_phase)];

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
__device__ void rxn_gpu_aqueous_equilibrium_calc_jac_contrib(double *rate_constants, double *state,
          double *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  double water;
  double forward_rate;
  double reverse_rate;

  //todo: check if works

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
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (-forward_rate) /
        state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water);
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (forward_rate) /
        state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
      }
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = 0; i_prod_ind < NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (reverse_rate) /
        state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water);
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
        (-reverse_rate) /
          state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
      (reverse_rate) /
         ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(i_react_dep));
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
         -( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep));
    }

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase)<0) {
      i_jac += NUM_REACT_ + NUM_PROD_;
      continue;
    }
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
          reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(i_react_dep) * water);
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd(&(J[JAC_ID_(i_jac++)]),
          -reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
    }

  }


  /*
  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase=0, i_jac = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If not aerosol water is present, no reaction occurs
    if (state[WATER_(i_phase)] < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + NUM_PROD_ + 1);
      continue;
    }

    // Slow down rates as water approaches the minimum value
    double water_adj = state[WATER_(i_phase)] -
                         MIN_WATER_ * SMALL_WATER_CONC_(i_phase);
    water_adj = ( water_adj > ZERO ) ? water_adj : ZERO;
    double water_scaling =
      2.0 / ( 1.0 + exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) - 1.0;
    double water_scaling_deriv =
      2.0 / ( SMALL_WATER_CONC_(i_phase) *
              ( exp(  water_adj / SMALL_WATER_CONC_(i_phase) ) + 2.0 +
                exp( -water_adj / SMALL_WATER_CONC_(i_phase) ) ) );

    // Get the lowest concentration to use in slowing rates
    double min_react_conc = 1.0e100;
    double min_prod_conc  = 1.0e100;
    int low_react_id = 0;
    int low_prod_id  = 0;

    // Calculate the forward rate (M/s)
    double forward_rate = RATE_CONST_FORWARD_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      forward_rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / state[WATER_(i_phase)];
      if (min_react_conc > state[REACT_(i_phase*NUM_REACT_+i_react)]) {
        min_react_conc = state[REACT_(i_phase*NUM_REACT_+i_react)];
        low_react_id = i_react;
      }
    }

    // Calculate the reverse rate (M/s)
    double reverse_rate = RATE_CONST_REVERSE_;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      reverse_rate *= state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / state[WATER_(i_phase)];
      if (min_prod_conc > state[PROD_(i_phase*NUM_PROD_+i_prod)]) {
        min_prod_conc = state[PROD_(i_phase*NUM_PROD_+i_prod)];
        low_prod_id = NUM_REACT_ + i_prod;
      }
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) reverse_rate *=
            state[ACTIVITY_COEFF_(i_phase)];

    // Slow rates as concentrations become low
    double min_conc;
    int low_spec_id;
    if (forward_rate > reverse_rate) {
      min_conc = min_react_conc;
      low_spec_id = low_react_id;
    } else {
      min_conc = min_prod_conc;
      low_spec_id = low_prod_id;
    }
    min_conc -= SMALL_NUMBER_;
    if (min_conc <= ZERO) continue;
    double spec_scaling =
      2.0 / ( 1.0 + exp( -min_conc / SMALL_CONC_(i_phase) ) ) - 1.0;
    double spec_scaling_deriv =
      2.0 / ( SMALL_CONC_(i_phase) *
              ( exp(  min_conc / SMALL_CONC_(i_phase) ) + 2.0 +
                exp( -min_conc / SMALL_CONC_(i_phase) ) ) );

    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = 0; i_react_ind < NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
        if (low_spec_id == i_react_ind) {
          J[JAC_ID_(i_jac)] += (-forward_rate) * water_scaling *
                               spec_scaling_deriv;
        }
        J[JAC_ID_(i_jac++)] += (-forward_rate) * water_scaling * spec_scaling /
                state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * state[WATER_(i_phase)];
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
        if (low_spec_id == i_react_ind) {
          J[JAC_ID_(i_jac)] += (forward_rate) * water_scaling *
                               spec_scaling_deriv;
        }
        J[JAC_ID_(i_jac++)] += (forward_rate) * water_scaling * spec_scaling /
                state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) *
                state[WATER_(i_phase)];
      }
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = 0; i_prod_ind < NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        if (low_spec_id == NUM_REACT_ + i_prod_ind) {
          J[JAC_ID_(i_jac)] += (reverse_rate) * water_scaling *
                               spec_scaling_deriv;
        }
        J[JAC_ID_(i_jac++)] += (reverse_rate) * water_scaling * spec_scaling /
                state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * state[WATER_(i_phase)];
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        if (low_spec_id == NUM_REACT_ + i_prod_ind) {
          J[JAC_ID_(i_jac)] += (-reverse_rate) * water_scaling *
                               spec_scaling_deriv;
        }
        J[JAC_ID_(i_jac++)] += (-reverse_rate) * water_scaling * spec_scaling /
                state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) *
                state[WATER_(i_phase)];
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] +=
        ( ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
          * water_scaling * spec_scaling / state[WATER_(i_phase)] +
          ( forward_rate - reverse_rate ) * spec_scaling * water_scaling_deriv
        ) / MASS_FRAC_TO_M_(i_react_dep);
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] -=
        ( ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
          * water_scaling * spec_scaling / state[WATER_(i_phase)] +
          ( forward_rate - reverse_rate ) * spec_scaling * water_scaling_deriv
        ) / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep);
    }

  }

   */

}
#endif

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef SMALL_NUMBER_

#undef MIN_WATER_

#undef NUM_REACT_
#undef NUM_PROD_
#undef NUM_AERO_PHASE_
#undef A_
#undef C_
#undef RATE_CONST_REVERSE_
#undef RATE_CONST_FORWARD_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef PROD_
#undef WATER_
#undef ACTIVITY_COEFF_
#undef DERIV_ID_
#undef JAC_ID_
#undef MASS_FRAC_TO_M_
#undef SMALL_WATER_CONC_
#undef SMALL_CONC_
}