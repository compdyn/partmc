/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aqueous Equilibrium reaction solver functions
 *
*/
/** \file
 * \brief Aqueous Equilibrium reaction solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Small number
#define SMALL_NUMBER_ 1.0e-30

// Factor used to calculate minimum water concentration for aqueous
// phase equilibrium reactions
#define MIN_WATER_ 1.0e-4

#define NUM_REACT_ (int_data[0])
#define NUM_PROD_ (int_data[1])
#define NUM_AERO_PHASE_ (int_data[2])
#define A_ (float_data[0])
#define C_ (float_data[1])
#define RATE_CONST_REVERSE_ (float_data[2])
#define RATE_CONST_FORWARD_ (float_data[3])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 4
#define REACT_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PROD_(x) (int_data[NUM_INT_PROP_+NUM_REACT_*NUM_AERO_PHASE_+x]-1)
#define WATER_(x) (int_data[NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_)*NUM_AERO_PHASE_+x]-1)
#define ACTIVITY_COEFF_(x) (int_data[NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+1)*NUM_AERO_PHASE_+x]-1)
#define DERIV_ID_(x) (int_data[NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+2)*NUM_AERO_PHASE_+x])
#define JAC_ID_(x) (int_data[NUM_INT_PROP_+(2*(NUM_REACT_+NUM_PROD_)+2)*NUM_AERO_PHASE_+x])
#define MASS_FRAC_TO_M_(x) (float_data[NUM_FLOAT_PROP_+x])
#define SMALL_WATER_CONC_(x) (float_data[NUM_FLOAT_PROP_+NUM_REACT_+NUM_PROD_+x])
#define SMALL_CONC_(x) (float_data[NUM_FLOAT_PROP_+NUM_REACT_+NUM_PROD_+NUM_AERO_PHASE_+x])

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 */
void rxn_aqueous_equilibrium_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
          bool **jac_struct)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Loop over all the instances of the specified phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {

    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = i_phase*NUM_REACT_;
              i_react_ind < (i_phase+1)*NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][REACT_(i_react_ind)] = true;
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][REACT_(i_react_ind)] = true;
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = i_phase*NUM_PROD_;
              i_prod_ind < (i_phase+1)*NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][PROD_(i_prod_ind)] = true;
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][PROD_(i_prod_ind)] = true;
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
      jac_struct[REACT_(i_react_dep)][WATER_(i_phase)] = true;
    for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
      jac_struct[PROD_(i_prod_dep)][WATER_(i_phase)] = true;

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase)<0) continue;
    for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; ++i_react_dep)
      jac_struct[REACT_(i_react_dep)][ACTIVITY_COEFF_(i_phase)] = true;
    for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; ++i_prod_dep)
      jac_struct[PROD_(i_prod_dep)][ACTIVITY_COEFF_(i_phase)] = true;

  }

  return;
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_aqueous_equilibrium_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, int *rxn_int_data, double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  for (int i_phase = 0, i_deriv = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    for (int i_react = 0; i_react < NUM_REACT_; i_react++)
      DERIV_ID_(i_deriv++) = deriv_ids[REACT_(i_phase*NUM_REACT_+i_react)];
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++)
      DERIV_ID_(i_deriv++) = deriv_ids[PROD_(i_phase*NUM_PROD_+i_prod)];
  }

  // Update the Jacobian ids
  for (int i_phase = 0, i_jac = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {

    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = i_phase*NUM_REACT_;
              i_react_ind < (i_phase+1)*NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][REACT_(i_react_ind)];
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][REACT_(i_react_ind)];
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = i_phase*NUM_PROD_;
              i_prod_ind < (i_phase+1)*NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][PROD_(i_prod_ind)];
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][PROD_(i_prod_ind)];
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
      JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][WATER_(i_phase)];
    for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
      JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][WATER_(i_phase)];

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase)<0) {
      for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; ++i_react_dep)
        JAC_ID_(i_jac++) = -1;
      for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; ++i_prod_dep)
        JAC_ID_(i_jac++) = -1;
    } else {
      for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; ++i_react_dep)
        JAC_ID_(i_jac++) =
          jac_ids[REACT_(i_react_dep)][ACTIVITY_COEFF_(i_phase)];
      for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; ++i_prod_dep)
        JAC_ID_(i_jac++) =
          jac_ids[PROD_(i_prod_dep)][ACTIVITY_COEFF_(i_phase)];
    }

  }

  // Calculate a small concentration for aerosol-phase species based on the
  // integration tolerances to use during solving. TODO find a better place
  // to do this
  double *abs_tol = (double*) model_data->abs_tol;
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++ ) {
    SMALL_CONC_(i_phase) = 99999.0;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (SMALL_CONC_(i_phase) >
          abs_tol[REACT_(i_phase*NUM_REACT_+i_react)] / 100.0)
          SMALL_CONC_(i_phase) =
            abs_tol[REACT_(i_phase*NUM_REACT_+i_react)] / 100.0;
    }
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (SMALL_CONC_(i_phase) >
          abs_tol[PROD_(i_phase*NUM_PROD_+i_prod)] / 100.0)
          SMALL_CONC_(i_phase) =
            abs_tol[PROD_(i_phase*NUM_PROD_+i_prod)] / 100.0;
    }
  }

  // Calculate a small concentration for aerosol-phase water based on the
  // integration tolerances to use during solving. TODO find a better place
  // to do this
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    SMALL_WATER_CONC_(i_phase) = abs_tol[WATER_(i_phase)] / 10.0;
  }

  return;

}

/** \brief Update reaction data for new environmental conditions
 *
 * For Aqueous Equilibrium reaction this only involves recalculating the
 * forward rate constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_aqueous_equilibrium_update_env_state(double *rate_constants,
    ModelData *model_data, int *rxn_int_data, double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

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

  // Set forward rate constant (reverse rate constant is constant)
  rate_constants[0] = RATE_CONST_FORWARD_;

  return;
}

/** \brief Calculate the overall per-particle reaction rate [M/s]
 *
 * This function is called by the deriv and Jac functions to get the overall
 * reaction rate on a per-particle basis, trying to avoid floating-point
 * subtraction errors.
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param state State array
 * \param react_fact Product of all reactant concentrations [M^n_react)
 * \param prod_fact Product of all product concentrations [M^n_prod]
 * \param water Water concentration [ug/m3]
 * \param i_phase Index for the aerosol phase being calcualted
 */
#ifdef PMC_USE_SUNDIALS
realtype rxn_aqueous_equilibrium_calc_overall_rate(int *rxn_int_data, double *rxn_float_data,
          realtype *state, realtype react_fact, realtype prod_fact,
          realtype water, int i_phase)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  realtype rate         = ONE;
  realtype rc_forward   = RATE_CONST_FORWARD_;
  realtype rc_reverse   = RATE_CONST_REVERSE_;
  realtype react_fact_l = react_fact;
  realtype prod_fact_l  = prod_fact;
  realtype water_l      = water;

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
    realtype mod_react = ONE;
    for (int i_react = 1; i_react < NUM_REACT_; ++i_react) {
      mod_react *= (realtype) state[REACT_(i_phase*NUM_REACT_+i_react)] *
              (realtype) MASS_FRAC_TO_M_(i_react) / water_l;
    }
    realtype r1_eq = (prod_fact_l / mod_react) *
                        (rc_reverse / rc_forward);
    rate = ((realtype) state[REACT_(i_phase*NUM_REACT_)] *
            (realtype) MASS_FRAC_TO_M_(0) / water_l - r1_eq) *
           rc_forward * mod_react;
  } else {
    realtype mod_prod = ONE;
    for (int i_prod = 1; i_prod < NUM_PROD_; ++i_prod) {
      mod_prod *= (realtype) state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              (realtype) MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / water_l;
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) mod_prod *=
            (realtype) state[ACTIVITY_COEFF_(i_phase)];
    realtype p1_eq = (react_fact_l / mod_prod) *
                        (rc_forward / rc_reverse);
    rate = (p1_eq - (realtype) state[PROD_(i_phase*NUM_PROD_)] *
            (realtype) MASS_FRAC_TO_M_(NUM_REACT_) / water_l) *
           rc_reverse * mod_prod;
  }

  return (realtype) rate;
}
#endif

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param time_step Current time step of the itegrator (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_aqueous_equilibrium_calc_deriv_contrib(double *rate_constants,
          ModelData *model_data,
          realtype *deriv, int *rxn_int_data, double *rxn_float_data, double time_step)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state    = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;
  int cell_id      = model_data->grid_cell_id;

  // Update the forward rate constant for this grid cell
  RATE_CONST_FORWARD_ = rate_constants[0];

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If no aerosol water is present, no reaction occurs
    realtype water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Get the product of all reactants
    realtype react_fact = ONE;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      react_fact *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / water;
    }

    // Get the product of all products
    realtype prod_fact = ONE;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      prod_fact *= state[PROD_(i_phase*NUM_PROD_+i_prod)] *
              MASS_FRAC_TO_M_(NUM_REACT_+i_prod) / water;
    }
    if (ACTIVITY_COEFF_(i_phase)>=0) prod_fact *=
            state[ACTIVITY_COEFF_(i_phase)];

    realtype rate = rxn_aqueous_equilibrium_calc_overall_rate(
                        rxn_int_data, rxn_float_data,
                        state, react_fact, prod_fact, water, i_phase);
    if (rate == ZERO) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] -= (rate / MASS_FRAC_TO_M_(i_react)) *
        water;
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] += (rate / MASS_FRAC_TO_M_(NUM_REACT_+i_prod)) *
        water;
    }

  }

  return;

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param time_step Current time step of the itegrator (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_aqueous_equilibrium_calc_jac_contrib(double *rate_constants,
          ModelData *model_data,
          realtype *J, int *rxn_int_data, double *rxn_float_data, double time_step)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state    = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;
  int cell_id      = model_data->grid_cell_id;

  // Update the forward rate constant for this grid cell
  RATE_CONST_FORWARD_ = rate_constants[0];

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase=0, i_jac = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If not aerosol water is present, no reaction occurs
    realtype water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + NUM_PROD_ + 2);
      continue;
    }

    // Calculate the forward rate (M/s)
    realtype forward_rate = RATE_CONST_FORWARD_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      forward_rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              MASS_FRAC_TO_M_(i_react) / water;
    }

    // Calculate the reverse rate (M/s)
    realtype reverse_rate = RATE_CONST_REVERSE_;
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
        J[JAC_ID_(i_jac++)] += (-forward_rate) /
                state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water;
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
        J[JAC_ID_(i_jac++)] += (forward_rate) /
                state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
      }
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = 0; i_prod_ind < NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        J[JAC_ID_(i_jac++)] += (reverse_rate) /
                state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(i_react_dep) * water;
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        J[JAC_ID_(i_jac++)] += (-reverse_rate) /
                state[PROD_(i_phase*NUM_PROD_+i_prod_ind)] /
		MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] +=
         ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(i_react_dep);
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] -=
         ( forward_rate * (NUM_REACT_-1) - reverse_rate * (NUM_PROD_-1) )
         / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep);
    }

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase)<0) {
      i_jac += NUM_REACT_ + NUM_PROD_;
      continue;
    }
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] += reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(i_react_dep) * water;
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[JAC_ID_(i_jac++)] -= reverse_rate
            / state[ACTIVITY_COEFF_(i_phase)]
            / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water;
    }

  }

  return;

}
#endif

/** \brief Print the Aqueous Equilibrium reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_aqueous_equilibrium_print(int *rxn_int_data, double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nAqueous Equilibrium reaction\n");

  return;
}
