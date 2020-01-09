/* Copyright (C) 2015-2019 Matthew Dawson
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
#define WATER_CONC_ (float_data[3])
#define ACTIVITY_COEFF_VALUE_ (float_data[4])
#define RATE_CONST_FORWARD_ (rxn_env_data[0])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 5
#define NUM_ENV_PARAM_ 1
#define REACT_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define PROD_(x) \
  (int_data[NUM_INT_PROP_ + NUM_REACT_ * NUM_AERO_PHASE_ + x] - 1)
#define WATER_(x) \
  (int_data[NUM_INT_PROP_ + (NUM_REACT_ + NUM_PROD_) * NUM_AERO_PHASE_ + x] - 1)
#define ACTIVITY_COEFF_(x)                                                   \
  (int_data[NUM_INT_PROP_ + (NUM_REACT_ + NUM_PROD_ + 1) * NUM_AERO_PHASE_ + \
            x] -                                                             \
   1)
#define DERIV_ID_(x) \
  (int_data[NUM_INT_PROP_ + (NUM_REACT_ + NUM_PROD_ + 2) * NUM_AERO_PHASE_ + x])
#define JAC_ID_(x)          \
  (int_data[NUM_INT_PROP_ + \
            (2 * (NUM_REACT_ + NUM_PROD_) + 2) * NUM_AERO_PHASE_ + x])
#define MASS_FRAC_TO_M_(x) (float_data[NUM_FLOAT_PROP_ + x])
#define REACT_CONC_(x) \
  (float_data[NUM_FLOAT_PROP_ + NUM_REACT_ + NUM_PROD_ + x])
#define PROD_CONC_(x) \
  (float_data[NUM_FLOAT_PROP_ + 2 * NUM_REACT_ + NUM_PROD_ + x])
#define SMALL_WATER_CONC_(x) \
  (float_data[NUM_FLOAT_PROP_ + 2 * NUM_REACT_ + 2 * NUM_PROD_ + x])
#define SMALL_CONC_(x)                                           \
  (float_data[NUM_FLOAT_PROP_ + 2 * NUM_REACT_ + 2 * NUM_PROD_ + \
              NUM_AERO_PHASE_ + x])

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 */
void rxn_aqueous_equilibrium_get_used_jac_elem(int *rxn_int_data,
                                               double *rxn_float_data,
                                               bool **jac_struct) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Loop over all the instances of the specified phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = i_phase * NUM_REACT_;
         i_react_ind < (i_phase + 1) * NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][REACT_(i_react_ind)] = true;
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][REACT_(i_react_ind)] = true;
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = i_phase * NUM_PROD_;
         i_prod_ind < (i_phase + 1) * NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][PROD_(i_prod_ind)] = true;
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][PROD_(i_prod_ind)] = true;
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase * NUM_REACT_;
         i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
      jac_struct[REACT_(i_react_dep)][WATER_(i_phase)] = true;
    for (int i_prod_dep = i_phase * NUM_PROD_;
         i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
      jac_struct[PROD_(i_prod_dep)][WATER_(i_phase)] = true;

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase) < 0) continue;
    for (int i_react_dep = i_phase * NUM_REACT_;
         i_react_dep < (i_phase + 1) * NUM_REACT_; ++i_react_dep)
      jac_struct[REACT_(i_react_dep)][ACTIVITY_COEFF_(i_phase)] = true;
    for (int i_prod_dep = i_phase * NUM_PROD_;
         i_prod_dep < (i_phase + 1) * NUM_PROD_; ++i_prod_dep)
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
                                        int **jac_ids, int *rxn_int_data,
                                        double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  for (int i_phase = 0, i_deriv = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    for (int i_react = 0; i_react < NUM_REACT_; i_react++)
      DERIV_ID_(i_deriv++) = deriv_ids[REACT_(i_phase * NUM_REACT_ + i_react)];
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++)
      DERIV_ID_(i_deriv++) = deriv_ids[PROD_(i_phase * NUM_PROD_ + i_prod)];
  }

  // Update the Jacobian ids
  for (int i_phase = 0, i_jac = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = i_phase * NUM_REACT_;
         i_react_ind < (i_phase + 1) * NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][REACT_(i_react_ind)];
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][REACT_(i_react_ind)];
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = i_phase * NUM_PROD_;
         i_prod_ind < (i_phase + 1) * NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][PROD_(i_prod_ind)];
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][PROD_(i_prod_ind)];
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase * NUM_REACT_;
         i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
      JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][WATER_(i_phase)];
    for (int i_prod_dep = i_phase * NUM_PROD_;
         i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
      JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][WATER_(i_phase)];

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase) < 0) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; ++i_react_dep)
        JAC_ID_(i_jac++) = -1;
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; ++i_prod_dep)
        JAC_ID_(i_jac++) = -1;
    } else {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; ++i_react_dep)
        JAC_ID_(i_jac++) =
            jac_ids[REACT_(i_react_dep)][ACTIVITY_COEFF_(i_phase)];
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; ++i_prod_dep)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][ACTIVITY_COEFF_(i_phase)];
    }
  }

  // Calculate a small concentration for aerosol-phase species based on the
  // integration tolerances to use during solving. TODO find a better place
  // to do this
  double *abs_tol = (double *)model_data->abs_tol;
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    SMALL_CONC_(i_phase) = 99999.0;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (SMALL_CONC_(i_phase) >
          abs_tol[REACT_(i_phase * NUM_REACT_ + i_react)] / 100.0)
        SMALL_CONC_(i_phase) =
            abs_tol[REACT_(i_phase * NUM_REACT_ + i_react)] / 100.0;
    }
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (SMALL_CONC_(i_phase) >
          abs_tol[PROD_(i_phase * NUM_PROD_ + i_prod)] / 100.0)
        SMALL_CONC_(i_phase) =
            abs_tol[PROD_(i_phase * NUM_PROD_ + i_prod)] / 100.0;
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
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_aqueous_equilibrium_update_env_state(ModelData *model_data,
                                              int *rxn_int_data,
                                              double *rxn_float_data,
                                              double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate the equilibrium constant
  // (assumes reactant and product concentrations in M)
  double equil_const;
  if (C_ == 0.0) {
    equil_const = A_;
  } else {
    equil_const = A_ * exp(C_ * (1.0 / TEMPERATURE_K_ - 1.0 / 298.0));
  }

  // Set the forward rate constant
  RATE_CONST_FORWARD_ = equil_const * RATE_CONST_REVERSE_;

  return;
}

/** \brief Calculate the reaction rate for a set of conditions using the
 *         standard equation per mixing ratio of water [M_X/s*ug_H2O/m^3]
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param is_water_partial Flag indicating whether the calculation should
 *                         return the partial derivative d_rate/d_H2O
 * \param rate_forward [output] calculated forward rate
 * \param rate_reverse [output] calculated reverse rate
 * \return reaction rate per mixing ratio of water [M_X/s*ug_H2O/m^3]
 */
long double calc_standard_rate(int *rxn_int_data, double *rxn_float_data,
                               double *rxn_env_data, bool is_water_partial,
                               long double *rate_forward,
                               long double *rate_reverse) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  long double react_fact, prod_fact;
  long double water = WATER_CONC_;

  // Get the product of all reactants
  react_fact = (long double)REACT_CONC_(0) * MASS_FRAC_TO_M_(0);
  for (int i_react = 1; i_react < NUM_REACT_; i_react++) {
    react_fact *= REACT_CONC_(i_react) * MASS_FRAC_TO_M_(i_react) / water;
  }

  // Get the product of all product
  prod_fact = (long double)PROD_CONC_(0) * MASS_FRAC_TO_M_(NUM_REACT_);
  prod_fact *= (long double)ACTIVITY_COEFF_VALUE_;
  for (int i_prod = 1; i_prod < NUM_PROD_; i_prod++) {
    prod_fact *=
        PROD_CONC_(i_prod) * MASS_FRAC_TO_M_(NUM_REACT_ + i_prod) / water;
  }

  *rate_forward = RATE_CONST_FORWARD_ * react_fact;
  *rate_reverse = RATE_CONST_REVERSE_ * prod_fact;

  if (is_water_partial) {
    return *rate_forward * (NUM_REACT_ - 1) - *rate_reverse * (NUM_PROD_ - 1);
  } else {
    return *rate_forward - *rate_reverse;
  }
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step of the itegrator (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_aqueous_equilibrium_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase = 0, i_deriv = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // If no aerosol water is present, no reaction occurs
    long double water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Set the concentrations for all species and the activity coefficient
    for (int i_react = 0; i_react < NUM_REACT_; ++i_react)
      REACT_CONC_(i_react) = state[REACT_(i_phase * NUM_REACT_ + i_react)];
    for (int i_prod = 0; i_prod < NUM_PROD_; ++i_prod)
      PROD_CONC_(i_prod) = state[PROD_(i_phase * NUM_PROD_ + i_prod)];
    WATER_CONC_ = state[WATER_(i_phase)];
    if (ACTIVITY_COEFF_(i_phase) >= 0) {
      ACTIVITY_COEFF_VALUE_ = state[ACTIVITY_COEFF_(i_phase)];
    } else {
      ACTIVITY_COEFF_VALUE_ = 1.0;
    }

    // Get the rate using the standard calculation
    long double rate_forward, rate_reverse;
    long double rate =
        calc_standard_rate(rxn_int_data, rxn_float_data, rxn_env_data, false,
                           &rate_forward, &rate_reverse);
    if (rate == ZERO) {
      i_deriv += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv) < 0) {
        i_deriv++;
        continue;
      }
      time_derivative_add_value(time_deriv, DERIV_ID_(i_deriv),
                                -rate_forward / MASS_FRAC_TO_M_(i_react));
      time_derivative_add_value(time_deriv, DERIV_ID_(i_deriv++),
                                rate_reverse / MASS_FRAC_TO_M_(i_react));
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv) < 0) {
        i_deriv++;
        continue;
      }
      time_derivative_add_value(
          time_deriv, DERIV_ID_(i_deriv),
          rate_forward / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod));
      time_derivative_add_value(
          time_deriv, DERIV_ID_(i_deriv++),
          -rate_reverse / MASS_FRAC_TO_M_(NUM_REACT_ + i_prod));
    }
  }

  return;
}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param jac Reaction Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step of the itegrator (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_aqueous_equilibrium_calc_jac_contrib(ModelData *model_data,
                                              Jacobian jac, int *rxn_int_data,
                                              double *rxn_float_data,
                                              double *rxn_env_data,
                                              double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase = 0, i_jac = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // If not aerosol water is present, no reaction occurs
    realtype water = state[WATER_(i_phase)];
    if (water < MIN_WATER_ * SMALL_WATER_CONC_(i_phase)) {
      i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + NUM_PROD_ + 2);
      continue;
    }

    // Calculate the forward rate (M/s)
    realtype forward_rate = RATE_CONST_FORWARD_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      forward_rate *= state[REACT_(i_phase * NUM_REACT_ + i_react)] *
                      MASS_FRAC_TO_M_(i_react) / water;
    }

    // Calculate the reverse rate (M/s)
    realtype reverse_rate = RATE_CONST_REVERSE_;
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      reverse_rate *= state[PROD_(i_phase * NUM_PROD_ + i_prod)] *
                      MASS_FRAC_TO_M_(NUM_REACT_ + i_prod) / water;
    }
    if (ACTIVITY_COEFF_(i_phase) >= 0)
      reverse_rate *= state[ACTIVITY_COEFF_(i_phase)];

    // Add dependence on reactants for reactants and products (forward reaction)
    for (int i_react_ind = 0; i_react_ind < NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
        if (JAC_ID_(i_jac) < 0 || forward_rate == 0.0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(
            jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
            forward_rate / state[REACT_(i_phase * NUM_REACT_ + i_react_ind)] /
                MASS_FRAC_TO_M_(i_react_dep) * water);
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
        if (JAC_ID_(i_jac) < 0 || forward_rate == 0.0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(
            jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_PRODUCTION,
            forward_rate / state[REACT_(i_phase * NUM_REACT_ + i_react_ind)] /
                MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
      }
    }

    // Add dependence on products for reactants and products (reverse reaction)
    for (int i_prod_ind = 0; i_prod_ind < NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
        if (JAC_ID_(i_jac) < 0 || reverse_rate == 0.0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(
            jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_PRODUCTION,
            reverse_rate / state[PROD_(i_phase * NUM_PROD_ + i_prod_ind)] /
                MASS_FRAC_TO_M_(i_react_dep) * water);
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
        if (JAC_ID_(i_jac) < 0 || reverse_rate == 0.0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(
            jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
            reverse_rate / state[PROD_(i_phase * NUM_PROD_ + i_prod_ind)] /
                MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(i_jac), JACOBIAN_LOSS,
          -forward_rate / MASS_FRAC_TO_M_(i_react_dep) * (NUM_REACT_ - 1));
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_PRODUCTION,
          -reverse_rate / MASS_FRAC_TO_M_(i_react_dep) * (NUM_PROD_ - 1));
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac), JACOBIAN_PRODUCTION,
                         -forward_rate /
                             MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) *
                             (NUM_REACT_ - 1));
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
                         -reverse_rate /
                             MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) *
                             (NUM_PROD_ - 1));
    }

    // Add dependence on activity coefficients for reactants and products
    if (ACTIVITY_COEFF_(i_phase) < 0) {
      i_jac += NUM_REACT_ + NUM_PROD_;
      continue;
    }
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac++),
                         JACOBIAN_PRODUCTION,
                         reverse_rate / state[ACTIVITY_COEFF_(i_phase)] /
                             MASS_FRAC_TO_M_(i_react_dep) * water);
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
                         reverse_rate / state[ACTIVITY_COEFF_(i_phase)] /
                             MASS_FRAC_TO_M_(NUM_REACT_ + i_prod_dep) * water);
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
void rxn_aqueous_equilibrium_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nAqueous Equilibrium reaction\n");

  return;
}
