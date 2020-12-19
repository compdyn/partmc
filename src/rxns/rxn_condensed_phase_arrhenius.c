/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Condensed Phase Arrhenius reaction solver functions
 *
 */
/** \file
 * \brief Condensed Phase Arrhenius reaction solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ (int_data[0])
#define NUM_PROD_ (int_data[1])
#define NUM_AERO_PHASE_ (int_data[2])
#define A_ (float_data[0])
#define B_ (float_data[1])
#define C_ (float_data[2])
#define D_ (float_data[3])
#define E_ (float_data[4])
#define RATE_CONSTANT_ (rxn_env_data[0])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 5
#define NUM_ENV_PARAM_ 1
#define REACT_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define PROD_(x) \
  (int_data[NUM_INT_PROP_ + NUM_REACT_ * NUM_AERO_PHASE_ + x] - 1)
#define WATER_(x) \
  (int_data[NUM_INT_PROP_ + (NUM_REACT_ + NUM_PROD_) * NUM_AERO_PHASE_ + x] - 1)
#define DERIV_ID_(x) \
  (int_data[NUM_INT_PROP_ + (NUM_REACT_ + NUM_PROD_ + 1) * NUM_AERO_PHASE_ + x])
#define JAC_ID_(x)          \
  (int_data[NUM_INT_PROP_ + \
            (2 * (NUM_REACT_ + NUM_PROD_) + 1) * NUM_AERO_PHASE_ + x])
#define YIELD_(x) (float_data[NUM_FLOAT_PROP_ + x])
#define KGM3_TO_MOLM3_(x) (float_data[NUM_FLOAT_PROP_ + NUM_PROD_ + x])

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_condensed_phase_arrhenius_get_used_jac_elem(int *rxn_int_data,
                                                     double *rxn_float_data,
                                                     Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Loop over all the instances of the specified phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // Add dependence on reactants for reactants and products
    for (int i_react_ind = i_phase * NUM_REACT_;
         i_react_ind < (i_phase + 1) * NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        jacobian_register_element(jac, REACT_(i_react_dep),
                                  REACT_(i_react_ind));
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        jacobian_register_element(jac, PROD_(i_prod_dep), REACT_(i_react_ind));
    }

    // Add dependence on aerosol-phase water for reactants and products in
    // aqueous reactions
    if (WATER_(i_phase) >= 0) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        jacobian_register_element(jac, REACT_(i_react_dep), WATER_(i_phase));
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        jacobian_register_element(jac, PROD_(i_prod_dep), WATER_(i_phase));
    }
  }

  return;
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_condensed_phase_arrhenius_update_ids(ModelData *model_data,
                                              int *deriv_ids, Jacobian jac,
                                              int *rxn_int_data,
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
    // Add dependence on reactants for reactants and products
    for (int i_react_ind = i_phase * NUM_REACT_;
         i_react_ind < (i_phase + 1) * NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase * NUM_REACT_;
           i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jacobian_get_element_id(jac, REACT_(i_react_dep),
                                                   REACT_(i_react_ind));
      for (int i_prod_dep = i_phase * NUM_PROD_;
           i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jacobian_get_element_id(jac, PROD_(i_prod_dep),
                                                   REACT_(i_react_ind));
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase * NUM_REACT_;
         i_react_dep < (i_phase + 1) * NUM_REACT_; i_react_dep++)
      if (WATER_(i_phase) >= 0) {
        JAC_ID_(i_jac++) =
            jacobian_get_element_id(jac, REACT_(i_react_dep), WATER_(i_phase));
      } else {
        JAC_ID_(i_jac++) = -1;
      }
    for (int i_prod_dep = i_phase * NUM_PROD_;
         i_prod_dep < (i_phase + 1) * NUM_PROD_; i_prod_dep++)
      if (WATER_(i_phase) >= 0) {
        JAC_ID_(i_jac++) =
            jacobian_get_element_id(jac, PROD_(i_prod_dep), WATER_(i_phase));
      } else {
        JAC_ID_(i_jac++) = -1;
      }
  }

  return;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Condensed Phase Arrhenius reaction this only involves recalculating the
 * forward rate constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_condensed_phase_arrhenius_update_env_state(ModelData *model_data,
                                                    int *rxn_int_data,
                                                    double *rxn_float_data,
                                                    double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate the rate constant in (M or mol/m3)
  // k = A*exp(C/T) * (T/D)^B * (1+E*P)
  RATE_CONSTANT_ = A_ * exp(C_ / TEMPERATURE_K_) *
                   (B_ == 0.0 ? 1.0 : pow(TEMPERATURE_K_ / D_, B_)) *
                   (E_ == 0.0 ? 1.0 : (1.0 + E_ * PRESSURE_PA_));

  return;
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step of the itegrator (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_condensed_phase_arrhenius_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase = 0, i_deriv = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    long double unit_conv = 1.0;
    if (WATER_(i_phase) >= 0) {
      unit_conv = state[WATER_(i_phase)];  // convert from kg/m3->L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction
      // occurs
      if (unit_conv <= ZERO) {
        i_deriv += NUM_REACT_ + NUM_PROD_;
        continue;
      }
      unit_conv = 1.0 / unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    long double rate = RATE_CONSTANT_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      rate *= state[REACT_(i_phase * NUM_REACT_ + i_react)] *
              KGM3_TO_MOLM3_(i_react) * unit_conv;
    }

    // Reactant change
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv) < 0) {
        i_deriv++;
        continue;
      }
      time_derivative_add_value(time_deriv, DERIV_ID_(i_deriv++),
                                -rate / (KGM3_TO_MOLM3_(i_react) * unit_conv));
    }

    // Products change
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv) < 0) {
        i_deriv++;
        continue;
      }
      time_derivative_add_value(
          time_deriv, DERIV_ID_(i_deriv++),
          rate * YIELD_(i_prod) /
              (KGM3_TO_MOLM3_(NUM_REACT_ + i_prod) * unit_conv));
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
void rxn_condensed_phase_arrhenius_calc_jac_contrib(
    ModelData *model_data, Jacobian jac, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase = 0, i_jac = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    realtype unit_conv = 1.0;
    if (WATER_(i_phase) >= 0) {
      unit_conv = state[WATER_(i_phase)];  // convert from kg/m3->L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction
      // occurs
      if (unit_conv <= ZERO) {
        i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + 1);
        continue;
      }
      unit_conv = 1.0 / unit_conv;
    }

    // Add dependence on reactants for reactants and products
    for (int i_react_ind = 0; i_react_ind < NUM_REACT_; i_react_ind++) {
      // Calculate d_rate / d_react_i
      realtype rate = RATE_CONSTANT_;
      for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
        if (i_react == i_react_ind) {
          rate *= KGM3_TO_MOLM3_(i_react) * unit_conv;
        } else {
          rate *= state[REACT_(i_phase * NUM_REACT_ + i_react)] *
                  KGM3_TO_MOLM3_(i_react) * unit_conv;
        }
      }

      // Add the Jacobian elements
      //
      // For reactant dependence on reactants
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
        if (JAC_ID_(i_jac) < 0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
                           rate / (KGM3_TO_MOLM3_(i_react_dep) * unit_conv));
      }
      // For product dependence on reactants
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
        if (JAC_ID_(i_jac) < 0) {
          i_jac++;
          continue;
        }
        jacobian_add_value(
            jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_PRODUCTION,
            rate * YIELD_(i_prod_dep) /
                (KGM3_TO_MOLM3_(NUM_REACT_ + i_prod_dep) * unit_conv));
      }
    }

    // Add dependence on aerosol-phase water for reactants and products in
    // aqueous reactions
    if (WATER_(i_phase) < 0) {
      i_jac += NUM_REACT_ + NUM_PROD_;
      continue;
    }

    // Calculate the overall reaction rate (M/s or mol/m3/s)
    realtype rate = RATE_CONSTANT_;
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      rate *= state[REACT_(i_phase * NUM_REACT_ + i_react)] *
              KGM3_TO_MOLM3_(i_react) * unit_conv;
    }

    // Dependence of reactants on aerosol-phase water
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(i_jac++), JACOBIAN_LOSS,
          -(NUM_REACT_ - 1) * rate / KGM3_TO_MOLM3_(i_react_dep));
    }
    // Dependence of products on aerosol-phase water
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac) < 0) {
        i_jac++;
        continue;
      }
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_jac++),
                         JACOBIAN_PRODUCTION,
                         -(NUM_REACT_ - 1) * rate * YIELD_(i_prod_dep) /
                             KGM3_TO_MOLM3_(NUM_REACT_ + i_prod_dep));
    }
  }

  return;
}
#endif

/** \brief Print the Condensed Phase Arrhenius reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_condensed_phase_arrhenius_print(int *rxn_int_data,
                                         double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  int phase_jac_size = (NUM_REACT_ + 1) * (NUM_REACT_ + NUM_PROD_);

  printf("\n\nCondensed Phase Arrhenius reaction\n");

  printf("\n number of reactants:      %d", NUM_REACT_);
  printf("\n number of products:       %d", NUM_PROD_);
  printf("\n number of aerosol phases: %d", NUM_AERO_PHASE_);
  printf("\n A: %le, B: %le, C: %le, D: %le, E: %le", A_, B_, C_, D_, E_);
  printf("\n water state ids (by phase):");
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase)
    printf(" %d", WATER_(i_phase));
  printf("\n *** Reactants ***");
  for (int i_react = 0; i_react < NUM_REACT_; ++i_react) {
    printf("\n reactant %d", i_react);
    printf("\n   kg/m3 -> mol/m3: %le", KGM3_TO_MOLM3_(i_react));
    printf("\n   state id (by phase):");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase)
      printf(" %d", REACT_(i_phase * NUM_REACT_ + i_react));
    printf("\n   deriv id (by phase):");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase)
      printf(" %d", DERIV_ID_(i_phase * (NUM_REACT_ + NUM_PROD_) + i_react));
  }
  printf("\n *** Products ***");
  for (int i_prod = 0; i_prod < NUM_PROD_; ++i_prod) {
    printf("\n product %d", i_prod);
    printf("\n   kg/m3 -> mol/m3: %le", KGM3_TO_MOLM3_(NUM_REACT_ + i_prod));
    printf("\n   yield: %le", YIELD_(i_prod));
    printf("\n   state id (by phase):");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase)
      printf(" %d", PROD_(i_phase * NUM_PROD_ + i_prod));
    printf("\n   deriv id (by phase):");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase)
      printf(" %d", DERIV_ID_(i_phase * (NUM_REACT_ + NUM_PROD_) + NUM_REACT_ +
                              i_prod));
  }
  printf("\n *** Jac Ids (by phase) ***");
  for (int i_ind = 0; i_ind < NUM_REACT_; ++i_ind) {
    for (int i_react = 0; i_react < NUM_REACT_; ++i_react) {
      printf("\n  R->R");
      for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase) {
        printf(" Jac[%d][%d] = %d;", REACT_(i_phase * NUM_REACT_ + i_react),
               REACT_(i_phase * NUM_REACT_ + i_ind),
               JAC_ID_(i_phase * phase_jac_size +
                       i_ind * (NUM_REACT_ + NUM_PROD_) + i_react));
      }
    }
    for (int i_prod = 0; i_prod < NUM_PROD_; ++i_prod) {
      printf("\n  P->R");
      for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase) {
        printf(" Jac[%d][%d] = %d;", PROD_(i_phase * NUM_PROD_ + i_prod),
               REACT_(i_phase * NUM_REACT_ + i_ind),
               JAC_ID_(i_phase * phase_jac_size +
                       i_ind * (NUM_REACT_ + NUM_PROD_) + NUM_REACT_ + i_prod));
      }
    }
  }
  for (int i_react = 0; i_react < NUM_REACT_; ++i_react) {
    printf("\n  R->W");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase) {
      printf(" Jac[%d][%d] = %d;", REACT_(i_phase * NUM_REACT_ + i_react),
             WATER_(i_phase),
             JAC_ID_(i_phase * phase_jac_size +
                     NUM_REACT_ * (NUM_REACT_ + NUM_PROD_) + i_react));
    }
  }
  for (int i_prod = 0; i_prod < NUM_PROD_; ++i_prod) {
    printf("\n  P->W");
    for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; ++i_phase) {
      printf(
          " Jac[%d][%d] = %d;", PROD_(i_phase * NUM_PROD_ + i_prod),
          WATER_(i_phase),
          JAC_ID_(i_phase * phase_jac_size +
                  NUM_REACT_ * (NUM_REACT_ + NUM_PROD_) + NUM_REACT_ + i_prod));
    }
  }
  return;
}
