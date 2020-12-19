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
#include "../util.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

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
#define CONV_ float_data[6]
#define MW_ float_data[7]
#define NUM_AERO_PHASE_ int_data[0]
#define GAS_SPEC_ (int_data[1] - 1)
#define MFP_M_ rxn_env_data[0]
#define ALPHA_ rxn_env_data[1]
#define EQUIL_CONST_ rxn_env_data[2]
#define KGM3_TO_PPM_ rxn_env_data[3]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 8
#define NUM_ENV_PARAM_ 4
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 1 + NUM_AERO_PHASE_ + x]
#define PHASE_INT_LOC_(x) \
  (int_data[NUM_INT_PROP_ + 2 + 6 * NUM_AERO_PHASE_ + x] - 1)
#define PHASE_REAL_LOC_(x) \
  (int_data[NUM_INT_PROP_ + 2 + 7 * NUM_AERO_PHASE_ + x] - 1)
#define AERO_SPEC_(x) (int_data[PHASE_INT_LOC_(x)] - 1)
#define AERO_WATER_(x) (int_data[PHASE_INT_LOC_(x) + 1] - 1)
#define AERO_PHASE_ID_(x) (int_data[PHASE_INT_LOC_(x) + 2] - 1)
#define AERO_REP_ID_(x) (int_data[PHASE_INT_LOC_(x) + 3] - 1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[PHASE_INT_LOC_(x) + 4])
#define PHASE_JAC_ID_(x, s, e) \
  int_data[PHASE_INT_LOC_(x) + 5 + s * NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define SMALL_WATER_CONC_(x) (float_data[PHASE_REAL_LOC_(x)])
#define EFF_RAD_JAC_ELEM_(x, e) float_data[PHASE_REAL_LOC_(x) + 1 + e]
#define NUM_CONC_JAC_ELEM_(x, e) \
  float_data[PHASE_REAL_LOC_(x) + 1 + NUM_AERO_PHASE_JAC_ELEM_(x) + e]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_HL_phase_transfer_get_used_jac_elem(ModelData *model_data,
                                             int *rxn_int_data,
                                             double *rxn_float_data,
                                             Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  bool *aero_jac_elem =
      (bool *)malloc(sizeof(bool) * model_data->n_per_cell_state_var);
  if (aero_jac_elem == NULL) {
    printf(
        "\n\nERROR allocating space for 1D jacobian structure array for HL "
        "partitioning reaction\n\n");
    exit(1);
  }

  jacobian_register_element(jac, GAS_SPEC_, GAS_SPEC_);
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    jacobian_register_element(jac, AERO_SPEC_(i_aero_phase), GAS_SPEC_);
    jacobian_register_element(jac, GAS_SPEC_, AERO_SPEC_(i_aero_phase));
    jacobian_register_element(jac, AERO_SPEC_(i_aero_phase),
                              AERO_SPEC_(i_aero_phase));
    jacobian_register_element(jac, GAS_SPEC_, AERO_WATER_(i_aero_phase));
    jacobian_register_element(jac, AERO_SPEC_(i_aero_phase),
                              AERO_WATER_(i_aero_phase));

    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; i_elem++)
      aero_jac_elem[i_elem] = false;

    int n_jac_elem =
        aero_rep_get_used_jac_elem(model_data, AERO_REP_ID_(i_aero_phase),
                                   AERO_PHASE_ID_(i_aero_phase), aero_jac_elem);
    if (n_jac_elem > NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase)) {
      printf(
          "\n\nERROR Received more Jacobian elements than expected for HL "
          "partitioning reaction. Got %d, expected <= %d",
          n_jac_elem, NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase));
      exit(1);
    }
    int i_used_elem = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; i_elem++) {
      if (aero_jac_elem[i_elem] == true) {
        jacobian_register_element(jac, GAS_SPEC_, i_elem);
        jacobian_register_element(jac, AERO_SPEC_(i_aero_phase), i_elem);
        PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_used_elem) = i_elem;
        PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_used_elem) = i_elem;
        i_used_elem++;
      }
    }
    for (; i_used_elem < NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         ++i_used_elem) {
      PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_used_elem) = -1;
      PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_used_elem) = -1;
    }
    if (i_used_elem != n_jac_elem) {
      printf(
          "\n\nERROR Error setting used Jacobian elements in HL "
          "partitioning reaction %d %d\n\n",
          i_used_elem, n_jac_elem);
      rxn_HL_phase_transfer_print(rxn_int_data, rxn_float_data);
      exit(1);
    }
  }

  free(aero_jac_elem);

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
void rxn_HL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
                                      Jacobian jac, int *rxn_int_data,
                                      double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  DERIV_ID_(0) = deriv_ids[GAS_SPEC_];
  for (int i = 0; i < NUM_AERO_PHASE_; i++)
    DERIV_ID_(i + 1) = deriv_ids[AERO_SPEC_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  JAC_ID_(i_jac++) = jacobian_get_element_id(jac, GAS_SPEC_, GAS_SPEC_);
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    JAC_ID_(i_jac++) =
        jacobian_get_element_id(jac, AERO_SPEC_(i_aero_phase), GAS_SPEC_);
    JAC_ID_(i_jac++) =
        jacobian_get_element_id(jac, GAS_SPEC_, AERO_SPEC_(i_aero_phase));
    JAC_ID_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_(i_aero_phase),
                                               AERO_SPEC_(i_aero_phase));
    JAC_ID_(i_jac++) =
        jacobian_get_element_id(jac, GAS_SPEC_, AERO_WATER_(i_aero_phase));
    JAC_ID_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_(i_aero_phase),
                                               AERO_WATER_(i_aero_phase));
    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         i_elem++) {
      if (PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_elem) > 0) {
        PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_elem) = jacobian_get_element_id(
            jac, GAS_SPEC_, PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_elem));
      }
      if (PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_elem) > 0) {
        PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_elem) = jacobian_get_element_id(
            jac, AERO_SPEC_(i_aero_phase),
            PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_elem));
      }
    }
  }

  // Calculate a small concentration for aerosol-phase water based on the
  // integration tolerances to use during solving. TODO find a better place
  // to do this
  double *abs_tol = (double *)model_data->abs_tol;
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    SMALL_WATER_CONC_(i_aero_phase) = abs_tol[AERO_WATER_(i_aero_phase)] / 10.0;
  }

  return;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Phase Transfer reaction this only involves recalculating the rate
 * constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_HL_phase_transfer_update_env_state(ModelData *model_data,
                                            int *rxn_int_data,
                                            double *rxn_float_data,
                                            double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate the mass accomodation coefficient if the N* parameter
  // was provided, otherwise set it to 0.1 (per Zaveri 2008)
  ALPHA_ = 0.1;
  if (DELTA_H_ != 0.0 || DELTA_S_ != 0.0) {
    double del_G = DELTA_H_ - TEMPERATURE_K_ * DELTA_S_;
    ALPHA_ = exp(-del_G / (UNIV_GAS_CONST_ * TEMPERATURE_K_));
    ALPHA_ = ALPHA_ / (1.0 + ALPHA_);
  }

  // replaced by transition-regime rate equation
#if 0
  // Save c_rms * mass_acc for use in mass transfer rate calc
  MFP_M_ = PRE_C_AVG_ * sqrt(TEMPERATURE_K_) * mass_acc;
#endif

  // save the mean free path [m] for calculating condensation rates
  MFP_M_ = mean_free_path__m(DIFF_COEFF_, TEMPERATURE_K_, ALPHA_);

  // Calculate the Henry's Law equilibrium rate constant in units of
  // (ug_x/ug_H2O/ppm) where x is the aerosol-phase species. (A is in
  // units of M/Pa.)
  if (C_ == 0.0) {
    EQUIL_CONST_ = A_ * PRESSURE_PA_ * 1.0e-6 * MW_;
  } else {
    EQUIL_CONST_ = A_ * PRESSURE_PA_ * 1.0e-6 *
                   exp(C_ * (1.0 / TEMPERATURE_K_ - 1.0 / 298.0)) * MW_;
  }

  // Calculate the conversion from kg/m^3 -> ppm
  KGM3_TO_PPM_ = CONV_ * TEMPERATURE_K_ / PRESSURE_PA_;

  return;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being computed (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_HL_phase_transfer_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius__m(
        model_data,               // model data
        AERO_REP_ID_(i_phase),    // aerosol representation index
        AERO_PHASE_ID_(i_phase),  // aerosol phase index
        &radius,                  // particle effective radius (m)
        NULL);                    // partial derivative

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
        model_data,                // model data
        AERO_REP_ID_(i_phase),     // aerosol representation index
        AERO_PHASE_ID_(i_phase));  // aerosol phase index

    // Get the particle number concentration (#/m3) for per-particle mass
    // concentrations; otherwise set to 1
    realtype number_conc = ONE;
    if (aero_conc_type == 0) {
      aero_rep_get_number_conc__n_m3(
          model_data,               // model data
          AERO_REP_ID_(i_phase),    // aerosol representation index
          AERO_PHASE_ID_(i_phase),  // aerosol phase index
          &number_conc,             // particle number concentration
                                    // (#/m3)
          NULL);                    // partial derivative
    }

    // this was replaced with transition-regime rate equation
#if 0
    long double cond_rate =
        ((long double)1.0) / (radius * radius / (3.0 * DIFF_COEFF_) +
                              4.0 * radius / (3.0 * MFP_M_));
#endif

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    long double cond_rate =
        gas_aerosol_rxn_rate_constant(DIFF_COEFF_, MFP_M_, radius, ALPHA_);

    // Calculate the evaporation rate constant (1/s)
    long double evap_rate = cond_rate / (EQUIL_CONST_);

    // Calculate the evaporation and condensation rates (ppm/s)
    cond_rate *= state[GAS_SPEC_];
    evap_rate *= state[AERO_SPEC_(i_phase)] / state[AERO_WATER_(i_phase)];

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (DERIV_ID_(0) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_(0),
                                number_conc * evap_rate);
      time_derivative_add_value(time_deriv, DERIV_ID_(0),
                                -number_conc * cond_rate);
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (kg/m^3/s)
    if (DERIV_ID_(1 + i_phase) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_(1 + i_phase),
                                -evap_rate / KGM3_TO_PPM_);
      time_derivative_add_value(time_deriv, DERIV_ID_(1 + i_phase),
                                cond_rate / KGM3_TO_PPM_);
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
 * \param time_step Current time step being calculated (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_HL_phase_transfer_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                            int *rxn_int_data,
                                            double *rxn_float_data,
                                            double *rxn_env_data,
                                            realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius__m(
        model_data,                         // model data
        AERO_REP_ID_(i_phase),              // aerosol representation index
        AERO_PHASE_ID_(i_phase),            // aerosol phase index
        &radius,                            // particle effective radius (m)
        &(EFF_RAD_JAC_ELEM_(i_phase, 0)));  // partial derivative

    // Check the aerosol concentration type (per-particle or total per-phase
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
        model_data,                // model data
        AERO_REP_ID_(i_phase),     // aerosol representation index
        AERO_PHASE_ID_(i_phase));  // aerosol phase index

    // Get the particle number concentration (#/m3) for per-particle
    // concentrations
    realtype number_conc = ONE;
    if (aero_conc_type == 0) {
      aero_rep_get_number_conc__n_m3(
          model_data,                          // model data
          AERO_REP_ID_(i_phase),               // aerosol representation index
          AERO_PHASE_ID_(i_phase),             // aerosol phase index
          &number_conc,                        // particle number concentration
                                               // (#/m3)
          &(NUM_CONC_JAC_ELEM_(i_phase, 0)));  // partial derivative
    } else {
      for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_phase); ++i_elem)
        NUM_CONC_JAC_ELEM_(i_phase, i_elem) = ZERO;
    }

    // this was replaced with transition-regime rate equation
#if 0
    long double cond_rate = 1.0 / (radius * radius / (3.0 * DIFF_COEFF_) +
                                   4.0 * radius / (3.0 * MFP_M_));
#endif

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    long double cond_rate =
        gas_aerosol_rxn_rate_constant(DIFF_COEFF_, MFP_M_, radius, ALPHA_);

    // Calculate the evaporation rate constant (1/s)
    long double evap_rate = cond_rate / (EQUIL_CONST_);

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (JAC_ID_(0) >= 0)
      jacobian_add_value(jac, (unsigned int)JAC_ID_(0), JACOBIAN_LOSS,
                         number_conc * cond_rate);
    if (JAC_ID_(1 + i_phase * 5 + 1) >= 0)
      jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 5 + 1),
                         JACOBIAN_PRODUCTION,
                         number_conc * evap_rate / state[AERO_WATER_(i_phase)]);
    if (JAC_ID_(1 + i_phase * 5 + 3) >= 0)
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(1 + i_phase * 5 + 3), JACOBIAN_PRODUCTION,
          -number_conc * evap_rate * state[AERO_SPEC_(i_phase)] /
              state[AERO_WATER_(i_phase)] / state[AERO_WATER_(i_phase)]);

    // Change in the aerosol-phase species is condensation - evaporation
    // (kg/m^3/s)
    if (JAC_ID_(1 + i_phase * 5) >= 0)
      jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 5),
                         JACOBIAN_PRODUCTION, cond_rate / KGM3_TO_PPM_);
    if (JAC_ID_(1 + i_phase * 5 + 2) >= 0)
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(1 + i_phase * 5 + 2), JACOBIAN_LOSS,
          evap_rate / state[AERO_WATER_(i_phase)] / KGM3_TO_PPM_);
    if (JAC_ID_(1 + i_phase * 5 + 4) >= 0)
      jacobian_add_value(
          jac, (unsigned int)JAC_ID_(1 + i_phase * 5 + 4), JACOBIAN_LOSS,
          -evap_rate * state[AERO_SPEC_(i_phase)] / KGM3_TO_PPM_ /
              state[AERO_WATER_(i_phase)] / state[AERO_WATER_(i_phase)]);

    // Calculate the condensation and evaporation rates (ppm/s)
    cond_rate *= state[GAS_SPEC_];
    evap_rate *= state[AERO_SPEC_(i_phase)] / state[AERO_WATER_(i_phase)];

    // Add contributions from species used in aerosol property calculations

    // Calculate d_rate/d_effecive_radius and d_rate/d_number_concentration
    // ( This was replaced with transition-regime rate equation. )
#if 0
    long double d_rate_d_radius =
        -rate * cond_rate *
        (2.0 * radius / (3.0 * DIFF_COEFF_) + 4.0 / (3.0 * MFP_M_));
#endif
    long double d_cond_d_radius = d_gas_aerosol_rxn_rate_constant_d_radius(
                                      DIFF_COEFF_, MFP_M_, radius, ALPHA_) *
                                  state[GAS_SPEC_];
    long double d_evap_d_radius = d_cond_d_radius / state[GAS_SPEC_] /
                                  (EQUIL_CONST_)*state[AERO_SPEC_(i_phase)] /
                                  state[AERO_WATER_(i_phase)];

    // Loop through Jac elements and update
    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_phase); ++i_elem) {
      // Gas-phase species dependencies
      if (PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem) > 0) {
        // species involved in effective radius calculation
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
            JACOBIAN_PRODUCTION,
            number_conc * d_evap_d_radius * EFF_RAD_JAC_ELEM_(i_phase, i_elem));
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
            JACOBIAN_LOSS,
            number_conc * d_cond_d_radius * EFF_RAD_JAC_ELEM_(i_phase, i_elem));

        // species involved in number concentration
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
            JACOBIAN_PRODUCTION,
            evap_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
            JACOBIAN_LOSS, cond_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));
      }

      // Aerosol-phase species dependencies
      if (PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem) > 0) {
        // species involved in effective radius calculation
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
            JACOBIAN_LOSS,
            d_evap_d_radius / KGM3_TO_PPM_ *
                EFF_RAD_JAC_ELEM_(i_phase, i_elem));
        jacobian_add_value(
            jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
            JACOBIAN_PRODUCTION,
            d_cond_d_radius / KGM3_TO_PPM_ *
                EFF_RAD_JAC_ELEM_(i_phase, i_elem));
      }
    }
  }
  return;
}
#endif

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_HL_phase_transfer_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nPhase Transfer reaction\n");

  return;
}
