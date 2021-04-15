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
#include "../util.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Jacobian set indices
#define JAC_GAS 0
#define JAC_AERO 1

// Aerosol mass concentration types
#define PER_PARTICLE_MASS 0
#define TOTAL_PARTICLE_MASS 1

#define DELTA_H_ float_data[0]
#define DELTA_S_ float_data[1]
#define DIFF_COEFF_ float_data[2]
#define PRE_C_AVG_ float_data[3]
#define B1_ float_data[4]
#define B2_ float_data[5]
#define B3_ float_data[6]
#define B4_ float_data[7]
#define CONV_ float_data[8]
#define MW_ float_data[9]
#define NUM_AERO_PHASE_ int_data[0]
#define GAS_SPEC_ (int_data[1] - 1)
#define MFP_M_ rxn_env_data[0]
#define ALPHA_ rxn_env_data[1]
#define EQUIL_CONST_ rxn_env_data[2]
#define KGM3_TO_PPM_ rxn_env_data[3]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 10
#define NUM_ENV_PARAM_ 4
#define AERO_SPEC_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define AERO_ACT_ID_(x) (int_data[NUM_INT_PROP_ + NUM_AERO_PHASE_ + x] - 1)
#define AERO_PHASE_ID_(x) \
  (int_data[NUM_INT_PROP_ + 2 * (NUM_AERO_PHASE_) + x] - 1)
#define AERO_REP_ID_(x) \
  (int_data[NUM_INT_PROP_ + 3 * (NUM_AERO_PHASE_) + x] - 1)
#define DERIV_ID_(x) (int_data[NUM_INT_PROP_ + 4 * (NUM_AERO_PHASE_) + x])
#define GAS_ACT_JAC_ID_(x) \
  int_data[NUM_INT_PROP_ + 1 + 5 * (NUM_AERO_PHASE_) + x]
#define AERO_ACT_JAC_ID_(x) \
  int_data[NUM_INT_PROP_ + 1 + 6 * (NUM_AERO_PHASE_) + x]
#define JAC_ID_(x) (int_data[NUM_INT_PROP_ + 1 + 7 * (NUM_AERO_PHASE_) + x])
#define PHASE_INT_LOC_(x) \
  (int_data[NUM_INT_PROP_ + 2 + 10 * (NUM_AERO_PHASE_) + x] - 1)
#define PHASE_FLOAT_LOC_(x) \
  (int_data[NUM_INT_PROP_ + 2 + 11 * (NUM_AERO_PHASE_) + x] - 1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[PHASE_INT_LOC_(x)])
#define PHASE_JAC_ID_(x, s, e) \
  int_data[PHASE_INT_LOC_(x) + 1 + s * NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define EFF_RAD_JAC_ELEM_(x, e) float_data[PHASE_FLOAT_LOC_(x) + e]
#define NUM_CONC_JAC_ELEM_(x, e) \
  float_data[PHASE_FLOAT_LOC_(x) + NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define MASS_JAC_ELEM_(x, e) \
  float_data[PHASE_FLOAT_LOC_(x) + 2 * NUM_AERO_PHASE_JAC_ELEM_(x) + e]
#define MW_JAC_ELEM_(x, e) \
  float_data[PHASE_FLOAT_LOC_(x) + 3 * NUM_AERO_PHASE_JAC_ELEM_(x) + e]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_SIMPOL_phase_transfer_get_used_jac_elem(ModelData *model_data,
                                                 int *rxn_int_data,
                                                 double *rxn_float_data,
                                                 Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  bool *aero_jac_elem =
      (bool *)malloc(sizeof(bool) * model_data->n_per_cell_state_var);
  if (aero_jac_elem == NULL) {
    printf(
        "\n\nERROR allocating space for 1D Jacobian structure array for "
        "SIMPOL phase transfer reaction\n\n");
    exit(1);
  }

  jacobian_register_element(jac, GAS_SPEC_, GAS_SPEC_);
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    jacobian_register_element(jac, AERO_SPEC_(i_aero_phase), GAS_SPEC_);
    jacobian_register_element(jac, GAS_SPEC_, AERO_SPEC_(i_aero_phase));
    jacobian_register_element(jac, AERO_SPEC_(i_aero_phase),
                              AERO_SPEC_(i_aero_phase));

    if (AERO_ACT_ID_(i_aero_phase) > 0) {
      jacobian_register_element(jac, GAS_SPEC_, AERO_ACT_ID_(i_aero_phase));
      jacobian_register_element(jac, AERO_SPEC_(i_aero_phase),
                                AERO_ACT_ID_(i_aero_phase));
    }

    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem)
      aero_jac_elem[i_elem] = false;

    int n_jac_elem =
        aero_rep_get_used_jac_elem(model_data, AERO_REP_ID_(i_aero_phase),
                                   AERO_PHASE_ID_(i_aero_phase), aero_jac_elem);
    if (n_jac_elem > NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase)) {
      printf(
          "\n\nERROR Received more Jacobian elements than expected for SIMPOL "
          "partitioning reaction. Got %d, expected <= %d",
          n_jac_elem, NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase));
      exit(1);
    }
    int i_used_elem = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem) {
      if (aero_jac_elem[i_elem] == true) {
        jacobian_register_element(jac, GAS_SPEC_, i_elem);
        jacobian_register_element(jac, AERO_SPEC_(i_aero_phase), i_elem);
        PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_used_elem) = i_elem;
        PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_used_elem) = i_elem;
        ++i_used_elem;
      }
    }
    for (; i_used_elem < NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         ++i_used_elem) {
      PHASE_JAC_ID_(i_aero_phase, JAC_GAS, i_used_elem) = -1;
      PHASE_JAC_ID_(i_aero_phase, JAC_AERO, i_used_elem) = -1;
    }
    if (i_used_elem != n_jac_elem) {
      printf(
          "\n\nERROR setting used Jacobian elements in SIMPOL phase "
          "transfer reaction %d %d\n\n",
          i_used_elem, n_jac_elem);
      rxn_SIMPOL_phase_transfer_print(rxn_int_data, rxn_float_data);
      exit(1);
    }
  }

  free(aero_jac_elem);
  return;
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data for finding sub model ids
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_SIMPOL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
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
    if (AERO_ACT_ID_(i_aero_phase) > 0) {
      GAS_ACT_JAC_ID_(i_aero_phase) =
          jacobian_get_element_id(jac, GAS_SPEC_, AERO_ACT_ID_(i_aero_phase));
      AERO_ACT_JAC_ID_(i_aero_phase) = jacobian_get_element_id(
          jac, AERO_SPEC_(i_aero_phase), AERO_ACT_ID_(i_aero_phase));
    } else {
      GAS_ACT_JAC_ID_(i_aero_phase) = -1;
      AERO_ACT_JAC_ID_(i_aero_phase) = -1;
    }
    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_aero_phase);
         ++i_elem) {
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
void rxn_SIMPOL_phase_transfer_update_env_state(ModelData *model_data,
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
  MFP_M_ = PRE_C_AVG_ * sqrt(TEMPERATURE_K_) * mass_acc;  // [m/s]
#endif

  /// save the mean free path [m] for calculating condensation rates
  MFP_M_ = mean_free_path__m(DIFF_COEFF_, TEMPERATURE_K_, ALPHA_);

  // SIMPOL.1 vapor pressure [Pa]
  double vp = B1_ / TEMPERATURE_K_ + B2_ + B3_ * TEMPERATURE_K_ +
              B4_ * log(TEMPERATURE_K_);
  vp = 101325.0 * pow(10, vp);

  // Calculate the conversion from kg_x/m^3 -> ppm_x
  KGM3_TO_PPM_ = CONV_ * TEMPERATURE_K_ / PRESSURE_PA_;

  // Calculate the partitioning coefficient K_eq (ppm_x/kg_x*kg_tot)
  // such that for partitioning species X at equilibrium:
  //   [X]_gas = [X]_aero * activity_coeff_X * K_eq * MW_tot_aero / [tot]_aero
  // where 'tot' indicates all species within an aerosol phase combined
  // with []_gas in (ppm) and []_aero in (kg/m^3)
  EQUIL_CONST_ = vp              // Pa_x / (mol_x_aero/mol_tot_aero)
                 / PRESSURE_PA_  // 1/Pa_air
                 / MW_           // mol_x_aero/kg_x_aero
                 * 1.0e6;        // ppm_x / (Pa_x/P_air)
                                 // = ppm_x * mol_tot_aero / kg_x_aero

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
void rxn_SIMPOL_phase_transfer_calc_deriv_contrib(
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

    // Get the particle number concentration (#/m3)
    realtype number_conc;
    aero_rep_get_number_conc__n_m3(
        model_data,               // model data
        AERO_REP_ID_(i_phase),    // aerosol representation index
        AERO_PHASE_ID_(i_phase),  // aerosol phase index
        &number_conc,             // particle number conc (#/m3)
        NULL);                    // partial derivative

    // Get the total mass of the aerosol phase (kg/m3)
    realtype aero_phase_mass;
    aero_rep_get_aero_phase_mass__kg_m3(
        model_data,               // model data
        AERO_REP_ID_(i_phase),    // aerosol representation index
        AERO_PHASE_ID_(i_phase),  // aerosol phase index
        &aero_phase_mass,         // total aerosol-phase mass (kg/m3)
        NULL);                    // partial derivatives

    // Get the total mass of the aerosol phase (kg/mol)
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_avg_MW__kg_mol(
        model_data,               // model data
        AERO_REP_ID_(i_phase),    // aerosol representation index
        AERO_PHASE_ID_(i_phase),  // aerosol phase index
        &aero_phase_avg_MW,       // avg MW in the aerosol phase (kg/mol)
        NULL);                    // partial derivatives

    // This was replaced with the transition-regime condensation rate
    // equations
#if 0
    long double cond_rate =
        ((long double)1.0) / (radius * radius / (3.0 * DIFF_COEFF_) +
                              4.0 * radius / (3.0 * MFP_M_));
#endif

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (m3/#/s)
    long double cond_rate =
        gas_aerosol_rxn_rate_constant(DIFF_COEFF_, MFP_M_, radius, ALPHA_);

    // Calculate the evaporation rate constant (ppm_x*m^3/kg_x/s)
    long double evap_rate =
        cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW / aero_phase_mass);

    // Get the activity coefficient (if one exists)
    long double act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase) > -1) {
      act_coeff = state[AERO_ACT_ID_(i_phase)];
    }

    // Calculate aerosol-phase evaporation rate (ppm/s)
    evap_rate *= act_coeff;

    // Calculate the evaporation and condensation rates
    cond_rate *= state[GAS_SPEC_];
    evap_rate *= state[AERO_SPEC_(i_phase)];

    // per-particle mass concentration rates
    if (aero_conc_type == PER_PARTICLE_MASS) {
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

      // total-aerosol mass concentration rates
    } else {
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
                                  -number_conc * evap_rate / KGM3_TO_PPM_);
        time_derivative_add_value(time_deriv, DERIV_ID_(1 + i_phase),
                                  number_conc * cond_rate / KGM3_TO_PPM_);
      }
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
void rxn_SIMPOL_phase_transfer_calc_jac_contrib(ModelData *model_data,
                                                Jacobian jac, int *rxn_int_data,
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

    // Get the particle number concentration (#/m3)
    realtype number_conc;
    aero_rep_get_number_conc__n_m3(
        model_data,                          // model data
        AERO_REP_ID_(i_phase),               // aerosol representation index
        AERO_PHASE_ID_(i_phase),             // aerosol phase index
        &number_conc,                        // particle number conc (#/m3)
        &(NUM_CONC_JAC_ELEM_(i_phase, 0)));  // partial derivative

    // Get the total mass of the aerosol phase (kg/m3)
    realtype aero_phase_mass;
    aero_rep_get_aero_phase_mass__kg_m3(
        model_data,                      // model data
        AERO_REP_ID_(i_phase),           // aerosol representation index
        AERO_PHASE_ID_(i_phase),         // aerosol phase index
        &aero_phase_mass,                // total aerosol-phase mass (kg/m3)
        &(MASS_JAC_ELEM_(i_phase, 0)));  // partial derivatives

    // Get the total average MW of the aerosol phase (kg/mol)
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_avg_MW__kg_mol(
        model_data,                    // model data
        AERO_REP_ID_(i_phase),         // aerosol representation index
        AERO_PHASE_ID_(i_phase),       // aerosol phase index
        &aero_phase_avg_MW,            // avg MW in the aerosol phase (kg/mol)
        &(MW_JAC_ELEM_(i_phase, 0)));  // partial derivatives

    // This was replaced with the transition-regime condensation rate
    // equations
#if 0
    long double cond_rate =
        ((long double)1.0) / (radius * radius / (3.0 * DIFF_COEFF_) +
                              4.0 * radius / (3.0 * MFP_M_));
#endif

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (m3/#/s)
    long double cond_rate =
        gas_aerosol_rxn_rate_constant(DIFF_COEFF_, MFP_M_, radius, ALPHA_);

    // Calculate the evaporation rate constant (ppm_x*m^3/kg_x/s)
    long double evap_rate =
        cond_rate * (EQUIL_CONST_ * aero_phase_avg_MW / aero_phase_mass);

    // Get the activity coefficient (if one exists)
    long double act_coeff = 1.0;
    if (AERO_ACT_ID_(i_phase) > -1) {
      act_coeff = state[AERO_ACT_ID_(i_phase)];
    }

    // per-particle mass concentrations
    if (aero_conc_type == PER_PARTICLE_MASS) {
      // Change in the gas-phase is evaporation - condensation (ppm/s)
      if (JAC_ID_(1 + i_phase * 3 + 1) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3 + 1),
                           JACOBIAN_PRODUCTION,
                           number_conc * evap_rate * act_coeff);
      }
      if (JAC_ID_(0) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(0), JACOBIAN_LOSS,
                           number_conc * cond_rate);
      }

      // Change in the aerosol-phase species is condensation - evaporation
      // (kg/m^3/s)
      if (JAC_ID_(1 + i_phase * 3) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3),
                           JACOBIAN_PRODUCTION, cond_rate / KGM3_TO_PPM_);
      }
      if (JAC_ID_(1 + i_phase * 3 + 2) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3 + 2),
                           JACOBIAN_LOSS, evap_rate * act_coeff / KGM3_TO_PPM_);
      }

      // Activity coefficient contributions
      if (GAS_ACT_JAC_ID_(i_phase) > 0) {
        jacobian_add_value(
            jac, (unsigned int)GAS_ACT_JAC_ID_(i_phase), JACOBIAN_PRODUCTION,
            number_conc * evap_rate * state[AERO_SPEC_(i_phase)]);
      }
      if (AERO_ACT_JAC_ID_(i_phase) > 0) {
        jacobian_add_value(
            jac, (unsigned int)AERO_ACT_JAC_ID_(i_phase), JACOBIAN_LOSS,
            evap_rate / KGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)]);
      }

      // total-particle mass concentrations
    } else {
      // Change in the gas-phase is evaporation - condensation (ppm/s)
      if (JAC_ID_(1 + i_phase * 3 + 1) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3 + 1),
                           JACOBIAN_PRODUCTION,
                           number_conc * evap_rate * act_coeff);
      }
      if (JAC_ID_(0) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(0), JACOBIAN_LOSS,
                           number_conc * cond_rate);
      }

      // Change in the aerosol-phase species is condensation - evaporation
      // (kg/m^3/s)
      if (JAC_ID_(1 + i_phase * 3) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3),
                           JACOBIAN_PRODUCTION,
                           number_conc * cond_rate / KGM3_TO_PPM_);
      }
      if (JAC_ID_(1 + i_phase * 3 + 2) >= 0) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(1 + i_phase * 3 + 2),
                           JACOBIAN_LOSS,
                           number_conc * evap_rate * act_coeff / KGM3_TO_PPM_);
      }

      // Activity coefficient contributions
      if (GAS_ACT_JAC_ID_(i_phase) > 0) {
        jacobian_add_value(
            jac, (unsigned int)GAS_ACT_JAC_ID_(i_phase), JACOBIAN_PRODUCTION,
            number_conc * evap_rate * state[AERO_SPEC_(i_phase)]);
      }
      if (AERO_ACT_JAC_ID_(i_phase) > 0) {
        jacobian_add_value(jac, (unsigned int)AERO_ACT_JAC_ID_(i_phase),
                           JACOBIAN_LOSS,
                           number_conc * evap_rate / KGM3_TO_PPM_ *
                               state[AERO_SPEC_(i_phase)]);
      }
    }

    // Get the overall rates
    evap_rate *= act_coeff;
    cond_rate *= state[GAS_SPEC_];
    evap_rate *= state[AERO_SPEC_(i_phase)];

    // Calculate partial derivatives

    // this was replaced with the transition regime rate equations
#if 0
    realtype d_cond_d_radius =
        -(2.0 * radius / (3.0 * DIFF_COEFF_) + 4.0 / (3.0 * MFP_M_)) *
        cond_rate * cond_rate / state[GAS_SPEC_];
#endif
    realtype d_cond_d_radius = d_gas_aerosol_rxn_rate_constant_d_radius(
                                   DIFF_COEFF_, MFP_M_, radius, ALPHA_) *
                               state[GAS_SPEC_];
    realtype d_evap_d_radius = d_cond_d_radius / state[GAS_SPEC_] *
                               EQUIL_CONST_ * aero_phase_avg_MW /
                               aero_phase_mass * state[AERO_SPEC_(i_phase)];
    realtype d_evap_d_mass = -evap_rate / aero_phase_mass;
    realtype d_evap_d_MW = evap_rate / aero_phase_avg_MW;

    // per-particle mass concentrations
    if (aero_conc_type == PER_PARTICLE_MASS) {
      // Loop through Jac elements and update
      for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_phase);
           ++i_elem) {
        // Gas-phase species dependencies
        if (PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem) > 0) {
          // species involved in effective radius calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_radius *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_LOSS,
              number_conc * d_cond_d_radius *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));

          // species involved in number concentration
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              evap_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_LOSS, cond_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));

          // species involved in mass calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_mass * MASS_JAC_ELEM_(i_phase, i_elem));

          // species involved in average MW calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_MW * MW_JAC_ELEM_(i_phase, i_elem));
        }

        // Aerosol-phase species dependencies
        if (PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem) > 0) {
          // species involved in effective radius calculations
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

          // species involved in mass calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              d_evap_d_mass / KGM3_TO_PPM_ * MASS_JAC_ELEM_(i_phase, i_elem));

          // species involved in average MW calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              d_evap_d_MW / KGM3_TO_PPM_ * MW_JAC_ELEM_(i_phase, i_elem));
        }
      }

      // total-particle mass concentrations
    } else {
      // Loop through Jac elements and update
      for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_(i_phase);
           ++i_elem) {
        // Gas-phase species dependencies
        if (PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem) > 0) {
          // species involved in effective radius calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_radius *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_LOSS,
              number_conc * d_cond_d_radius *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));

          // species involved in number concentration
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              evap_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_LOSS, cond_rate * NUM_CONC_JAC_ELEM_(i_phase, i_elem));

          // species involved in mass calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_mass * MASS_JAC_ELEM_(i_phase, i_elem));

          // species involved in average MW calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_GAS, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_evap_d_MW * MW_JAC_ELEM_(i_phase, i_elem));
        }

        // Aerosol-phase species dependencies
        if (PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem) > 0) {
          // species involved in effective radius calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              number_conc * d_evap_d_radius / KGM3_TO_PPM_ *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_PRODUCTION,
              number_conc * d_cond_d_radius / KGM3_TO_PPM_ *
                  EFF_RAD_JAC_ELEM_(i_phase, i_elem));

          // species involved in number concentration
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              evap_rate / KGM3_TO_PPM_ * NUM_CONC_JAC_ELEM_(i_phase, i_elem));
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_PRODUCTION,
              cond_rate / KGM3_TO_PPM_ * NUM_CONC_JAC_ELEM_(i_phase, i_elem));

          // species involved in mass calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              number_conc * d_evap_d_mass / KGM3_TO_PPM_ *
                  MASS_JAC_ELEM_(i_phase, i_elem));

          // species involved in average MW calculations
          jacobian_add_value(
              jac, (unsigned int)PHASE_JAC_ID_(i_phase, JAC_AERO, i_elem),
              JACOBIAN_LOSS,
              number_conc * d_evap_d_MW / KGM3_TO_PPM_ *
                  MW_JAC_ELEM_(i_phase, i_elem));
        }
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
void rxn_SIMPOL_phase_transfer_print(int *rxn_int_data,
                                     double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nSIMPOL.1 Phase Transfer reaction\n");
  printf("\ndelta H: %le delta S: %le diffusion coeff: %le Pre C_avg: %le",
         DELTA_H_, DELTA_S_, DIFF_COEFF_, PRE_C_AVG_);
  printf("\nB1: %le B2: %le B3:%le B4: %le", B1_, B2_, B3_, B4_);
  printf("\nconv: %le MW: %le", CONV_, MW_);
  printf("\nNumber of aerosol phases: %d", NUM_AERO_PHASE_);
  printf("\nGas-phase species id: %d", GAS_SPEC_);
  printf("\nGas-phase derivative id: %d", DERIV_ID_(0));
  printf("\ndGas/dGas Jac id: %d", JAC_ID_(0));
  printf("\n*** Aerosol phase data ***");
  for (int i = 0; i < NUM_AERO_PHASE_; ++i) {
    printf("\n  Aerosol species id: %d", AERO_SPEC_(i));
    printf("\n  Activity coefficient id: %d", AERO_ACT_ID_(i));
    printf("\n  Aerosol phase id: %d", AERO_PHASE_ID_(i));
    printf("\n  Aerosol representation id: %d", AERO_REP_ID_(i));
    printf("\n  Aerosol species derivative id: %d", DERIV_ID_(i + 1));
    printf("\n  dGas/dAct coeff  Jac id: %d", GAS_ACT_JAC_ID_(i));
    printf("\n  dAero/dAct coeff Jac id: %d", AERO_ACT_JAC_ID_(i));
    printf("\n  dAero/dGas  Jac id: %d", JAC_ID_(1 + i * 3));
    printf("\n  dGas/dAero  Jac id: %d", JAC_ID_(2 + i * 3));
    printf("\n  dAero/dAero Jac id: %d", JAC_ID_(3 + i * 3));
    printf("\n  Number of aerosol-phase species Jac elements: %d",
           NUM_AERO_PHASE_JAC_ELEM_(i));
    printf("\n  dGas/dx ids:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %d", PHASE_JAC_ID_(i, JAC_GAS, j));
    printf("\n  dAero/dx ids:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %d", PHASE_JAC_ID_(i, JAC_AERO, j));
    printf("\n  Effective radius Jac elem:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %le", EFF_RAD_JAC_ELEM_(i, j));
    printf("\n  Number concentration Jac elem:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %le", NUM_CONC_JAC_ELEM_(i, j));
    printf("\n  Aerosol mass Jac elem:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %le", MASS_JAC_ELEM_(i, j));
    printf("\n  Average MW Jac elem:");
    for (int j = 0; j < NUM_AERO_PHASE_JAC_ELEM_(i); ++j)
      printf(" %le", MW_JAC_ELEM_(i, j));
  }

  return;
}
