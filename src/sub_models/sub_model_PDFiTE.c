/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * PDFiTE Activity sub model solver functions
 *
 */
/** \file
 * \brief PDFiTE Activity sub model solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../sub_models.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_PHASE_ (int_data[0])
#define GAS_WATER_ID_ (int_data[1] - 1)
#define NUM_ION_PAIRS_ (int_data[2])
#define INT_DATA_SIZE_ (int_data[3])
#define FLOAT_DATA_SIZE_ (int_data[4])
#define PPM_TO_RH_ (sub_model_env_data[0])
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 0
#define NUM_ENV_PARAM_ 1
#define PHASE_ID_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define PAIR_INT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_ + NUM_PHASE_ + x] - 1)
#define PAIR_FLOAT_PARAM_LOC_(x) \
  (int_data[NUM_INT_PROP_ + NUM_PHASE_ + NUM_ION_PAIRS_ + x] - 1)
#define ION_PAIR_ACT_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)])
#define NUM_CATION_(x) (int_data[PAIR_INT_PARAM_LOC_(x) + 1])
#define NUM_ANION_(x) (int_data[PAIR_INT_PARAM_LOC_(x) + 2])
#define CATION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x) + 3])
#define ANION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x) + 4])
#define NUM_INTER_(x) (int_data[PAIR_INT_PARAM_LOC_(x) + 5])
#define JAC_WATER_ID_(p, x) (int_data[PAIR_INT_PARAM_LOC_(x) + 6 + p])
#define JAC_CATION_ID_(p, x, y) \
  (int_data[PAIR_INT_PARAM_LOC_(x) + 6 + NUM_PHASE_ + p * NUM_ION_PAIRS_ + y])
#define JAC_ANION_ID_(p, x, y)                                               \
  (int_data[PAIR_INT_PARAM_LOC_(x) + 6 + (1 + NUM_ION_PAIRS_) * NUM_PHASE_ + \
            p * NUM_ION_PAIRS_ + y])
#define NUM_B_(x, y)                     \
  (int_data[PAIR_INT_PARAM_LOC_(x) + 6 + \
            (1 + 2 * NUM_ION_PAIRS_) * NUM_PHASE_ + y])
#define INTER_SPEC_ID_(x, y)                                             \
  (int_data[PAIR_INT_PARAM_LOC_(x) + 6 +                                 \
            (1 + 2 * NUM_ION_PAIRS_) * NUM_PHASE_ + NUM_INTER_(x) + y] - \
   1)
#define INTER_SPEC_LOC_(x, y)                                                  \
  (int_data[PAIR_INT_PARAM_LOC_(x) + 6 +                                       \
            (1 + 2 * NUM_ION_PAIRS_) * NUM_PHASE_ + 2 * (NUM_INTER_(x)) + y] - \
   1)
#define CATION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define ANION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x) + 1])
#define CATION_N_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x) + 2])
#define ANION_N_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x) + 3])
#define MIN_RH_(x, y) (float_data[INTER_SPEC_LOC_(x, y)])
#define MAX_RH_(x, y) (float_data[INTER_SPEC_LOC_(x, y) + 1])
#define B_Z_(x, y, z) (float_data[INTER_SPEC_LOC_(x, y) + 2 + z])

/** \brief Flag Jacobian elements used by this sub model
 *
 * activity sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_struct A matrix of flags for needed Jac elements
 */
void sub_model_PDFiTE_get_used_jac_elem(int *sub_model_int_data,
                                        double *sub_model_float_data,
                                        bool **jac_struct) {
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; ++i_ion_pair) {
      jac_struct[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                [GAS_WATER_ID_] = true;
      for (int j_ion_pair = 0; j_ion_pair < NUM_ION_PAIRS_; ++j_ion_pair) {
        jac_struct[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                  [PHASE_ID_(i_phase) + CATION_ID_(j_ion_pair)] = true;
        jac_struct[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                  [PHASE_ID_(i_phase) + ANION_ID_(j_ion_pair)] = true;
      }
    }
  }
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * activity sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param deriv_ids Indices for state array variables on the solver state array
 * \param jac_ids Indices for Jacobian elements in the sparse data array
 */
void sub_model_PDFiTE_update_ids(int *sub_model_int_data,
                                 double *sub_model_float_data, int *deriv_ids,
                                 int **jac_ids) {
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; ++i_ion_pair) {
      JAC_WATER_ID_(i_phase, i_ion_pair) =
          jac_ids[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                 [GAS_WATER_ID_];
      for (int j_ion_pair = 0; j_ion_pair < NUM_ION_PAIRS_; ++j_ion_pair) {
        JAC_CATION_ID_(i_phase, i_ion_pair, j_ion_pair) =
            jac_ids[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                   [PHASE_ID_(i_phase) + CATION_ID_(j_ion_pair)];
        JAC_ANION_ID_(i_phase, i_ion_pair, j_ion_pair) =
            jac_ids[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)]
                   [PHASE_ID_(i_phase) + ANION_ID_(j_ion_pair)];
      }
    }
  }
}

/** \brief Update sub model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 */
void sub_model_PDFiTE_update_env_state(int *sub_model_int_data,
                                       double *sub_model_float_data,
                                       double *sub_model_env_data,
                                       ModelData *model_data) {
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate PPM_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  double t_steam = 373.15;  // steam temperature (K)
  double a = 1.0 - t_steam / TEMPERATURE_K_;

  a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a;
  double water_vp = 101325.0 * exp(a);  // (Pa)

  PPM_TO_RH_ = PRESSURE_PA_ / water_vp / 1.0e6;  // (1/ppm)

  return;
}

/** \brief Perform the sub-model calculations for the current model state
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data including the current state and
 *                   environmental conditions
 */
void sub_model_PDFiTE_calculate(int *sub_model_int_data,
                                double *sub_model_float_data,
                                double *sub_model_env_data,
                                ModelData *model_data) {
  double *state = model_data->grid_cell_state;
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate the water activity---i.e., relative humidity (0-1)
  long double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];

  // Keep a_w within 0-1
  // TODO Filter =( try to remove
  if (a_w < 0.0) a_w = 0.0;
  if (a_w > 1.0) a_w = 1.0;

  // Calculate ion_pair activity coefficients in each phase
  for (int i_phase = 0; i_phase < NUM_PHASE_; i_phase++) {
    // Initialize omega' (defined below)
    long double omega_prime = 0.0;

    // Calculate the number of moles of each ion
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; i_ion_pair++) {
      // N (mol_i/m3) = c_i (ug/m3) / 10^6 (ug/g) / MW_i (ug/umol)
      CATION_N_(i_ion_pair) =
          state[PHASE_ID_(i_phase) + CATION_ID_(i_ion_pair)] /
          CATION_MW_(i_ion_pair) / 1000000.0;
      ANION_N_(i_ion_pair) = state[PHASE_ID_(i_phase) + ANION_ID_(i_ion_pair)] /
                             ANION_MW_(i_ion_pair) / 1000000.0;

    }  // Loop on primary ion_pair

    // Calculate the activity coefficient
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; i_ion_pair++) {
      // If there are no interactions, the remaining ion pairs will not
      // have activity calculations (they only participate in interactions)
      if (NUM_INTER_(i_ion_pair) == 0) break;

      // Calculate omega for this ion_pair
      // (eq. 15 in \cite{Topping2009}) as the sum of
      //   v_e^C = 2 * ( v_cation + v_anion) * N_cation * N_anion
      // across all other ion pairs (not i_ion_pair)
      // where v_x is the stoichiometric coefficient for species x in
      // the other ion_pair and N_x is its concentration.
      long double omega = 0.0;
      for (int j_ion_pair = 0; j_ion_pair < NUM_ION_PAIRS_; ++j_ion_pair) {
        if (i_ion_pair == j_ion_pair) continue;
        omega += (long double)2.0 *
                 (NUM_CATION_(j_ion_pair) + NUM_ANION_(j_ion_pair)) *
                 CATION_N_(j_ion_pair) * ANION_N_(j_ion_pair);
      }

      // Initialize ln(gamma)
      long double ln_gamma = 0.0;

      // Add contributions from each interacting ion_pair
      for (int i_inter = 0; i_inter < NUM_INTER_(i_ion_pair); i_inter++) {
        // Only include interactions in the correct RH range
        // where the range is in (minRH, maxRH] except when a_w = 0.0
        // where the range is in [0.0, maxRH]
        if ((a_w <= MIN_RH_(i_ion_pair, i_inter) ||
             a_w > MAX_RH_(i_ion_pair, i_inter)) &&
            !(a_w <= 0.0 && MIN_RH_(i_ion_pair, i_inter) <= 0.0))
          continue;

        // Get the ion_pair id of the interacting species
        int j_ion_pair = INTER_SPEC_ID_(i_ion_pair, i_inter);

        // Calculate ln_gamma_inter
        long double ln_gamma_inter = 0.0;
        for (int i_B = 0; i_B < NUM_B_(i_ion_pair, i_inter); i_B++) {
          ln_gamma_inter += B_Z_(i_ion_pair, i_inter, i_B) * pow(a_w, i_B);
        }

        // If this is the "self" interaction, ln_gamma_inter is ln(gamma_0A)
        // (eq. 15 in \cite{Topping2009})
        if (i_ion_pair == j_ion_pair) {
          // Add contribution to ln(gamma_A) from ln(gamma_0A)
          ln_gamma += ln_gamma_inter;
        }
        // ... otherwise it is d(ln(gamma_A))/g(N_B,M N_B,x)
        // (eq. 15 in \cite{Topping2009})
        else {
          // Add contribution to ln(gamma_A) from interacting ion_pair.
          // When omega == 0, N_cation or N_anion must be 0, so
          // skip to avoid divide by zero errors
          if (omega > 0.0) {
            ln_gamma += ln_gamma_inter * CATION_N_(j_ion_pair) *
                        ANION_N_(j_ion_pair) / omega;
          }
        }

      }  // Loop on interacting ion_pairs

      // Set the ion_pair activity
      state[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)] = exp(ln_gamma);

    }  // Loop on primary ion_pairs

  }  // Loop on aerosol phases
}

/** \brief Add contributions to the Jacobian from derivates calculated using the
 * output of this sub model
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current time step in [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_PDFiTE_get_jac_contrib(int *sub_model_int_data,
                                      double *sub_model_float_data,
                                      double *sub_model_env_data,
                                      ModelData *model_data, realtype *J,
                                      double time_step) {
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate the water activity---i.e., relative humidity (0-1)
  double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];

  // Keep a_w within 0-1
  // TODO Filter =( try to remove
  if (a_w < 0.0) a_w = 0.0;
  if (a_w > 1.0) a_w = 1.0;

  // Calculate ion_pair activity coefficients in each phase
  for (int i_phase = 0; i_phase < NUM_PHASE_; i_phase++) {
    // Initialize omega' (defined below)
    long double omega_prime = 0.0;

    // Calculate the number of moles of each ion
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; i_ion_pair++) {
      // N (mol_i/m3) = c_i (ug/m3) / 10^6 (ug/g) / MW_i (ug/umol)
      CATION_N_(i_ion_pair) =
          state[PHASE_ID_(i_phase) + CATION_ID_(i_ion_pair)] /
          CATION_MW_(i_ion_pair) / 1000000.0;
      ANION_N_(i_ion_pair) = state[PHASE_ID_(i_phase) + ANION_ID_(i_ion_pair)] /
                             ANION_MW_(i_ion_pair) / 1000000.0;

    }  // Loop on primary ion_pair

    // Calculate the activity coefficient
    for (int i_ion_pair = 0; i_ion_pair < NUM_ION_PAIRS_; i_ion_pair++) {
      // If there are no interactions, the remaining ion pairs will not
      // have activity calculations (they only participate in interactions)
      if (NUM_INTER_(i_ion_pair) == 0) break;

      // Calculate omega for this ion_pair
      // (eq. 15 in \cite{Topping2009}) as the sum of
      //   v_e^C = 2 * ( v_cation + v_anion) * N_cation * N_anion
      // across all other ion pairs (not i_ion_pair)
      // where v_x is the stoichiometric coefficient for species x in
      // the other ion_pair and N_x is its concentration.
      long double omega = 0.0;
      for (int j_ion_pair = 0; j_ion_pair < NUM_ION_PAIRS_; ++j_ion_pair) {
        if (i_ion_pair == j_ion_pair) continue;
        omega += (long double)2.0 *
                 (NUM_CATION_(j_ion_pair) + NUM_ANION_(j_ion_pair)) *
                 CATION_N_(j_ion_pair) * ANION_N_(j_ion_pair);
      }

      // Initialize ln(gamma)
      long double ln_gamma = 0.0;

      // Add contributions from each interacting ion_pair
      for (int i_inter = 0; i_inter < NUM_INTER_(i_ion_pair); i_inter++) {
        // Only include interactions in the correct RH range
        // where the range is in (minRH, maxRH] except when a_w = 0.0
        // where the range is in [0.0, maxRH]
        if ((a_w <= MIN_RH_(i_ion_pair, i_inter) ||
             a_w > MAX_RH_(i_ion_pair, i_inter)) &&
            !(a_w <= 0.0 && MIN_RH_(i_ion_pair, i_inter) <= 0.0))
          continue;

        // Get the ion_pair id of the interacting species
        int j_ion_pair = INTER_SPEC_ID_(i_ion_pair, i_inter);

        // Calculate ln_gamma_inter
        long double ln_gamma_inter = 0.0;
        for (int i_B = 0; i_B < NUM_B_(i_ion_pair, i_inter); i_B++) {
          ln_gamma_inter += B_Z_(i_ion_pair, i_inter, i_B) * pow(a_w, i_B);
        }

        // If this is the "self" interaction, ln_gamma_inter is ln(gamma_0A)
        // (eq. 15 in \cite{Topping2009})
        if (i_ion_pair == j_ion_pair) {
          // Add contribution to ln(gamma_A) from ln(gamma_0A)
          ln_gamma += ln_gamma_inter;
        }
        // ... otherwise it is d(ln(gamma_A))/g(N_B,M N_B,x)
        // (eq. 15 in \cite{Topping2009})
        else {
          // Add contribution to ln(gamma_A) from interacting ion_pair.
          // When omega == 0, N_cation or N_anion must be 0, so
          // skip to avoid divide by zero errors
          if (omega > 0.0) {
            ln_gamma += ln_gamma_inter * CATION_N_(j_ion_pair) *
                        ANION_N_(j_ion_pair) / omega;
          }
        }

      }  // Loop on interacting ion_pairs

      long double gamma_i = exp(ln_gamma);

      // Loop through the ion pairs to set the partial derivatives
      for (int i_inter = 0; i_inter < NUM_INTER_(i_ion_pair); i_inter++) {
        // Only include interactions in the correct RH range
        // where the range is in (minRH, maxRH] except when a_w = 0.0
        // where the range is in [0.0, maxRH]
        if ((a_w <= MIN_RH_(i_ion_pair, i_inter) ||
             a_w > MAX_RH_(i_ion_pair, i_inter)) &&
            !(a_w <= 0.0 && MIN_RH_(i_ion_pair, i_inter) <= 0.0))
          continue;

        // Get the ion_pair id of the interacting species
        int j_ion_pair = INTER_SPEC_ID_(i_ion_pair, i_inter);

        // Calculate ln_gamma_inter and dln_gamma_inter_d_water
        long double ln_gamma_inter = B_Z_(i_ion_pair, i_inter, 0);
        long double d_ln_gamma_inter_d_water = 0.0;
        for (int i_B = 1; i_B < NUM_B_(i_ion_pair, i_inter); i_B++) {
          ln_gamma_inter += B_Z_(i_ion_pair, i_inter, i_B) * pow(a_w, i_B);
          d_ln_gamma_inter_d_water +=
              B_Z_(i_ion_pair, i_inter, i_B) * i_B * pow(a_w, i_B - 1);
        }
        d_ln_gamma_inter_d_water *= PPM_TO_RH_;

        // If this is the "self" interaction, ln_gamma_inter is ln(gamma_0A)
        // (eq. 15 in \cite{Topping2009})
        if (i_ion_pair == j_ion_pair) {
          // Add the water contribution to ln(gamma_0A)
          J[JAC_WATER_ID_(i_phase, i_ion_pair)] +=
              gamma_i * d_ln_gamma_inter_d_water;

        }
        // ... otherwise it is d(ln(gamma_A))/g(N_B,M N_B,x)
        // (eq. 15 in \cite{Topping2009})
        else {
          // Add contribution to ln(gamma_A) from interacting ion_pair.
          // When omega == 0, N_cation or N_anion must be 0, so
          // skip to avoid divide by zero errors
          if (omega > 0.0) {
            // d_gamma / d_cation
            J[JAC_CATION_ID_(i_phase, i_ion_pair, j_ion_pair)] +=
                gamma_i * ln_gamma_inter * ANION_N_(j_ion_pair) /
                CATION_MW_(j_ion_pair) / omega * 1.0e-6;

            // d_gamma / d_anion
            J[JAC_ANION_ID_(i_phase, i_ion_pair, j_ion_pair)] +=
                gamma_i * ln_gamma_inter * CATION_N_(j_ion_pair) /
                ANION_MW_(j_ion_pair) / omega * 1.0e-6;

            // d_gamma / d_water
            J[JAC_WATER_ID_(i_phase, i_ion_pair)] +=
                gamma_i * CATION_N_(j_ion_pair) * ANION_N_(j_ion_pair) / omega *
                d_ln_gamma_inter_d_water;

            // Add the contributions from other ion pairs to omega
            for (int k_ion_pair = 0; k_ion_pair < NUM_ION_PAIRS_;
                 ++k_ion_pair) {
              // The ion pair whose activity is being calculated is not
              // included in omega
              if (k_ion_pair == i_ion_pair) continue;

              // d_gamma / d_cation
              J[JAC_CATION_ID_(i_phase, i_ion_pair, k_ion_pair)] +=
                  -gamma_i * ln_gamma_inter * CATION_N_(j_ion_pair) *
                  ANION_N_(j_ion_pair) / (omega * omega) * 2.0 *
                  (NUM_CATION_(k_ion_pair) + NUM_ANION_(k_ion_pair)) *
                  ANION_N_(k_ion_pair) / CATION_MW_(k_ion_pair) * 1.0e-6;

              // d_gamma / d_anion
              J[JAC_ANION_ID_(i_phase, i_ion_pair, k_ion_pair)] +=
                  -gamma_i * ln_gamma_inter * CATION_N_(j_ion_pair) *
                  ANION_N_(j_ion_pair) / (omega * omega) * 2.0 *
                  (NUM_CATION_(k_ion_pair) + NUM_ANION_(k_ion_pair)) *
                  CATION_N_(k_ion_pair) / ANION_MW_(k_ion_pair) * 1.0e-6;
            }
          }
        }

      }  // Loop over ion pairs for partial derivatives

    }  // Loop on primary ion_pairs

  }  // Loop on aerosol phases
}
#endif

/** \brief Print the PDFiTE Activity sub model parameters
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 */
void sub_model_PDFiTE_print(int *sub_model_int_data,
                            double *sub_model_float_data) {
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  printf("\n\n*** PD-FiTE activity sub model ***\n");
#if 0
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
#endif
  printf("\n number of phases:     %d", NUM_PHASE_);
  printf("\n gas-phase water id:   %d", GAS_WATER_ID_);
  printf("\n number of ion pairs:  %d", NUM_ION_PAIRS_);
  printf("\n ** Phase data **");
  printf("\n   (Ids for ions and activity coefficients are relative");
  printf("\n    to the state id of aerosol-phase water for each phase.)");
  for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
    printf("\n   phase %d: aerosol-phase water state id: %d", i_phase,
           PHASE_ID_(i_phase));
  }
  printf("\n ** Ion pair data **");
  for (int i_pair = 0; i_pair < NUM_ION_PAIRS_; ++i_pair) {
    printf("\n Ion pair %d", i_pair);
    printf("\n   State ids (relative to aerosol-phase water)");
    printf("\n     Activity coeff: %d", ION_PAIR_ACT_ID_(i_pair));
    printf("\n     Cation:         %d", CATION_ID_(i_pair));
    printf("\n     Anion:          %d", ANION_ID_(i_pair));
    printf("\n   Jacobian aerosol water id (by phase):");
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
      printf(" %d", JAC_WATER_ID_(i_phase, i_pair));
    printf("\n   Cation stoichiometry:      %d", NUM_CATION_(i_pair));
    printf("\n   Anion stoichiometry:       %d", NUM_ANION_(i_pair));
    printf("\n   Cation MW (kg/mol):        %le", CATION_MW_(i_pair));
    ;
    printf("\n   Anion MW (kg/mol):         %le", ANION_MW_(i_pair));
    ;
    printf("\n   Cation concentration (mol m-3): %le", CATION_N_(i_pair));
    printf("\n   Anion concentration (mol m-3):  %le", ANION_N_(i_pair));
    printf("\n   Number of ion-pair interactions: %d", NUM_INTER_(i_pair));
    printf("\n   * Jacobian Omega ids *");
    for (int i_oip = 0; i_oip < NUM_ION_PAIRS_; ++i_oip) {
      printf("\n  Ion pair index: %d", i_oip);
      printf("\n    Anion  (by phase)");
      for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
        printf(" %d", JAC_ANION_ID_(i_phase, i_pair, i_oip));
      printf("\n    Cation (by phase)");
      for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
        printf(" %d", JAC_CATION_ID_(i_phase, i_pair, i_oip));
    }
    printf("\n   * Ion-pair interaction data *");
    for (int i_inter = 0; i_inter < NUM_INTER_(i_pair); ++i_inter) {
      printf("\n   Interaction %d:", i_inter);
      printf("\n     Ion pair index: %d", INTER_SPEC_ID_(i_pair, i_inter));
      printf("\n     Min RH: %le", MIN_RH_(i_pair, i_inter));
      printf("\n     Max RH: %le", MAX_RH_(i_pair, i_inter));
      printf("\n     Jacobian ids (by phase)");
      printf("\n       Cation:");
      for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
        printf(" %d", JAC_CATION_ID_(i_phase, i_pair,
                                     INTER_SPEC_ID_(i_pair, i_inter)));
      printf("\n       Anion: ");
      for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
        printf(" %d",
               JAC_ANION_ID_(i_phase, i_pair, INTER_SPEC_ID_(i_pair, i_inter)));
      printf("\n    ");
      for (int i_B = 0; i_B < NUM_B_(i_pair, i_inter); ++i_B)
        printf(" B%d: %le", i_B, B_Z_(i_pair, i_inter, i_B));
    }
  }
  printf("\n*** end PD-FiTE activity sub model ***");
}
