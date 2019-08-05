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
#define GAS_WATER_ID_ (int_data[1]-1)
#define NUM_ION_PAIRS_ (int_data[2])
#define INT_DATA_SIZE_ (int_data[3])
#define FLOAT_DATA_SIZE_ (int_data[4])
#define PPM_TO_RH_ (float_data[0])
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 1
#define PHASE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PAIR_INT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PAIR_FLOAT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIRS_+x]-1)
#define ION_PAIR_ACT_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)])
#define NUM_CATION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define NUM_ANION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+2])
#define CATION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+3])
#define ANION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+4])
#define NUM_INTER_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+5])
#define NUM_B_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+6+y])
#define INTER_SPEC_ID_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+6+NUM_INTER_(x)+y]-1)
#define INTER_SPEC_LOC_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+6+2*(NUM_INTER_(x))+y]-1)
#define CATION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define ANION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define CATION_N_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define ANION_N_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3])
#define MIN_RH_(x,y) (float_data[INTER_SPEC_LOC_(x,y)])
#define MAX_RH_(x,y) (float_data[INTER_SPEC_LOC_(x,y)+1])
#define B_Z_(x,y,z) (float_data[INTER_SPEC_LOC_(x,y)+2+z])


/** \brief Flag Jacobian elements used by this sub model
 *
 * activity sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_struct A matrix of flags for needed Jac elements
 */
void sub_model_PDFiTE_get_used_jac_elem(int *sub_model_int_data,
    double *sub_model_float_data, bool **jac_struct)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
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
    double *sub_model_float_data, int *deriv_ids, int **jac_ids)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}

/** \brief Get the id of a parameter in the condensed data block
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param identifiers For the PDFiTE model, the identifer is the index of the
 *                    ion pair for which the activity coefficient is needed
 *                    ( id = i_phase*num_ion_pair+i_ion_pair )
 * \param parameter_id Parameter id for the requested activity coefficient if
 *                     found
 */
void sub_model_PDFiTE_get_parameter_id(int *sub_model_int_data,
    double *sub_model_float_data, void *identifiers, int *parameter_id)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}

/** \brief Update sub model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param env_data Pointer to the environmental state array
 */
void sub_model_PDFiTE_update_env_state(int *sub_model_int_data,
    double *sub_model_float_data, double *env_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate PPM_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  double t_steam = 373.15; 				// steam temperature (K)
  double a = 1.0 - t_steam/TEMPERATURE_K_;

  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  double water_vp = 101325.0 * exp(a); 			// (Pa)

  PPM_TO_RH_ = PRESSURE_PA_ / water_vp / 1.0e6;		// (1/ppm)
}

/** \brief Perform the sub-model calculations for the current model state
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param model_data Pointer to the model data including the current state and
 *                   environmental conditions
 */
void sub_model_PDFiTE_calculate(int *sub_model_int_data,
    double *sub_model_float_data, ModelData *model_data)
{
  double *state = model_data->state;
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate the water activity---i.e., relative humidity (0-1)
  double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];

  // Keep a_w within 0-1
  // TODO Filter =( try to remove
  if (a_w<0.0) a_w = 0.0;
  if (a_w>1.0) a_w = 1.0;

  // Calculate ion_pair activity coefficients in each phase
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {

    // Initialize omega' (defined below)
    double omega_prime = 0.0;

    // Calculate the number of moles of each ion and omega' for the phase
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIRS_; i_ion_pair++) {

      // N (mol_i/m3) = c_i (ug/m3) / 10^6 (ug/g) / MW_i (ug/umol)
      CATION_N_(i_ion_pair) = state[PHASE_ID_(i_phase)+CATION_ID_(i_ion_pair)]
	      /CATION_MW_(i_ion_pair)/1000000.0;
      ANION_N_(i_ion_pair) = state[PHASE_ID_(i_phase)+ANION_ID_(i_ion_pair)]
	      /ANION_MW_(i_ion_pair)/1000000.0;

      // Calculate omega' (eq. 14 in \cite Topping 2009) as
      //   v_e^C = 2 * ( v_cation + v_anion) * N_cation * N_anion
      // where v_x is the stoichiometric coefficient for species x in
      // the ion_pair and N_x is its concentration summed across all
      // ion_pairs.
      // For each activity coefficient
      //   omega = omega' - omega_x
      // where omega_x is the contribution of ion_pair x to omega'
      omega_prime += 2.0 * ( NUM_CATION_(i_ion_pair) +
          NUM_ANION_(i_ion_pair) ) * CATION_N_(i_ion_pair) *
          ANION_N_(i_ion_pair);

    } // Loop on primary ion_pair

    // Calculate the activity coefficient
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIRS_; i_ion_pair++) {

      // If there are no interactions, the remaining ion pairs will not
      // have activity calculations (they only participate in interactions)
      if (NUM_INTER_(i_ion_pair)==0) break;

      // Calculate omega for this ion_pair
      // (eq. 15 in \cite{Topping2009})
      double omega = omega_prime - 2.0 * ( NUM_CATION_(i_ion_pair) +
      NUM_ANION_(i_ion_pair) ) * CATION_N_(i_ion_pair) *
      ANION_N_(i_ion_pair);

      // Initialize ln(gamma)
      double ln_gamma = 0.0;

      // Add contributions from each interacting ion_pair
      for (int i_inter=0; i_inter<NUM_INTER_(i_ion_pair); i_inter++) {

        // Only include interactions in the correct RH range
        // where the range is in (minRH, maxRH] except when a_w = 0.0
        // where the range is in [0.0, maxRH]
        if ((a_w<=MIN_RH_(i_ion_pair, i_inter) ||
             a_w>MAX_RH_(i_ion_pair, i_inter))
            &&
            !(a_w<=0.0 &&
             MIN_RH_(i_ion_pair, i_inter)<=0.0)
           ) continue;

        // Get the ion_pair id of the interacting species
        int j_ion_pair = INTER_SPEC_ID_(i_ion_pair, i_inter);

        // Calculate ln_gamma_inter
        double ln_gamma_inter = 0.0;
        for (int i_B=0; i_B<NUM_B_(i_ion_pair, i_inter); i_B++) {
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
          if (omega>0.0) {
            ln_gamma += ln_gamma_inter * CATION_N_(j_ion_pair) *
                     ANION_N_(j_ion_pair) / omega;
          }

        }

      } // Loop on interacting ion_pairs

      // Set the ion_pair activity
      state[PHASE_ID_(i_phase) + ION_PAIR_ACT_ID_(i_ion_pair)] =
              exp(ln_gamma);

    } // Loop on primary ion_pairs

  } // Loop on aerosol phases
}

// TODO finish adding J contributions
/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_row Pointer to the Jacobian row to modify
 */
void sub_model_PDFiTE_get_jac_contrib(int *sub_model_int_data,
    double *sub_model_float_data, double *jac_row)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}

/** \brief Print the PDFiTE Activity sub model parameters
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 */
void sub_model_PDFiTE_print(int *sub_model_int_data,
    double *sub_model_float_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  printf("\n\nPDFiTE Activity sub model\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef NUM_PHASE_
#undef GAS_WATER_ID_
#undef NUM_ION_PAIRS_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
#undef PPM_TO_RH_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef PHASE_ID_
#undef PAIR_INT_PARAM_LOC_
#undef PAIR_FLOAT_PARAM_LOC_
#undef ION_PAIR_ACT_ID_
#undef NUM_CATION_
#undef NUM_ANION_
#undef CATION_ID_
#undef ANION_ID_
#undef NUM_INTER_
#undef NUM_B_
#undef INTER_SPEC_ID_
#undef INTER_SPEC_LOC_
#undef CATION_MW_
#undef ANION_MW_
#undef CATION_N_
#undef ANION_N_
#undef MIN_RH_
#undef MAX_RH_
#undef B_Z_
