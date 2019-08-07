/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * ZSR Aerosol Water sub model solver functions
 *
*/
/** \file
 * \brief ZSR Aerosol Water sub model solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../sub_models.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_state[0]
#define PRESSURE_PA_ env_state[1]

#define SMALL_NUMBER_ 1.0e-30

#define ACT_TYPE_JACOBSON 1
#define ACT_TYPE_EQSAM 2

#define NUM_PHASE_ (int_data[0])
#define GAS_WATER_ID_ (int_data[1]-1)
#define NUM_ION_PAIR_ (int_data[2])
#define INT_DATA_SIZE_ (int_data[3])
#define FLOAT_DATA_SIZE_ (int_data[4])
#define PPM_TO_RH_ (float_data[0])
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 1
#define PHASE_ID_(p) (int_data[NUM_INT_PROP_+p]-1)
#define PAIR_INT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PAIR_FLOAT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIR_+x]-1)
#define TYPE_(x) (int_data[PAIR_INT_PARAM_LOC_(x)])
#define JACOB_NUM_CATION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define JACOB_NUM_ANION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+2])
#define JACOB_CATION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+3])
#define JACOB_ANION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+4])
#define JACOB_NUM_Y_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+5])
#define JACOB_GAS_WATER_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+p]
#define JACOB_CATION_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+NUM_PHASE_+p]
#define JACOB_ANION_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+2*NUM_PHASE_+p]
#define EQSAM_NUM_ION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define EQSAM_GAS_WATER_JAC_ID_(p,x) (int_data[PAIR_INT_PARAM_LOC_(x)+2+p])
#define EQSAM_ION_ID_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+2+NUM_PHASE_+y])
#define EQSAM_ION_JAC_ID_(p,x,y) int_data[PAIR_INT_PARAM_LOC_(x)+2+NUM_PHASE_+EQSAM_NUM_ION_(x)+y*NUM_PHASE_+p]
#define JACOB_low_RH_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define JACOB_CATION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define JACOB_ANION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define JACOB_Y_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])
#define EQSAM_NW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define EQSAM_ZW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define EQSAM_ION_PAIR_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define EQSAM_ION_MW_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])

// Update types (These must match values in sub_model_UNIFAC.F90)
// (none right now)

/** \brief Flag Jacobian elements used by this sub model
 *
 * ZSR aerosol water sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_struct A matrix of flags for needed Jac elements
 */
void sub_model_ZSR_aerosol_water_get_used_jac_elem(int *sub_model_int_data,
    double *sub_model_float_data, bool **jac_struct)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Loop through the dependent species - aerosol water and set all Jacobian
  // elements for each
  for (int i_phase=0; i_phase < NUM_PHASE_; ++i_phase) {

    // Flag the gas-phase water species
    jac_struct[PHASE_ID_(i_phase)][GAS_WATER_ID_] = true;

    // Flag elements for each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; ++i_ion_pair) {

      // Flag aerosol elements by calculation type
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Flag the anion and cation Jacobian elements
          jac_struct[PHASE_ID_(i_phase)]
                    [PHASE_ID_(i_phase)+JACOB_CATION_ID_(i_ion_pair)] = true;
          jac_struct[PHASE_ID_(i_phase)]
                    [PHASE_ID_(i_phase)+JACOB_ANION_ID_(i_ion_pair)] = true;
          break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Flag the ion Jacobian elements
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
            jac_struct[PHASE_ID_(i_phase)]
                      [PHASE_ID_(i_phase)+
                          EQSAM_ION_ID_(i_ion_pair,i_ion)] = true;
          }
          break;
      }
    }
  }
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * ZSR aerosol water sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param deriv_ids Indices for state array variables on the solver state array
 * \param jac_ids Indices for Jacobian elements in the sparse data array
 */
void sub_model_ZSR_aerosol_water_update_ids(int *sub_model_int_data,
    double *sub_model_float_data, int *deriv_ids, int **jac_ids)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Loop through the dependent species - aerosol water and set all Jacobian
  // elements for each
  for (int i_phase=0; i_phase < NUM_PHASE_; ++i_phase) {

    // Flag elements for each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; ++i_ion_pair) {

      // Flag aerosol elements by calculation type
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Save the gas-phase water species
          JACOB_GAS_WATER_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)][GAS_WATER_ID_];

          // Save the cation and anion Jacobian elements
          JACOB_CATION_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)]
                     [PHASE_ID_(i_phase)+JACOB_CATION_ID_(i_ion_pair)];
          JACOB_ANION_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)]
                     [PHASE_ID_(i_phase)+JACOB_ANION_ID_(i_ion_pair)];
          break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Save the gas-phase water species
          EQSAM_GAS_WATER_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)][GAS_WATER_ID_];

          // Save the ion Jacobian elements
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
            EQSAM_ION_JAC_ID_(i_phase,i_ion_pair,i_ion) =
                jac_ids[PHASE_ID_(i_phase)]
                       [PHASE_ID_(i_phase)+EQSAM_ION_ID_(i_ion_pair,i_ion)];
          }
          break;
      }
    }
  }
}

/** \brief Get the id of a parameter in the condensed data block
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param identifiers For the ZSR model, the identifer is just the id
 *                    on the state array of the phase for which water is being
 *                    calculated (not necessarily the state id of the water
 *                    species).
 * \param parameter_id Parameter id for the requested aerosol-phase water if
 *                     found
 */
void sub_model_ZSR_aerosol_water_get_parameter_id(int *sub_model_int_data,
    double *sub_model_float_data, void *identifiers, int *parameter_id)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}
/** \brief Update sub model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param env_state Pointer to the environmental state array
 */
void sub_model_ZSR_aerosol_water_update_env_state(int *sub_model_int_data,
    double *sub_model_float_data, double *env_state)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate PPM_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  double t_steam = 373.15; 		// steam temperature (K)
  double a = 1.0 - t_steam/TEMPERATURE_K_;

  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  double water_vp = 101325.0 * exp(a); 		// (Pa)

  PPM_TO_RH_ = PRESSURE_PA_ / water_vp / 1.0e6;	// (1/ppm)
}

/** \brief Do pre-derivative calculations
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param model_data Pointer to the model data, including the state array
 */
void sub_model_ZSR_aerosol_water_calculate(int *sub_model_int_data,
    double *sub_model_float_data, ModelData *model_data)
{
  double *state = model_data->state;
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate the water activity---i.e., relative humidity (0-1)
  double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];

  // Calculate the total aerosol water for each instance of the aerosol phase
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {
    double *water = &(state[PHASE_ID_(i_phase)]);
    *water = 0.0;

    // Get the contribution from each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; i_ion_pair++) {

      double molality;
      double j_aw, e_aw;
      double conc;

      // Determine which type of activity calculation should be used
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Determine whether to use the minimum RH in the calculation
          j_aw = a_w>JACOB_low_RH_(i_ion_pair) ? a_w :
                  JACOB_low_RH_(i_ion_pair);

          // Calculate the molality of the pure binary ion pair solution
	  molality = 0.0;
          for (int i_order=0; i_order<JACOB_NUM_Y_(i_ion_pair); i_order++)
		  molality += JACOB_Y_(i_ion_pair, i_order) *
                          pow(j_aw,i_order);
          molality *= molality; // (mol/kg)

	  // Calculate the water associated with this ion pair
          double cation = state[PHASE_ID_(i_phase) +
                  JACOB_CATION_ID_(i_ion_pair)] /
		  JACOB_NUM_CATION_(i_ion_pair) /
                  JACOB_CATION_MW_(i_ion_pair);
          double anion = state[PHASE_ID_(i_phase) +
                  JACOB_ANION_ID_(i_ion_pair)] /
		  JACOB_NUM_ANION_(i_ion_pair) /
                  JACOB_ANION_MW_(i_ion_pair); // (umol/m3)
	  conc = (cation>anion ? anion : cation);
          conc = (conc>0.0 ? conc : 0.0);
          *water += conc / molality * 1000.0; // (ug/m3)

	  break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Keep the water activity within the range specified in EQSAM
          e_aw = a_w > 0.99 ? 0.99 : a_w;
          e_aw = e_aw < 0.001 ? 0.001 : e_aw;

	  // Calculate the molality of the ion pair
	  molality = (EQSAM_NW_(i_ion_pair) * 55.51 * 18.01 /
                    EQSAM_ION_PAIR_MW_(i_ion_pair) * (1.0/e_aw-1.0));
	  molality = pow(molality, EQSAM_ZW_(i_ion_pair)); // (mol/kg)

	  // Calculate the water associated with this ion pair
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
	    conc = state[PHASE_ID_(i_phase)+EQSAM_ION_ID_(i_ion_pair,i_ion)];
            conc = (conc>0.0 ? conc : 0.0);
            *water += conc / EQSAM_ION_MW_(i_ion_pair,i_ion) /
                    molality * 1000.0; // (ug/m3);
	  }

	  break;
      }
    }
    *water = (*water > SMALL_NUMBER_) ? *water : 0.0;
  }
}

// TODO finish adding J contributions
/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_row Pointer to the Jacobian row to modify
 */
void sub_model_ZSR_aerosol_water_get_jac_contrib(int *sub_model_int_data,
    double *sub_model_float_data, double *jac_row)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}

/** \brief Print the ZSR Aerosol Water sub model parameters
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 */
void sub_model_ZSR_aerosol_water_print(int *sub_model_int_data,
    double *sub_model_float_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  printf("\n\nZSR aerosol water sub model\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
}
