/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * ZSR Aerosol Water reaction solver functions
 *
*/
/** \file
 * \brief ZSR Aerosol Water reaction solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

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
#define PHASE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PAIR_INT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PAIR_FLOAT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIR_+x]-1)
#define TYPE_(x) (int_data[PAIR_INT_PARAM_LOC_(x)])
#define JACOB_NUM_CATION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define JACOB_NUM_ANION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+2])
#define JACOB_CATION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+3])
#define JACOB_ANION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+4])
#define JACOB_NUM_Y_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+5])
#define EQSAM_NUM_ION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define EQSAM_ION_ID_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+2+y])
#define JACOB_low_RH_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define JACOB_CATION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define JACOB_ANION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define JACOB_Y_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])
#define EQSAM_NW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define EQSAM_ZW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define EQSAM_ION_PAIR_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define EQSAM_ION_MW_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])


/** \brief Flag Jacobian elements used by this reaction
 *
 * ZSR aerosol water reactions are assumed to be at equilibrium
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_get_used_jac_elem(void *rxn_data,
          bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * ZSR aerosol water reactions are assumed to be at equilibrium
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_update_env_state(double *env_data,
          void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Calculate PPM_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  double t_steam = 373.15; 		// steam temperature (K)
  double a = 1.0 - t_steam/TEMPERATURE_K_;

  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  double water_vp = 101325.0 * exp(a); 		// (Pa)

  PPM_TO_RH_ = PRESSURE_PA_ / water_vp / 1.0e6;	// (1/ppm)

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_pre_calc(ModelData *model_data, void *rxn_data)
{
  double *state = model_data->state;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being computed (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void * rxn_ZSR_aerosol_water_calc_deriv_contrib(ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void * rxn_ZSR_aerosol_water_calc_jac_contrib(ModelData *model_data,
          realtype *J, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the ZSR Aerosol Water reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_ZSR_aerosol_water_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nZSR Aerosol Water reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef SMALL_NUMBER_

#undef ACT_TYPE_JACOBSON
#undef ACT_TYPE_EQSAM

#undef NUM_PHASE_
#undef GAS_WATER_ID_
#undef NUM_ION_PAIR_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
#undef PPM_TO_RH_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef PHASE_ID_
#undef PAIR_INT_PARAM_LOC_
#undef PAIR_FLOAT_PARAM_LOC_
#undef TYPE_
#undef JACOB_NUM_CATION_
#undef JACOB_NUM_ANION_
#undef JACOB_CATION_ID_
#undef JACOB_ANION_ID_
#undef JACOB_NUM_Y_
#undef EQSAM_NUM_ION_
#undef EQSAM_ION_ID_
#undef JACOB_low_RH_
#undef JACOB_CATION_MW_
#undef JACOB_ANION_MW_
#undef JACOB_Y_
#undef EQSAM_NW_
#undef EQSAM_ZW_
#undef EQSAM_ION_PAIR_MW_
#undef EQSAM_ION_MW_
