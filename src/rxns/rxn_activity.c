/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Activity reaction solver functions
 *
*/
/** \file
 * \brief Activity reaction solver functions
*/
#ifdef PMC_USE_SUNDIALS

#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define ACT_TYPE_JACOBSON 1
#define ACT_TYPE_EQSAM 2

#define _NUM_PHASE_ (int_data[0])
#define _GAS_WATER_ID_ (int_data[1]-1)
#define _NUM_ION_PAIR_ (int_data[2])
#define _INT_DATA_SIZE_ (int_data[3])
#define _FLOAT_DATA_SIZE_ (int_data[4])
#define _ppm_TO_RH_ (float_data[0])
#define _NUM_INT_PROP_ 5
#define _NUM_REAL_PROP_ 1
#define _PHASE_ID_(x) (int_data[_NUM_INT_PROP_+x]-1)
#define _PAIR_INT_PARAM_LOC_(x) (int_data[_NUM_INT_PROP_+_NUM_PHASE_+x]-1)
#define _PAIR_FLOAT_PARAM_LOC_(x) (int_data[_NUM_INT_PROP_+_NUM_PHASE_+_NUM_ION_PAIR_+x]-1)
#define _TYPE_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)])
#define _JACOB_NUM_CATION_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+1])
#define _JACOB_NUM_ANION_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+2])
#define _JACOB_CATION_ID_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+3])
#define _JACOB_ANION_ID_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+4])
#define _JACOB_NUM_Y_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+5])
#define _EQSAM_NUM_ION_(x) (int_data[_PAIR_INT_PARAM_LOC_(x)+1])
#define _EQSAM_ION_ID_(x,y) (int_data[_PAIR_INT_PARAM_LOC_(x)+2+y])
#define _JACOB_low_RH_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)])
#define _JACOB_CATION_MW_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+1])
#define _JACOB_ANION_MW_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+2])
#define _JACOB_Y_(x,y) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+3+y])
#define _EQSAM_NW_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)])
#define _EQSAM_ZW_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+1])
#define _EQSAM_ION_PAIR_MW_(x) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+2])
#define _EQSAM_ION_MW_(x,y) (float_data[_PAIR_FLOAT_PARAM_LOC_(x)+3+y])


/** \brief Flag Jacobian elements used by this reaction
 *
 * activity reactions are assumed to be at equilibrium
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}  

/** \brief Update the time derivative and Jacbobian array indices
 *
 * activity reactions are assumed to be at equilibrium
 *
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate _ppm_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  realtype t_steam = 373.15; 				// steam temperature (K)
  realtype a = 1.0 - t_steam/_TEMPERATURE_K_;

  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  realtype water_vp = 101325.0 * exp(a); 			// (Pa)
  
  _ppm_TO_RH_ = _PRESSURE_PA_ / water_vp / 1.0e6;		// (1/ppm)

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_pre_calc(ModelData *model_data, void *rxn_data)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the water activity---i.e., relative humidity (0-1)
  realtype a_w = _ppm_TO_RH_ * state[_GAS_WATER_ID_];

  // Calculate the total aerosol water for each instance of the aerosol phase
  for (int i_phase=0; i_phase<_NUM_PHASE_; i_phase++) {
    realtype *water = &(state[_PHASE_ID_(i_phase)]);
    *water = 0.0;

    // Get the contribution from each ion pair
    for (int i_ion_pair=0; i_ion_pair<_NUM_ION_PAIR_; i_ion_pair++) {
      
      realtype molality;
      realtype j_aw;

      // Determine which type of activity calculation should be used
      switch (_TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Determine whether to use the minimum RH in the calculation
          j_aw = a_w>_JACOB_low_RH_(i_ion_pair) ? a_w : _JACOB_low_RH_(i_ion_pair);

          // Calculate the molality of the pure binary ion pair solution
	  molality = 0.0;
          for (int i_order=0; i_order<_JACOB_NUM_Y_(i_ion_pair); i_order++) 
		  molality += _JACOB_Y_(i_ion_pair, i_order) * pow(j_aw,i_order);
          molality *= molality; // (mol/kg)

	  // Calculate the water associated with this ion pair
          realtype cation = state[_PHASE_ID_(i_phase)+_JACOB_CATION_ID_(i_ion_pair)]
		  /_JACOB_NUM_CATION_(i_ion_pair)/_JACOB_CATION_MW_(i_ion_pair);
          realtype anion = state[_PHASE_ID_(i_phase)+_JACOB_ANION_ID_(i_ion_pair)]
		  /_JACOB_NUM_ANION_(i_ion_pair)/_JACOB_ANION_MW_(i_ion_pair); // (umol/m3)
	  *water += (cation>anion ? anion : cation) / molality * 1000.0; // (ug/m3)
          
	  break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

	  // Calculate the molality of the ion pair
	  molality = (_EQSAM_NW_(i_ion_pair) * 55.51 * 18.01 / _EQSAM_ION_PAIR_MW_(i_ion_pair) *
		  (1.0/a_w-1.0));
	  molality = pow(molality, _EQSAM_ZW_(i_ion_pair)); // (mol/kg)

	  // Calculate the water associated with this ion pair
	  for (int i_ion=0; i_ion<_EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
            *water += state[_PHASE_ID_(i_phase)+_EQSAM_ION_ID_(i_ion_pair,i_ion)]
		    /_EQSAM_ION_MW_(i_ion_pair,i_ion) / molality * 1000.0; // (ug/m3);
	  }

	  break;
      }
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
		void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param state Pointer to the state array
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_calc_jac_contrib(ModelData *model_data, realtype *J,
		void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the Activity reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nActivity reaction\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Return the reaction rate for the current conditions
 *
 * activity reactions are assumed to be at equilibrium
 *
 * \param rxn_data Pointer to the reaction data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \param rate Pointer to a double value to store the calculated rate
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_activity_get_rate(void *rxn_data, realtype *state, realtype *env, realtype *rate)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *rate = 0.0;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef ACT_TYPE_JACOBSON
#undef ACT_TYPE_EQSAM

#undef _NUM_PHASE_
#undef _GAS_WATER_ID_
#undef _NUM_ION_PAIR_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_
#undef _ppm_TO_RH_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_ID_
#undef _PAIR_INT_PARAM_LOC_
#undef _PAIR_FLOAT_PARAM_LOC_
#undef _TYPE_
#undef _JACOB_NUM_CATION_
#undef _JACOB_NUM_ANION_
#undef _JACOB_CATION_ID_
#undef _JACOB_ANION_ID_
#undef _JACOB_NUM_Y_
#undef _EQSAM_NUM_ION_
#undef _EQSAM_ION_ID_
#undef _JACOB_low_RH_
#undef _JACOB_CATION_MW_
#undef _JACOB_ANION_MW_
#undef _JACOB_Y_
#undef _EQSAM_NW_
#undef _EQSAM_ZW_
#undef _EQSAM_ION_PAIR_MW_
#undef _EQSAM_ION_MW_

#endif
