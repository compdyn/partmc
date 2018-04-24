/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * PDFiTE Activity reaction solver functions
 *
*/
/** \file
 * \brief PDFiTE Activity reaction solver functions
*/
#ifdef PMC_USE_SUNDIALS

#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define _NUM_PHASE_ (int_data[0])
#define _GAS_WATER_ID_ (int_data[1]-1)
#define _NUM_ELECTROLYTES_ (int_data[2])
#define _INT_DATA_SIZE_ (int_data[3])
#define _FLOAT_DATA_SIZE_ (int_data[4])
#define _ppm_TO_RH_ (float_data[0])
#define _NUM_INT_PROP_ 5
#define _NUM_REAL_PROP_ 1
#define _PHASE_ID_(x) (int_data[_NUM_INT_PROP_+x]-1)
#define _ELECT_INT_PARAM_LOC_(x) (int_data[_NUM_INT_PROP_+_NUM_PHASE_+x]-1)
#define _ELECT_FLOAT_PARAM_LOC_(x) (int_data[_NUM_INT_PROP_+_NUM_PHASE_+_NUM_ELECTROLYTES_+x]-1)
#define _ELECTROLYTE_ID_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)]-1)
#define _ELECTROLYTE_ACT_ID_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)+1]-1)
#define _NUM_CATION_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)+2])
#define _NUM_ANION_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)+3])
#define _CATION_ID_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)+4])
#define _ANION_ID_(x) (int_data[_ELECT_INT_PARAM_LOC_(x)+5])
#define _NUM_B_(x,y) (int_data[_ELECT_INT_PARAM_LOC_(x)+6+y])
#define _INTER_SPEC_ID_(x,y) (int_data[_ELECT_INT_PARAM_LOC_(x)+6+_NUM_ELECTROLYTES_+y]-1)
#define _B0_LOC_(x,y) (int_data[_ELECT_INT_PARAM_LOC_(x)+6+2*(_NUM_ELECTROLYTES_)+y])
#define _ELECTROLYTE_MW_(x) (float_data[_ELECT_FLOAT_PARAM_LOC_(x)])
#define _CATION_MW_(x) (float_data[_ELECT_FLOAT_PARAM_LOC_(x)+1])
#define _ANION_MW_(x) (float_data[_ELECT_FLOAT_PARAM_LOC_(x)+2])
#define _CATION_N_(x) (float_data[_ELECT_FLOAT_PARAM_LOC_(x)+3])
#define _ANION_N_(x) (float_data[_ELECT_FLOAT_PARAM_LOC_(x)+4])
#define _B_z_(x,y,z) (float_data[_B0_LOC_(x,y)+z])


/** \brief Flag Jacobian elements used by this reaction
 *
 * activity reactions are assumed to be at equilibrium
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_PDFiTE_activity_get_used_jac_elem(void *rxn_data, bool **jac_struct)
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
void * rxn_PDFiTE_activity_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data)
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
void * rxn_PDFiTE_activity_update_env_state(realtype *env_data, void *rxn_data)
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
void * rxn_PDFiTE_activity_pre_calc(ModelData *model_data, void *rxn_data)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the water activity---i.e., relative humidity (0-1)
  realtype a_w = _ppm_TO_RH_ * state[_GAS_WATER_ID_];

  // Calculate electrolyte activity coefficients in each phase
  for (int i_phase=0; i_phase<_NUM_PHASE_; i_phase++) {

    // Initialize omega' (defined below)
    realtype omega_prime = 0.0;
    
    // Calculate the number of moles of each ion and omega' for the phase
    for (int i_electrolyte=0; i_electrolyte<_NUM_ELECTROLYTES_; i_electrolyte++) {
      
      // N (mol_i/m3) = c_i (ug/m3) / 10^6 (ug/g) / MW_i (ug/umol) 
      _CATION_N_(i_electrolyte) = state[_PHASE_ID_(i_phase)+_CATION_ID_(i_electrolyte)]
	      /_CATION_MW_(i_electrolyte)/1000000.0;
      _ANION_N_(i_electrolyte) = state[_PHASE_ID_(i_phase)+_ANION_ID_(i_electrolyte)]
	      /_ANION_MW_(i_electrolyte)/1000000.0;

      // Calculate omega' (eq. 14 in \cite{Topping2009}) as
      //   v_e^C = 2 * ( v_cation + v_anion) * N_cation * N_anion
      // where v_x is the stoichiometric coefficient for species x in
      // the electrolyte and N_x is its concentration summed across all
      // electrolytes.
      // For each activity coefficient
      //   omega = omega' - omega_x
      // where omega_x is the contribution of electrolyte x to omega'
      omega_prime += 2.0 * ( _NUM_CATION_(i_electrolyte) + 
          _NUM_ANION_(i_electrolyte) ) * _CATION_N_(i_electrolyte) * 
          _ANION_N_(i_electrolyte);
    
    } // Loop on primary electrolyte

    // Calculate the activity coefficient
    for (int i_electrolyte=0; i_electrolyte<_NUM_ELECTROLYTES_; i_electrolyte++) {

      // Initialize ln(gamma)
      realtype ln_gamma = 0.0;

      // Add contributions from each interacting electrolyte
      for (int i_inter=0; i_inter<_NUM_ELECTROLYTES_; i_inter++) {

        // Get the electrolyte id of the interacting species
        int j_electrolyte = _INTER_SPEC_ID_(i_electrolyte, i_inter);

        // Calculate omega for this electrolyte
        // (eq. 15 in \cite{Topping2009})
        realtype omega = omega_prime - 2.0 * ( _NUM_CATION_(i_electrolyte) + 
        _NUM_ANION_(i_electrolyte) ) * _CATION_N_(i_electrolyte) * 
        _ANION_N_(i_electrolyte);

        // Calculate ln_gamma_inter
        realtype ln_gamma_inter = 0.0;
        for (int i_B=0; i_B<_NUM_B_(i_electrolyte, i_inter); i_B++) {
          ln_gamma_inter += _B_z_(i_electrolyte, i_inter, i_B) * pow(a_w, i_B);
        }

        // If this is the "self" interaction, ln_gamma_inter is ln(gamma_0A)
        // (eq. 15 in \cite{Topping2009})
        if (i_electrolyte == j_electrolyte) {
         
          // Add contribution to ln(gamma_A) from ln(gamma_0A)
          ln_gamma += ln_gamma_inter;
        }
        // ... otherwise it is d(ln(gamma_A))/g(N_B,M N_B,x)
        // (eq. 15 in \cite{Topping2009})
        else {
      
          // Add contribution to ln(gamma_A) from interacting electrolyte
          ln_gamma += ln_gamma_inter * _CATION_N_(j_electrolyte) *
                   _ANION_N_(j_electrolyte) / omega;

        }

      } // Loop on interacting electrolytes

      // Set the electrolyte activity
      state[_ELECTROLYTE_ACT_ID_(i_electrolyte)] = exp(ln_gamma);

    } // Loop on primary electrolytes

  } // Loop on aerosol phases

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
void * rxn_PDFiTE_activity_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
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
void * rxn_PDFiTE_activity_calc_jac_contrib(ModelData *model_data, realtype *J,
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
void * rxn_PDFiTE_activity_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the PDFiTE Activity reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_PDFiTE_activity_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nPDFiTE Activity reaction\n");
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
void * rxn_PDFiTE_activity_get_rate(void *rxn_data, realtype *state, realtype *env, realtype *rate)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *rate = 0.0;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _NUM_PHASE_
#undef _GAS_WATER_ID_
#undef _NUM_ELECTROLYTES_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_
#undef _ppm_TO_RH_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_ID_
#undef _ELECT_INT_PARAM_LOC_
#undef _ELECT_FLOAT_PARAM_LOC_
#undef _ELECTROLYTE_ID_
#undef _ELECTROLYTE_ACT_ID_
#undef _NUM_CATION_
#undef _NUM_ANION_
#undef _CATION_ID_
#undef _ANION_ID_
#undef _NUM_B_
#undef _B0_LOC_
#undef _ELECTROLYTE_MW_
#undef _CATION_MW_
#undef _ANION_MW_
#undef _CATION_N_
#undef _ANION_N_
#undef _B_z_

#endif
