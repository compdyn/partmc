/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Single particle aerosol representation functions
 *
 */
/** \file
 * \brief Single particle aerosol representation functions
 */
#ifdef PMC_USE_SUNDIALS

#include "../aero_rep_solver.h"

// TODO Lookup environmental indicies during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define _NUM_PHASE_ int_data[0]
#define _RADIUS_ float_data[0]
#define _NUMBER_CONC_ float_data[1]
#define _NUM_INT_PARAM_ 1
#define _NUM_FLOAT_PARAM_ 2
#define _PHASE_STATE_ID_(x) (int_data[_NUM_INT_PARAM_+x]-1)
#define _PHASE_MODEL_DATA_ID_(x) (int_data[_NUM_INT_PARAM_+_NUM_PHASE_+x]-1)
#define _AERO_PHASE_MASS_(x) (float_data[_NUM_FLOAT_PARAM_+x])
#define _AERO_PHASE_AVG_MW_(x) (float_data[_NUM_FLOAT_PARAM_+_NUM_PHASE_+x])
#define _INT_DATA_SIZE_ (_NUM_INT_PARAM_+2*_NUM_PHASE_)
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PARAM_+2*_NUM_PHASE_)

// Update types (These must match values in aero_rep_single_particle.F90)
#define UPDATE_RADIUS 0
#define UPDATE_NUMBER_CONC 1

/** \brief Flag elements on the state array used by this aerosol representation
 *
 * The single particle aerosol representation functions do not use state array values
 *
 * \param aero_rep_data A pointer to the aerosol representation data
 * \param state_flags Array of flags indicating state array elements used
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_single_particle_get_dependencies(void *aero_rep_data, bool *state_flags)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update aerosol representation data for new environmental conditions
 *
 * The single particle aerosol representation does not use environmental conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_update_env_state(double *env_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update aerosol representation data for a new state
 *
 * Updates the mass and average MW for each aerosol phase in the particle
 *
 * \param model_data Pointer to the model data, include the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_update_state(ModelData *model_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the total aerosol phase masses
  for (int i_phase=0; i_phase<_NUM_PHASE_; i_phase++) {

    // Get a pointer to the phase on the state array
    realtype *state_var = (realtype*) (model_data->state);
    state_var += _PHASE_STATE_ID_(i_phase);

    // Get the mass and average MW
    aero_phase_get_mass(model_data, _PHASE_MODEL_DATA_ID_(i_phase), state_var, 
               &(_AERO_PHASE_MASS_(i_phase)), &(_AERO_PHASE_AVG_MW_(i_phase)));
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the effective particle radius
 *
 * The single particle radius is set by the aerosol model prior to solving the chemistry. 
 * Thus, all dr/dy are zero. Also, there is only one set of particles in the single particle
 * representation, so the phase index is not used.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param radius Effective particle radius (m)
 * \param partial_deriv dr/dy where y are species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_get_effective_radius(int aero_phase_idx, double *radius, 
		double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *radius = _RADIUS_;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the particle number concentration
 *
 * This single particle number concentration is set by the aerosol model prior to solving the chemistry.
 * Thus, all dn/dy are zero. Also, there is only one set of particles in the single particle representation,
 * so the phase index is not used.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param number_conc Particle number concentration (#/cm^3)
 * \param partial_deriv dn/dy where y are the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_get_number_conc(int aero_phase_idx, double *number_conc, 
		double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *number_conc = _NUMBER_CONC_;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the type of aerosol concentration type used.
 *
 * Single particle concentrations are per-particle.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_conc_type Pointer to int that will hold the concentration type code
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_get_aero_conc_type(int aero_phase_idx, int *aero_conc_type, 
		void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *aero_conc_type = 0;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}
  
/** Get the total mass in an aerosol phase
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_phase_mass Total mass in the aerosol phase (ug/m^3)
 * \param aero_phase_avg_MW Average molecular weight in the aerosol phase (kg/mol)
 * \param partial_deriv dn/dy where y are the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_get_aero_phase_mass(int aero_phase_idx, double *aero_phase_mass,
		double *aero_phase_avg_MW, double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *aero_phase_mass = _AERO_PHASE_MASS_(aero_phase_idx);
  *aero_phase_avg_MW = _AERO_PHASE_AVG_MW_(aero_phase_idx);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update the aerosol representation data
 *
 *  The single particle aerosol representation has two update types:
 *
 *  UPDATE_RADIUS : where the update data should point to a single floating-point
 *  variable holding the new particle radius
 *
 *  UPDATE_NUMBER_CONC : where the update data should point to a single floating-point
 *  variable holding the new particle number concentration
 *
 * \param update_type The type of update to perform
 * \param update_data Pointer to the data required for the update
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_single_particle_update_data(int update_type, void *update_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  switch (update_type) {
    case UPDATE_RADIUS :
      _RADIUS_ = *((realtype*)update_data);
      break;
    case UPDATE_NUMBER_CONC:
      _NUMBER_CONC_ = *((realtype*)update_data);
      break;
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the Single Particle reaction parameters
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_single_particle_print(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nSingle particle aerosol representation\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Advance the aerosol representation data pointer to the next aerosol representation
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_single_particle_skip(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _NUM_PHASE_
#undef _RADIUS_
#undef _NUMBER_CONC_
#undef _NUM_INT_PARAM_
#undef _NUM_FLOAT_PARAM_
#undef _PHASE_STATE_ID_
#undef _PHASE_MODEL_DATA_ID_
#undef _AERO_PHASE_MASS_
#undef _AERO_PHASE_AVG_MW_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
