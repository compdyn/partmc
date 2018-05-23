/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aerosol phase-specific functions for use by the solver
 *
 */
/** \file
 * \brief Aerosol phase functions
 */
#include "phlex_solver.h"

// TODO move all shared constants to a common header file
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

#define _NUM_STATE_VAR_ (int_data[0])
#define _NUM_INT_PROP_ 1
#define _NUM_FLOAT_PROP_ 0
#define _SPEC_TYPE_(x) (int_data[_NUM_INT_PROP_+x])
#define _MW_(x) (float_data[_NUM_FLOAT_PROP_+x])
#define _DENSITY_(x) (float_data[_NUM_FLOAT_PROP_+_NUM_STATE_VAR_+x])
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+_NUM_STATE_VAR_)
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_+2*_NUM_STATE_VAR_)

#ifdef PMC_USE_SUNDIALS

/** \brief Get the mass and average MW in an aerosol phase
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to use in the calculation
 * \param state_var Pointer the aerosol phase on the state variable array
 * \param mass Pointer to hold total aerosol phase mass (ug/m^3 or ug/particle)
 * \param MW Pointer to hold average MW of the aerosol phase (kg/mol)
 * \return A pointer to a set of partial derivatives dm/dy, or a NULL pointer if no
 *         partial derivatives exist
 */
void * aero_phase_get_mass(ModelData *model_data, int aero_phase_idx, realtype *state_var,
                realtype *mass, realtype *MW)
{

  // Set up a pointer for the partial derivatives
  void *partial_deriv = NULL;

  // Get the requested aerosol phase data
  int *int_data = (int*) aero_phase_find(model_data, aero_phase_idx);
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Sum the mass and MW
  *mass = 0.0;
  *MW = 0.0;
  for (int i_spec=0; i_spec<_NUM_STATE_VAR_; i_spec++) {
    if (_SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        _SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        _SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      *mass += state_var[i_spec];
      *MW += state_var[i_spec] / _MW_(i_spec);
    }
  }
  *MW = *mass / *MW;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Get the volume of an aerosol phase
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to use in the calculation
 * \param state_var Pointer to the aerosol phase on the state variable array
 * \param volume Pointer to hold the aerosol phase volume (m^3/m^3 or m^3/particle)
 * \return A pointer to a set of partial derivatives dv/dy, or a NULL pointer if no
 *         partial derivatives exist
 */
void * aero_phase_get_volume(ModelData *model_data, int aero_phase_idx, realtype *state_var,
                realtype *volume)
{

  // Set up a pointer for the partial derivatives
  void *partial_deriv = NULL;

  // Get the requested aerosol phase data
  int *int_data = (int*) aero_phase_find(model_data, aero_phase_idx);
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Sum the mass and MW
  *volume = 0.0;
  for (int i_spec=0; i_spec<_NUM_STATE_VAR_; i_spec++) {
    if (_SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        _SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        _SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      *volume += state_var[i_spec] / 1.0e9 / _DENSITY_(i_spec);;
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Find an aerosol phase in the list
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the desired aerosol phase
 * \return A pointer to the requested aerosol phase
 */
void * aero_phase_find(ModelData *model_data, int aero_phase_idx)
{

  // Get the number of aerosol phases
  int *aero_phase_data = (int*) (model_data->aero_phase_data);
  int n_aero_phase = *(aero_phase_data++);

  // Loop through the aerosol phases to find the one requested
  for (int i_aero_phase=0; i_aero_phase<aero_phase_idx; i_aero_phase++) {

    // Advance the pointer to the next aerosol phase
    aero_phase_data = aero_phase_skip((void*) aero_phase_data);

  }

  return (void*) aero_phase_data;

}

/** \brief Skip over an aerosol phase
 *
 * \param aero_phase_data Pointer to the aerosol phase to skip over
 * \return The aero_phase_data pointer advanced by the size of the aerosol phase
 */
void * aero_phase_skip(void *aero_phase_data)
{
  int *int_data = (int*) aero_phase_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#endif

/** \brief Add condensed data to the condensed data block for aerosol phases
 *
 * \param n_int_param Number of integer parameters
 * \param n_float_param Number of floating-point parameters
 * \param int_param Pointer to the integer parameter array
 * \param float_param Pointer to the floating-point parameter array
 * \param solver_data Pointer to the solver data
 */
void aero_phase_add_condensed_data(int n_int_param, int n_float_param,
              int *int_param, double *float_param, void *solver_data)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *aero_phase_data = (int*) (model_data->nxt_aero_phase);

#ifdef PMC_USE_SUNDIALS

  // Add the integer parameters
  for (; n_int_param>0; n_int_param--) *(aero_phase_data++) = *(int_param++);

  // Add the floating-point parameters
  realtype *flt_ptr = (realtype*) aero_phase_data;
  for (; n_float_param>0; n_float_param--) *(flt_ptr++) = (realtype) *(float_param++);

  // Set the pointer for the next free space in aero_phase_data;
  model_data->nxt_aero_phase = (void*) flt_ptr;

#endif
}

#undef _NUM_STATE_VAR_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _SPEC_TYPE_
#undef _MW_
#undef _DENSITY_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

