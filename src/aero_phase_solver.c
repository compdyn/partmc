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
#include "aero_phase_solver.h"

// TODO move all shared constants to a common header file
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

#define NUM_STATE_VAR_ (int_data[0])
#define NUM_INT_PROP_ 1
#define NUM_FLOAT_PROP_ 0
#define SPEC_TYPE_(x) (int_data[NUM_INT_PROP_+x])
#define MW_(x) (float_data[NUM_FLOAT_PROP_+x])
#define DENSITY_(x) (float_data[NUM_FLOAT_PROP_+NUM_STATE_VAR_+x])
#define INT_DATA_SIZE_ (NUM_INT_PROP_+NUM_STATE_VAR_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+2*NUM_STATE_VAR_)

/** \brief Get the mass and average MW in an aerosol phase
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to use in the calculation
 * \param state_var Pointer the aerosol phase on the state variable array
 * \param mass Pointer to hold total aerosol phase mass
 *             (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$ or
 *              \f$\mbox{\si{\micro\gram\per particle}}\f$)
 * \param MW Pointer to hold average MW of the aerosol phase
 *           (\f$\mbox{\si{\kilogram\per\mol}}\f$)
 * \return A pointer to a set of partial derivatives \f$\frac{dm}{dy}\f$, or a
 *         NULL pointer if no partial derivatives exist
 */
void * aero_phase_get_mass(ModelData *model_data, int aero_phase_idx,
        double *state_var, double *mass, double *MW)
{

  // Set up a pointer for the partial derivatives
  void *partial_deriv = NULL;

  // Get the requested aerosol phase data
  int *int_data = (int*) aero_phase_find(model_data, aero_phase_idx);
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Sum the mass and MW
  *mass = 0.0;
  *MW = 0.0;
  for (int i_spec=0; i_spec<NUM_STATE_VAR_; i_spec++) {
    if (SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      *mass += state_var[i_spec];
      *MW += state_var[i_spec] / MW_(i_spec);
    }
  }
  *MW = *mass / *MW;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the volume of an aerosol phase
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to use in the calculation
 * \param state_var Pointer to the aerosol phase on the state variable array
 * \param volume Pointer to hold the aerosol phase volume
 *               (\f$\mbox{\si{\cubic\metre\per\cubic\metre}}\f$ or
 *                \f$\mbox{\si{\cubic\metre\per particle}}\f$)
 * \return A pointer to a set of partial derivatives \f$\frac{dv}{dy}\f$, or
 *         a NULL pointer if no partial derivatives exist
 */
void * aero_phase_get_volume(ModelData *model_data, int aero_phase_idx,
          double *state_var, double *volume)
{

  // Set up a pointer for the partial derivatives
  void *partial_deriv = NULL;

  // Get the requested aerosol phase data
  int *int_data = (int*) aero_phase_find(model_data, aero_phase_idx);
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Sum the mass and MW
  *volume = 0.0;
  for (int i_spec=0; i_spec<NUM_STATE_VAR_; i_spec++) {
    if (SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      *volume += state_var[i_spec] / 1.0e9 / DENSITY_(i_spec);;
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \return The aero_phase_data pointer advanced by the size of the aerosol
 *         phase
 */
void * aero_phase_skip(void *aero_phase_data)
{
  int *int_data = (int*) aero_phase_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

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
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *aero_phase_data = (int*) (model_data->nxt_aero_phase);

  // Add the integer parameters
  for (; n_int_param>0; n_int_param--) *(aero_phase_data++) = *(int_param++);

  // Add the floating-point parameters
  double *flt_ptr = (double*) aero_phase_data;
  for (; n_float_param>0; n_float_param--)
          *(flt_ptr++) = (double) *(float_param++);

  // Set the pointer for the next free space in aero_phase_data;
  model_data->nxt_aero_phase = (void*) flt_ptr;
}

/** \brief Print the aerosol phase data
 * \param solver_data Pointer to the solver data
 */
void aero_phase_print_data(void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *aero_phase_data = (int*) (model_data->aero_phase_data);

  // Get the number of aerosol phases
  int n_aero_phase = *(aero_phase_data++);

  // Loop through the aerosol phases and print their data
  // advancing the aero_phase_data pointer each time
  for (int i_aero_phase=0; i_aero_phase<n_aero_phase; i_aero_phase++) {
    int *int_data = (int*) aero_phase_data;
    double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

    printf("\n\nAerosol Phase %d\n\n", i_aero_phase);
    printf("\nint_data");
    for (int i=0; i<INT_DATA_SIZE_; i++) printf(" %d", int_data[i]);
    printf("\nfloat_data");
    for (int i=0; i<FLOAT_DATA_SIZE_; i++) printf(" %le", float_data[i]);

    aero_phase_data = (int*) &(float_data[FLOAT_DATA_SIZE_]);
  }
}

#undef NUM_STATE_VAR_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef SPEC_TYPE_
#undef MW_
#undef DENSITY_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_

