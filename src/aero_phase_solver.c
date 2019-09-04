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
#include <stdio.h>
#include <stdlib.h>
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

/** \brief Flag Jacobian elements used in calculations of mass and volume
 *
 * \param model_data Pointer to the model data(state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to find elements for
 * \param state_var_id Index in the state array for this aerosol phase
 * \param jac_struct 1D array of flags indicating potentially non-zero
 *                   Jacobian elements. (The dependent variable should have
 *                   been chosen by the calling function.)
 * \return Number of Jacobian elements flagged
 */
int aero_phase_get_used_jac_elem(ModelData *model_data, int aero_phase_idx,
          int state_var_id, bool *jac_struct)
{

  // Get the requested aerosol phase data
  int *int_data      = model_data->aero_phase_int_ptrs[aero_phase_idx];
  double *float_data = model_data->aero_phase_float_ptrs[aero_phase_idx];

  int num_flagged_elem = 0;

  for (int i_spec=0; i_spec<NUM_STATE_VAR_; i_spec++) {
    if (SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      jac_struct[state_var_id+i_spec] = true;
      num_flagged_elem++;
    }
  }

  return num_flagged_elem;
}

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
 * \param jac_elem_mass When not NULL, a pointer to an array whose length is the
 *                 number of Jacobian elements used in calculations of mass and
 *                 volume of this aerosol phase returned by
 *                 \c aero_phase_get_used_jac_elem and whose contents will be
 *                 set to the partial derivatives of mass by concentration
 *                 \f$\frac{dm}{dy_i}\f$ of each component species \f$y_i\f$.
 * \param jac_elem_MW When not NULL, a pointer to an array whose length is the
 *                 number of Jacobian elements used in calculations of mass and
 *                 volume of this aerosol phase returned by
 *                 \c aero_phase_get_used_jac_elem and whose contents will be
 *                 set to the partial derivatives of total molecular weight by
 *                 concentration \f$\frac{dMW}{dy_i}\f$ of each component
 *                 species \f$y_i\f$.
 */
void aero_phase_get_mass(ModelData *model_data, int aero_phase_idx,
        double *state_var, double *mass, double *MW, double *jac_elem_mass,
        double *jac_elem_MW)
{

  // Get the requested aerosol phase data
  int *int_data      = model_data->aero_phase_int_ptrs[aero_phase_idx];
  double *float_data = model_data->aero_phase_float_ptrs[aero_phase_idx];

  // Sum the mass and MW
  long double l_mass = 0.0;
  long double moles = 0.0;
  int i_jac = 0;
  for (int i_spec=0; i_spec<NUM_STATE_VAR_; i_spec++) {
    if (SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      l_mass += state_var[i_spec];
      moles += state_var[i_spec] / (long double) MW_(i_spec);
      if (jac_elem_mass) jac_elem_mass[i_jac] = 1.0L;
      if (jac_elem_MW)   jac_elem_MW[i_jac]   = 1.0L / MW_(i_spec);
      i_jac++;
    }
  }
  *MW = (double) l_mass / moles;
  if (jac_elem_MW) {
    for (int j_jac=0; j_jac<i_jac; j_jac++) {
      jac_elem_MW[j_jac] = (moles - jac_elem_MW[j_jac] * l_mass)
                           / (moles * moles);
    }
  }
  *mass = (double) l_mass;
}

/** \brief Get the volume of an aerosol phase
 *
 * \param model_data Pointer to the model data (state, env, aero_phase)
 * \param aero_phase_idx Index of the aerosol phase to use in the calculation
 * \param state_var Pointer to the aerosol phase on the state variable array
 * \param volume Pointer to hold the aerosol phase volume
 *               (\f$\mbox{\si{\cubic\metre\per\cubic\metre}}\f$ or
 *                \f$\mbox{\si{\cubic\metre\per particle}}\f$)
 * \param jac_elem When not NULL, a pointer to an array whose length is the
 *                 number of Jacobian elements used in calculations of mass and
 *                 volume of this aerosol phase returned by
 *                 \c aero_phase_get_used_jac_elem and whose contents will be
 *                 set to the partial derivatives of total phase volume by
 *                 concentration \f$\frac{dv}{dy_i}\f$ of each component
 *                 species \f$y_i\f$.
 */
void aero_phase_get_volume(ModelData *model_data, int aero_phase_idx,
          double *state_var, double *volume, double *jac_elem)
{

  // Get the requested aerosol phase data
  int *int_data      = model_data->aero_phase_int_ptrs[aero_phase_idx];
  double *float_data = model_data->aero_phase_float_ptrs[aero_phase_idx];

  // Sum the mass and MW
  *volume = 0.0;
  int i_jac = 0;
  for (int i_spec=0; i_spec<NUM_STATE_VAR_; i_spec++) {
    if (SPEC_TYPE_(i_spec)==CHEM_SPEC_VARIABLE ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_CONSTANT ||
        SPEC_TYPE_(i_spec)==CHEM_SPEC_PSSA) {
      *volume += state_var[i_spec] * 1.0e-9 / DENSITY_(i_spec);
      if (jac_elem) jac_elem[i_jac++] = 1.0e-9 / DENSITY_(i_spec);
    }
  }

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
  int *aero_phase_int_data      = model_data->nxt_aero_phase_int;
  double *aero_phase_float_data = model_data->nxt_aero_phase_float;

  // Save the pointers to this aerosol phases data
  model_data->aero_phase_int_ptrs[model_data->n_added_aero_phases] =
    aero_phase_int_data;
  model_data->aero_phase_float_ptrs[model_data->n_added_aero_phases] =
    aero_phase_float_data;
  ++(model_data->n_added_aero_phases);

  // Add the integer parameters
  for (; n_int_param>0; --n_int_param)
    *(aero_phase_int_data++) = *(int_param++);

  // Add the floating-point parameters
  for (; n_float_param>0; --n_float_param)
    *(aero_phase_float_data++) = (double) *(float_param++);

  // Set the pointers for the next free space in aero_phase_data
  model_data->nxt_aero_phase_int   = aero_phase_int_data;
  model_data->nxt_aero_phase_float = aero_phase_float_data;;
}

/** \brief Print the aerosol phase data
 * \param solver_data Pointer to the solver data
 */
void aero_phase_print_data(void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of aerosol phases
  int n_aero_phase = *(model_data->aero_phase_int_data);

  // Loop through the aerosol phases and print their data
  // advancing the aero_phase_data pointer each time
  for (int i_aero_phase=0; i_aero_phase<n_aero_phase; i_aero_phase++) {
    int *int_data      = model_data->aero_phase_int_ptrs[i_aero_phase];
    double *float_data = model_data->aero_phase_float_ptrs[i_aero_phase];

    printf("\n\nAerosol Phase %d\n\n", i_aero_phase);
  }
  fflush(stdout);
}

#undef NUM_STATE_VAR_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef SPEC_TYPE_
#undef MW_
#undef DENSITY_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_

