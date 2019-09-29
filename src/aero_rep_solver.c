/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aerosol representation-specific functions for use by the solver
 *
 */
/** \file
 * \brief Aerosol representation functions
 */
#include "aero_rep_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aero_reps.h"

// Aerosol representations (Must match parameters defined in
// pmc_aero_rep_factory
#define AERO_REP_SINGLE_PARTICLE 1
#define AERO_REP_MODAL_BINNED_MASS 2

/** \brief Flag Jacobian elements used to calculated mass, volume, etc.
 *
 * \param model_data A pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param jac_struct 1D array of flags indicating potentially non-zero
 *                   Jacobian elements. (The dependent variable should have
 *                   been chosen by the calling function.)
 * \return Number of Jacobian elements flagged
 */
int aero_rep_get_used_jac_elem(ModelData *model_data, int aero_rep_idx,
                               int aero_phase_idx, bool *jac_struct) {
  int num_flagged_elem = 0;

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle radius and set of partial derivatives
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      num_flagged_elem = aero_rep_modal_binned_mass_get_used_jac_elem(
          model_data, aero_phase_idx, aero_rep_int_data, aero_rep_float_data,
          jac_struct);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      num_flagged_elem = aero_rep_single_particle_get_used_jac_elem(
          model_data, aero_phase_idx, aero_rep_int_data, aero_rep_float_data,
          jac_struct);
      break;
  }

  return num_flagged_elem;
}
/** \brief Get state array elements used by aerosol representation functions
 *
 * \param model_data A pointer to the model data
 * \param state_flags An array of flags the length of the state array
 *                    indicating species used
 */
void aero_rep_get_dependencies(ModelData *model_data, bool *state_flags) {
  // Get the number of aerosol representations
  int n_aero_rep = model_data->n_aero_rep;

  // Loop through the aerosol representations to determine the Jacobian elements
  // used advancing the aero_rep_data pointer each time
  for (int i_aero_rep = 0; i_aero_rep < n_aero_rep; i_aero_rep++) {
    // Get pointers to the aerosol data
    int *aero_rep_int_data =
            &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[i_aero_rep]]);
    double *aero_rep_float_data =
            &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[i_aero_rep]]);

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_int_data++);

    // Call the appropriate function
    switch (aero_rep_type) {
      case AERO_REP_MODAL_BINNED_MASS:
        aero_rep_modal_binned_mass_get_dependencies(
            aero_rep_int_data, aero_rep_float_data, state_flags);
        break;
      case AERO_REP_SINGLE_PARTICLE:
        aero_rep_single_particle_get_dependencies(
            aero_rep_int_data, aero_rep_float_data, state_flags);
        break;
    }
  }
}

/** \brief Update the aerosol representations for new environmental conditions
 *
 * \param model_data Pointer to the model data
 */
void aero_rep_update_env_state(ModelData *model_data) {
  // Get the number of aerosol representations
  int n_aero_rep = model_data->n_aero_rep;

  // Loop through the aerosol representations to update the environmental
  // conditions advancing the aero_rep_data pointer each time
  for (int i_aero_rep = 0; i_aero_rep < n_aero_rep; i_aero_rep++) {
    // Get pointers to the aerosol data
    int *aero_rep_int_data =
            &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[i_aero_rep]]);
    double *aero_rep_float_data =
            &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[i_aero_rep]]);
    double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                      model_data->aero_rep_env_idx[i_aero_rep]]);

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_int_data++);

    // Call the appropriate function
    switch (aero_rep_type) {
      case AERO_REP_MODAL_BINNED_MASS:
        aero_rep_modal_binned_mass_update_env_state(
            model_data, aero_rep_int_data, aero_rep_float_data,
            aero_rep_env_data);
        break;
      case AERO_REP_SINGLE_PARTICLE:
        aero_rep_single_particle_update_env_state(model_data, aero_rep_int_data,
                                                  aero_rep_float_data,
                                                  aero_rep_env_data);
        break;
    }
  }
}

/** \brief Update the aerosol representations for a new state
 *
 * \param model_data Pointer to the model data
 */
void aero_rep_update_state(ModelData *model_data) {
  // Get the number of aerosol representations
  int n_aero_rep = model_data->n_aero_rep;

  // Loop through the aerosol representations to update the state
  // advancing the aero_rep_data pointer each time
  for (int i_aero_rep = 0; i_aero_rep < n_aero_rep; i_aero_rep++) {
    // Get pointers to the aerosol data
    int *aero_rep_int_data =
            &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[i_aero_rep]]);
    double *aero_rep_float_data =
            &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[i_aero_rep]]);
    double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                      model_data->aero_rep_env_idx[i_aero_rep]]);

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_int_data++);

    // Call the appropriate function
    switch (aero_rep_type) {
      case AERO_REP_MODAL_BINNED_MASS:
        aero_rep_modal_binned_mass_update_state(model_data, aero_rep_int_data,
                                                aero_rep_float_data,
                                                aero_rep_env_data);
        break;
      case AERO_REP_SINGLE_PARTICLE:
        aero_rep_single_particle_update_state(model_data, aero_rep_int_data,
                                              aero_rep_float_data,
                                              aero_rep_env_data);
        break;
    }
  }
}

/** \brief Get the effective particle radius, \f$r_{eff}\f$ (m)
 *
 * Calculates effective particle radius \f$r_{eff}\f$ (m), as well as the set of
 * \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param radius Pointer to hold effective particle radius (m)
 * \param partial_deriv Pointer to the set of partial derivatives to be
 *                      calculated \f$\frac{\partial r_{eff}}{\partial y}\f$,
 *                      or a NULL pointer if no partial derivatives are needed
 */
void aero_rep_get_effective_radius(ModelData *model_data, int aero_rep_idx,
                                   int aero_phase_idx, double *radius,
                                   double *partial_deriv) {
  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                    model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle radius and set of partial derivatives
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      aero_rep_modal_binned_mass_get_effective_radius(
          model_data, aero_phase_idx, radius, partial_deriv, aero_rep_int_data,
          aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      aero_rep_single_particle_get_effective_radius(
          model_data, aero_phase_idx, radius, partial_deriv, aero_rep_int_data,
          aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return;
}

/** \brief Get the particle number concentration \f$n\f$
 * (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * Calculates particle number concentration, \f$n\f$
 * (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$), as well as the set of
 * \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param number_conc Pointer to hold calculated number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \param partial_deriv Pointer to the set of partial derivatives to be
 *                      calculated \f$\frac{\partial n}{\partial y}\f$, or a
 *                      NULL pointer if no partial derivatives are required
 */
void aero_rep_get_number_conc(ModelData *model_data, int aero_rep_idx,
                              int aero_phase_idx, double *number_conc,
                              double *partial_deriv) {
  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                    model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      aero_rep_modal_binned_mass_get_number_conc(
          model_data, aero_phase_idx, number_conc, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      aero_rep_single_particle_get_number_conc(
          model_data, aero_phase_idx, number_conc, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return;
}

/** \brief Check whether aerosol concentrations are per-particle or total for
 * each phase
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \return 0 for per-particle; 1 for total for each phase
 */
int aero_rep_get_aero_conc_type(ModelData *model_data, int aero_rep_idx,
                                int aero_phase_idx) {
  // Initialize the aerosol concentration type
  int aero_conc_type = 0;

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                    model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the type of aerosol concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      aero_rep_modal_binned_mass_get_aero_conc_type(
          aero_phase_idx, &aero_conc_type, aero_rep_int_data,
          aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      aero_rep_single_particle_get_aero_conc_type(
          aero_phase_idx, &aero_conc_type, aero_rep_int_data,
          aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return aero_conc_type;
}

/** \brief Get the total mass of an aerosol phase in this representation
 **        \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * Calculates total aerosol phase mass, \f$m\f$
 * (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$), as well as the set of
 * \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param aero_phase_mass Pointer to hold calculated aerosol-phase mass,
 *                        \f$m\f$
 *                        (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 * \param partial_deriv Pointer to the set of partial derivatives to be
 *                      calculated \f$\frac{\partial m}{\partial y}\f$, or a
 *                      NULL pointer if no partial derivatives are needed
 */
void aero_rep_get_aero_phase_mass(ModelData *model_data, int aero_rep_idx,
                                  int aero_phase_idx, double *aero_phase_mass,
                                  double *partial_deriv) {
  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                    model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      aero_rep_modal_binned_mass_get_aero_phase_mass(
          model_data, aero_phase_idx, aero_phase_mass, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      aero_rep_single_particle_get_aero_phase_mass(
          model_data, aero_phase_idx, aero_phase_mass, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
  }
}

/** \brief Get the average molecular weight of an aerosol phase in this
 **        representation \f$m\f$
 *(\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * Calculates total aerosol phase mass, \f$m\f$
 * (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$), as well as the set of
 * \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param aero_phase_avg_MW Pointer to hold calculated average MW in the
 *                          aerosol phase (\f$\mbox{\si{\kilogram\per\mole}}\f$)
 * \param partial_deriv Pointer to the set of partial derivatives to be
 *                      calculated \f$\frac{\partial m}{\partial y}\f$, or a
 *                      NULL pointer if no partial derivatives are needed
 */
void aero_rep_get_aero_phase_avg_MW(ModelData *model_data, int aero_rep_idx,
                                    int aero_phase_idx,
                                    double *aero_phase_avg_MW,
                                    double *partial_deriv) {
  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                    model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS:
      aero_rep_modal_binned_mass_get_aero_phase_avg_MW(
          model_data, aero_phase_idx, aero_phase_avg_MW, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE:
      aero_rep_single_particle_get_aero_phase_avg_MW(
          model_data, aero_phase_idx, aero_phase_avg_MW, partial_deriv,
          aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
  }
}

/** \brief Add condensed data to the condensed data block for aerosol
 * representations
 *
 * \param aero_rep_type Aerosol representation type
 * \param n_int_param Number of integer parameters
 * \param n_float_param Number of floating-point parameters
 * \param n_env_param Number of environment-dependent parameters
 * \param int_param Pointer to integer parameter array
 * \param float_param Pointer to floating-point parameter array
 * \param solver_data Pointer to solver data
 */
void aero_rep_add_condensed_data(int aero_rep_type, int n_int_param,
          int n_float_param, int n_env_param, int *int_param,
          double *float_param, void *solver_data)
{
  ModelData *model_data = (ModelData*)
          &(((SolverData*)solver_data)->model_data);

  // Get pointers to the aerosol representation data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[model_data->n_added_aero_reps]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[model_data->n_added_aero_reps]]);

  //Save next indices by adding lengths
  model_data->aero_rep_int_indices[model_data->n_added_aero_reps+1] =
          (n_int_param+1) + model_data->aero_rep_int_indices[model_data->n_added_aero_reps];//+1 is type
  model_data->aero_rep_float_indices[model_data->n_added_aero_reps+1] =
          n_float_param + model_data->aero_rep_float_indices[model_data->n_added_aero_reps];
  model_data->aero_rep_env_idx[model_data->n_added_aero_reps+1] =
          model_data->aero_rep_env_idx[model_data->n_added_aero_reps] + n_env_param;
  ++(model_data->n_added_aero_reps);

  // Add the reaction type
  *(aero_rep_int_data++) = aero_rep_type;

  // Add integer parameters
  for (; n_int_param>0; --n_int_param) *(aero_rep_int_data++) = *(int_param++);

  // Add floating-point parameters
  for (; n_float_param>0; --n_float_param)
    *(aero_rep_float_data++) = (double) *(float_param++);

  model_data->n_aero_rep_env_data += n_env_param;

}

/** \brief Update aerosol representation data
 *
 * \param cell_id Id of the grid cell to update
 * \param aero_rep_id Aerosol representation id (or 0 if unknown)
 * \param update_aero_rep_type Aerosol representation type to update
 * \param update_data Pointer to data needed for update
 * \param solver_data Pointer to solver data
 */
void aero_rep_update_data(int cell_id, int *aero_rep_id,
                          int update_aero_rep_type, void *update_data,
                          void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Point to the environment-dependent data for the grid cell
  model_data->grid_cell_aero_rep_env_data = &(
      model_data->aero_rep_env_data[cell_id * model_data->n_aero_rep_env_data]);

  // Get the number of aerosol representations
  int n_aero_rep = model_data->n_aero_rep;

  // Loop through the aerosol representations advancing the pointer each time
  for (; (*aero_rep_id) < n_aero_rep; (*aero_rep_id)++) {
    // Get pointers to the aerosol data
    int *aero_rep_int_data =
            &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[*aero_rep_id]]);
    double *aero_rep_float_data =
            &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[*aero_rep_id]]);
    double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
                                      model_data->aero_rep_env_idx[*aero_rep_id]]);

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_int_data++);

    bool found = false;

    // Find aerosol representations of the correct type
    if (aero_rep_type == update_aero_rep_type) {
      switch (aero_rep_type) {
        case AERO_REP_MODAL_BINNED_MASS:
          found = aero_rep_modal_binned_mass_update_data(
              (void *)update_data, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
          break;
        case AERO_REP_SINGLE_PARTICLE:
          found = aero_rep_single_particle_update_data(
              (void *)update_data, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
          break;
      }
      if (found) return;
    }
  }
}

/** \brief Print the aerosol representation data
 *
 * \param solver_data Pointer to the solver data
 */
void aero_rep_print_data(void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Get the number of aerosol representations
  int n_aero_rep = model_data->n_aero_rep;

  printf(
      "\n\nAerosol representation data\n\nnumber of aerosol "
      "representations: %d\n\n",
      n_aero_rep);

  // Loop through the aerosol representations advancing the pointer each time
  for (int i_aero_rep = 0; i_aero_rep < n_aero_rep; i_aero_rep++) {
    // Get pointers to the aerosol data
    int *aero_rep_int_data =
            &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[i_aero_rep]]);
    double *aero_rep_float_data =
            &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[i_aero_rep]]);

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_int_data++);

    // Call the appropriate printing function
    switch (aero_rep_type) {
      case AERO_REP_MODAL_BINNED_MASS:
        aero_rep_modal_binned_mass_print(aero_rep_int_data,
                                         aero_rep_float_data);
        break;
      case AERO_REP_SINGLE_PARTICLE:
        aero_rep_single_particle_print(aero_rep_int_data, aero_rep_float_data);
        break;
    }
  }
  fflush(stdout);
}

/** \brief Free an update data object
 *
 * \param update_data Object to free
 */
void aero_rep_free_update_data(void *update_data) { free(update_data); }
