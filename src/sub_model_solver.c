/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Sub model-specific functions for use by the solver
 */
/** \file
 * \brief Sub model solver functions
 */
#include "sub_model_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include "sub_models.h"

// Sub model types (Must match parameters in pmc_sub_model_factory)
#define SUB_MODEL_UNIFAC 1
#define SUB_MODEL_ZSR_AEROSOL_WATER 2
#define SUB_MODEL_PDFITE 3

/** \brief Get the Jacobian elements used by a particular sub model
 *
 * \param model_data A pointer to the model data
 * \param jac Jacobian
 */
void sub_model_get_used_jac_elem(ModelData *model_data, Jacobian *jac) {
  // Set up local Jacobian to collect used elements
  Jacobian local_jac;
  if (jacobian_initialize_empty(
          &local_jac, (unsigned int)model_data->n_per_cell_state_var) != 1) {
    printf(
        "\n\nERROR allocating sub-model Jacobian structure for sub-model "
        "interdepenedence\n\n");
    exit(EXIT_FAILURE);
  }

  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_get_used_jac_elem(sub_model_int_data,
                                           sub_model_float_data, &local_jac);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_get_used_jac_elem(sub_model_int_data,
                                           sub_model_float_data, &local_jac);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, &local_jac);
        break;
    }
  }

  // Build the sparse Jacobian
  if (jacobian_build_matrix(&local_jac) != 1) {
    printf("\n\nERROR building sparse Jacobian for sub models\n\n");
    exit(EXIT_FAILURE);
  }

  // Add registered elements to sub-model Jacobian
  for (unsigned int i_ind = 0; i_ind < model_data->n_per_cell_state_var;
       ++i_ind) {
    for (unsigned int i_elem = jacobian_column_pointer_value(local_jac, i_ind);
         i_elem < jacobian_column_pointer_value(local_jac, i_ind + 1);
         ++i_elem) {
      unsigned int i_dep = jacobian_row_index(local_jac, i_elem);
      jacobian_register_element(jac, i_dep, i_ind);
    }
  }

  // Add elements for sub-model interdependence and save number of
  // interdepenedent elements
  model_data->n_mapped_params = 0;
  for (unsigned int i_ind = 0; i_ind < model_data->n_per_cell_state_var;
       ++i_ind) {
    for (unsigned int i_elem = jacobian_column_pointer_value(local_jac, i_ind);
         i_elem < jacobian_column_pointer_value(local_jac, i_ind + 1);
         ++i_elem) {
      unsigned int i_dep = jacobian_row_index(local_jac, i_elem);
      if (i_dep != i_ind && model_data->var_type[i_ind] != CHEM_SPEC_VARIABLE) {
        for (unsigned int j_ind = 0; j_ind < model_data->n_per_cell_state_var;
             ++j_ind) {
          if (jacobian_get_element_id(local_jac, i_ind, j_ind) != -1) {
            jacobian_register_element(jac, i_dep, j_ind);
            ++(model_data->n_mapped_params);
          }
        }
      }
    }
  }

  jacobian_free(&local_jac);
}

/** \brief Set the map for sub-model interdependence
 *
 * Uses the indices provided in \c jac_struct along with individual
 * calls to the sub-model \c get_used_jac_elem() functions to set up a map
 * to account for the dependence of sub model calculations on the results
 * of other sub model calculations. The sub model priorities used to
 * assemble the array of sub models, must result in independent sub-model
 * calculations appearing before dependent sub-model calculations in
 * the sub-model array.
 *
 * \param model_data Pointer to the model data
 * \param jac Jacobian
 */
void sub_model_set_jac_map(ModelData *model_data, Jacobian jac) {
  // Allocate the map
  model_data->jac_map_params =
      malloc(sizeof(JacMap) * model_data->n_mapped_params);
  if (model_data->jac_map_params == NULL) {
    printf("\n\nError allocating sub model Jacobian map\n\n");
    exit(EXIT_FAILURE);
  }

  // Set up an index for map elements;
  int i_map = 0;

  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    // Set up a local Jacobian for individual sub-model Jacobian elements
    Jacobian local_jac;
    if (jacobian_initialize_empty(
            &local_jac, (unsigned int)model_data->n_per_cell_state_var) != 1) {
      printf(
          "\n\nERROR allocating sub-model Jacobian structure for sub-model "
          "interdepenedence\n\n");
      exit(EXIT_FAILURE);
    }

    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_get_used_jac_elem(sub_model_int_data,
                                           sub_model_float_data, &local_jac);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_get_used_jac_elem(sub_model_int_data,
                                           sub_model_float_data, &local_jac);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, &local_jac);
        break;
    }

    // Build the Jacobian
    if (jacobian_build_matrix(&local_jac) != 1) {
      printf(
          "\n\nERROR building sub-model Jacobian structure for sub-model "
          "interdependence\n\n");
      exit(EXIT_FAILURE);
    }

    // Check for dependence on sub-model calculations and set mapping elements
    // if necessary
    for (unsigned int i_ind = 0; i_ind < model_data->n_per_cell_state_var;
         ++i_ind) {
      for (unsigned int i_elem =
               jacobian_column_pointer_value(local_jac, i_ind);
           i_elem < jacobian_column_pointer_value(local_jac, i_ind + 1);
           ++i_elem) {
        unsigned int i_dep = jacobian_row_index(local_jac, i_elem);
        if (i_dep != i_ind &&
            model_data->var_type[i_ind] != CHEM_SPEC_VARIABLE) {
          for (unsigned int j_ind = 0; j_ind < model_data->n_per_cell_state_var;
               ++j_ind) {
            if (jacobian_get_element_id(jac, i_ind, j_ind) != -1) {
              model_data->jac_map_params[i_map].solver_id =
                  jacobian_get_element_id(jac, i_dep, j_ind);
              model_data->jac_map_params[i_map].rxn_id =
                  jacobian_get_element_id(jac, i_dep, i_ind);
              model_data->jac_map_params[i_map].param_id =
                  jacobian_get_element_id(jac, i_ind, j_ind);
              ++i_map;
            }
          }
        }
      }
    }

    // free the local Jacobian
    jacobian_free(&local_jac);
  }

  if (i_map != model_data->n_mapped_params) {
    printf("\n\nError mapping sub-model Jacobian elements %d %d\n\n", i_map,
           model_data->n_mapped_params);
    exit(EXIT_FAILURE);
  }
}

/** \brief Update the time derivative and Jacobian array ids
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Ids for state variables on the time derivative array
 * \param jac Jacobian
 */
void sub_model_update_ids(ModelData *model_data, int *deriv_ids, Jacobian jac) {
  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_update_ids(sub_model_int_data, sub_model_float_data,
                                    deriv_ids, jac);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_update_ids(sub_model_int_data, sub_model_float_data,
                                    deriv_ids, jac);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_update_ids(
            sub_model_int_data, sub_model_float_data, deriv_ids, jac);
        break;
    }
  }

  // Set the sub model interdependence Jacobian map
  sub_model_set_jac_map(model_data, jac);
}

/** \brief Update sub model data for a new environmental state
 * \param model_data Pointer to the model data with updated env state
 */
void sub_model_update_env_state(ModelData *model_data) {
  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models to update the environmental conditions
  // advancing the sub_model_data pointer each time
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);
    double *sub_model_env_data =
        &(model_data->grid_cell_sub_model_env_data
              [model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_update_env_state(sub_model_int_data,
                                          sub_model_float_data,
                                          sub_model_env_data, model_data);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_update_env_state(sub_model_int_data,
                                          sub_model_float_data,
                                          sub_model_env_data, model_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_update_env_state(
            sub_model_int_data, sub_model_float_data, sub_model_env_data,
            model_data);
        break;
    }
  }
}

/** \brief Perform the sub model calculations for the current model state
 * \param model_data Pointer to the model data
 */
void sub_model_calculate(ModelData *model_data) {
  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models to trigger their calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);
    double *sub_model_env_data =
        &(model_data->grid_cell_sub_model_env_data
              [model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_calculate(sub_model_int_data, sub_model_float_data,
                                   sub_model_env_data, model_data);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_calculate(sub_model_int_data, sub_model_float_data,
                                   sub_model_env_data, model_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_calculate(sub_model_int_data,
                                              sub_model_float_data,
                                              sub_model_env_data, model_data);
        break;
    }
  }
}

/** \brief Calculate the Jacobian constributions from sub model calculations
 *
 * \param model_data Pointer to the model data
 * \param J_data Pointer to sub-model Jacobian data
 * \param time_step Current time step [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_get_jac_contrib(ModelData *model_data, double *J_data,
                               realtype time_step) {
  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models to trigger their Jacobian calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);
    double *sub_model_env_data =
        &(model_data->grid_cell_sub_model_env_data
              [model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_get_jac_contrib(
            sub_model_int_data, sub_model_float_data, sub_model_env_data,
            model_data, J_data, (double)time_step);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_get_jac_contrib(
            sub_model_int_data, sub_model_float_data, sub_model_env_data,
            model_data, J_data, (double)time_step);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_get_jac_contrib(
            sub_model_int_data, sub_model_float_data, sub_model_env_data,
            model_data, J_data, (double)time_step);
        break;
    }
  }

  // Account for sub-model interdependence
  for (int i_map = 0; i_map < model_data->n_mapped_params; ++i_map)
    J_data[model_data->jac_map_params[i_map].solver_id] +=
        J_data[model_data->jac_map_params[i_map].rxn_id] *
        J_data[model_data->jac_map_params[i_map].param_id];
}
#endif

/** \brief Add condensed data to the condensed data block for sub models
 *
 * \param sub_model_type Sub model type
 * \param n_int_param Number of integer parameters
 * \param n_float_param Number of floating-point parameters
 * \param n_env_param Number of environment-dependent parameters
 * \param int_param Pointer to integer parameter array
 * \param float_param Pointer to floating-point parameter array
 * \param solver_data Pointer to solver data
 */
void sub_model_add_condensed_data(int sub_model_type, int n_int_param,
                                  int n_float_param, int n_env_param,
                                  int *int_param, double *float_param,
                                  void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Get pointers to the sub model data
  int *sub_model_int_data = &(
      model_data->sub_model_int_data
          [model_data->sub_model_int_indices[model_data->n_added_sub_models]]);
  double *sub_model_float_data =
      &(model_data->sub_model_float_data[model_data->sub_model_float_indices
                                             [model_data->n_added_sub_models]]);

  // Save next indices by adding lengths
  model_data->sub_model_int_indices[model_data->n_added_sub_models + 1] =
      (n_int_param + 1) +
      model_data
          ->sub_model_int_indices[model_data->n_added_sub_models];  //+1 is type
  model_data->sub_model_float_indices[model_data->n_added_sub_models + 1] =
      n_float_param +
      model_data->sub_model_float_indices[model_data->n_added_sub_models];
  model_data->sub_model_env_idx[model_data->n_added_sub_models + 1] =
      model_data->sub_model_env_idx[model_data->n_added_sub_models] +
      n_env_param;
  ++(model_data->n_added_sub_models);

  // Add the sub model type
  *(sub_model_int_data++) = sub_model_type;

  // Add integer parameters
  for (; n_int_param > 0; --n_int_param)
    *(sub_model_int_data++) = *(int_param++);

  // Add floating-point parameters
  for (; n_float_param > 0; --n_float_param)
    *(sub_model_float_data++) = (double)*(float_param++);

  model_data->n_sub_model_env_data += n_env_param;
}

/** \brief Update sub-model data
 *
 * Update data for one or more sub-models. Sub-models of a certain type are
 * passed a void pointer to updated data that must be in the format specified
 * by the sub-model type. This data could be used to find specific sub-models
 * of the specified type and change some model parameter(s).
 *
 * \param cell_id Id of the grid cell to update
 * \param sub_model_id Index of the sub model to update (or 0 if unknown)
 * \param update_sub_model_type Type of the sub-model
 * \param update_data Pointer to updated data to pass to the sub-model
 * \param solver_data Pointer to solver data
 */
void sub_model_update_data(int cell_id, int *sub_model_id,
                           int update_sub_model_type, void *update_data,
                           void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Point to the environment-dependent data for the grid cell
  model_data->grid_cell_sub_model_env_data =
      &(model_data
            ->sub_model_env_data[cell_id * model_data->n_sub_model_env_data]);

  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models advancing the sub_model_data pointer each time
  for (; (*sub_model_id) < n_sub_model; (*sub_model_id)++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[*sub_model_id]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[*sub_model_id]]);
    double *sub_model_env_data = &(
        model_data->grid_cell_sub_model_env_data[model_data->sub_model_env_idx[(
            *sub_model_id)]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    bool found = false;

    // Skip sub-models of other types
    if (sub_model_type != update_sub_model_type) continue;

    // Currently there are no sub-models with update data functions
  }
}

/** \brief Print the sub model data
 * \param solver_data Pointer to the solver data
 */
void sub_model_print_data(void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Get the number of sub models
  int n_sub_model = model_data->n_sub_model;

  // Loop through the sub models to print their data
  // advancing the sub_model_data pointer each time
  for (int i_sub_model = 0; i_sub_model < n_sub_model; i_sub_model++) {
    int *sub_model_int_data =
        &(model_data->sub_model_int_data
              [model_data->sub_model_int_indices[i_sub_model]]);
    double *sub_model_float_data =
        &(model_data->sub_model_float_data
              [model_data->sub_model_float_indices[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE:
        sub_model_PDFiTE_print(sub_model_int_data, sub_model_float_data);
        break;
      case SUB_MODEL_UNIFAC:
        sub_model_UNIFAC_print(sub_model_int_data, sub_model_float_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER:
        sub_model_ZSR_aerosol_water_print(sub_model_int_data,
                                          sub_model_float_data);
        break;
    }
  }
  fflush(stdout);
}
