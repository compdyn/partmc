/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Sub model-specific functions for use by the solver
 */
/** \file
 * \brief Sub model solver functions
 */
#include <stdio.h>
#include <stdlib.h>
#include "sub_model_solver.h"
#include "sub_models.h"
#include <stdio.h>
#include <stdlib.h>

// Sub model types (Must match parameters in pmc_sub_model_factory)
#define SUB_MODEL_UNIFAC 1
#define SUB_MODEL_ZSR_AEROSOL_WATER 2
#define SUB_MODEL_PDFITE 3

/** \brief Get the Jacobian elements used by a particular sub model
 *
 * \param model_data A pointer to the model data
 * \param jac_struct A 2D array of flags indicating which Jacobian elements
 *                   may be used
 */
void sub_model_get_used_jac_elem(ModelData *model_data, bool **jac_struct)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct);
        break;
    }
  }

  // Account for sub-model interdependence
  int n_map_elem = 0;
  for (int i_dep=0; i_dep < model_data->n_per_cell_state_var; ++i_dep)
    for (int i_ind=0; i_ind < model_data->n_per_cell_state_var; ++i_ind)
      if (jac_struct[i_dep][i_ind] == true &&
          model_data->var_type[i_ind] != CHEM_SPEC_VARIABLE)
        for (int j_ind; j_ind < model_data->n_per_cell_state_var; ++j_ind)
          if (jac_struct[i_ind][j_ind]>=0) {
            jac_struct[i_dep][j_ind] = true;
            ++n_map_elem;
          }
  model_data->n_mapped_params = n_map_elem;

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
 * \param jac_ids Jacobian indices for sub-models
 */
void sub_model_set_jac_map(ModelData *model_data, int **jac_ids)
{

  // Allocate the map
  model_data->jac_map_params = malloc(sizeof(JacMap) *
                                      model_data->n_mapped_params);
  if (model_data->jac_map_params==NULL) {
    printf("\n\nError allocating sub model Jacobian map\n\n");
    EXIT_FAILURE;
  }

  // Set up a local structure array for individual sub-model Jacobian elements
  bool **jac_struct_local = malloc(sizeof(bool*) * model_data->n_per_cell_state_var);
  if (jac_struct_local==NULL) {
    printf("\n\nError allocating space for sub model Jac structure array\n\n");
    EXIT_FAILURE;
  }
  for (int i_var=0; i_var < model_data->n_per_cell_state_var; ++i_var) {
    jac_struct_local[i_var] = malloc(sizeof(bool) * model_data->n_per_cell_state_var);
    if (jac_struct_local[i_var]==NULL) {
      printf("\n\nError allocating space for sub model Jac structure array\n\n");
      EXIT_FAILURE;
    }
  }

  // Set up an index for map elements;
  int i_map = 0;

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Reset the local structure array
    for (int i_dep=0; i_dep < model_data->n_per_cell_state_var; ++i_dep)
      for (int i_ind=0; i_ind < model_data->n_per_cell_state_var; ++i_ind)
        jac_struct_local[i_dep][i_ind] = false;

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct_local);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct_local);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_get_used_jac_elem(
            sub_model_int_data, sub_model_float_data, jac_struct_local);
        break;
    }

    // Check for dependence on sub-model calculations and set mapping elements
    // if necessary
    for (int i_dep=0; i_dep < model_data->n_per_cell_state_var; ++i_dep)
      for (int i_ind=0; i_ind < model_data->n_per_cell_state_var; ++i_ind)
        if (jac_struct_local[i_dep][i_ind] == true &&
            model_data->var_type[i_ind] != CHEM_SPEC_VARIABLE) {
          for (int j_ind; j_ind < model_data->n_per_cell_state_var; ++j_ind) {
            if (jac_ids[i_ind][j_ind]>=0) {
              model_data->jac_map_params[i_map].solver_id = jac_ids[i_dep][j_ind];
              model_data->jac_map_params[i_map].rxn_id    = jac_ids[i_dep][i_ind];
              model_data->jac_map_params[i_map].param_id  = jac_ids[i_ind][j_ind];
              ++i_map;
            }
          }
        }

  }

  if (i_map != model_data->n_mapped_params) {
    printf("\n\nError mapping sub-model Jacobian elements %d %d\n\n",
           i_map, model_data->n_mapped_params);
    EXIT_FAILURE;
  }

  // free allocated memory
  for (int i_var=0; i_var < model_data->n_per_cell_state_var; ++i_var)
    free(jac_struct_local[i_var]);
  free(jac_struct_local);

}
/** \brief Update the time derivative and Jacobian array ids
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Ids for state variables on the time derivative array
 * \param jac_ids Ids for the state variables on the Jacobian array
 */
void sub_model_update_ids(ModelData *model_data, int *deriv_ids, int **jac_ids)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models and get their Jacobian elements
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_update_ids(
            sub_model_int_data, sub_model_float_data, deriv_ids, jac_ids);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_update_ids(
            sub_model_int_data, sub_model_float_data, deriv_ids, jac_ids);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_update_ids(
            sub_model_int_data, sub_model_float_data, deriv_ids, jac_ids);
        break;
    }
  }

  // Set the sub model interdependence Jacobian map
  sub_model_set_jac_map(model_data, jac_ids);
}

/** \brief Update sub model data for a new environmental state
 * \param model_data Pointer to the model data with updated env state
 */
void sub_model_update_env_state(ModelData *model_data)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to update the environmental conditions
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];
    double *sub_model_env_data = &(model_data->grid_cell_sub_model_env_data[
                                     model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_update_env_state(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_update_env_state(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
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
void sub_model_calculate(ModelData *model_data)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to trigger their calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];
    double *sub_model_env_data = &(model_data->grid_cell_sub_model_env_data[
                                     model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_calculate(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_calculate(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_calculate(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data);
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
    realtype time_step)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to trigger their Jacobian calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];
    double *sub_model_env_data = &(model_data->grid_cell_sub_model_env_data[
                                     model_data->sub_model_env_idx[i_sub_model]]);

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data, J_data, (double) time_step);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data, J_data, (double) time_step);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, sub_model_env_data,
                  model_data, J_data, (double) time_step);
        break;
    }
  }

  // Account for sub-model interdependence
  for (int i_map=0; i_map < model_data->n_mapped_params; ++i_map)
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
          int n_float_param, int n_env_param, int *int_param,
          double *float_param, void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *sub_model_int_data      = model_data->nxt_sub_model_int;
  double *sub_model_float_data = model_data->nxt_sub_model_float;
  int sub_model_env_idx        = model_data->nxt_sub_model_env;

  // Save the pointers to this sub model's data
  model_data->sub_model_int_ptrs[model_data->n_added_sub_models]   = sub_model_int_data;
  model_data->sub_model_float_ptrs[model_data->n_added_sub_models] = sub_model_float_data;
  model_data->sub_model_env_idx[model_data->n_added_sub_models]    = sub_model_env_idx;
  ++(model_data->n_added_sub_models);

  // Add the sub model type
  *(sub_model_int_data++) = sub_model_type;

  // Add integer parameters
  for (; n_int_param>0; --n_int_param) *(sub_model_int_data++) = *(int_param++);

  // Add floating-point parameters
  for (; n_float_param>0; --n_float_param)
          *(sub_model_float_data++) = *(float_param++);

  // Set the pointers for the next free space in the sub model data arrays
  model_data->nxt_sub_model_int     = sub_model_int_data;
  model_data->nxt_sub_model_float   = sub_model_float_data;
  model_data->nxt_sub_model_env     = sub_model_env_idx + n_env_param;
  model_data->n_sub_model_env_data += n_env_param;
}

/** \brief Update sub-model data
 *
 * Update data for one or more sub-models. Sub-models of a certain type are
 * passed a void pointer to updated data that must be in the format specified
 * by the sub-model type. This data could be used to find specific sub-models
 * of the specified type and change some model parameter(s).
 *
 * \param update_sub_model_type Type of the sub-model
 * \param update_data Pointer to updated data to pass to the sub-model
 * \param solver_data Pointer to solver data
 */
void sub_model_update_data(int update_sub_model_type, void *update_data,
          void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Skip sub-models of other types
    if (sub_model_type!=update_sub_model_type) continue;

    // Currently there are no sub-models with update data functions
  }
}

/** \brief Print the sub model data
 * \param solver_data Pointer to the solver data
 */
void sub_model_print_data(void *solver_data)
{
  ModelData *model_data = (ModelData*)
          &(((SolverData*)solver_data)->model_data);

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to print their data
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_print(sub_model_int_data,
            sub_model_float_data);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_print(sub_model_int_data,
            sub_model_float_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_print(
                  sub_model_int_data, sub_model_float_data);
        break;
    }
  }
  fflush(stdout);
}

