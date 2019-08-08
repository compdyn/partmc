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
}

/** \brief Get a pointer to a calcualted sub model parameter
 * \param solver_data Pointer to the solver data
 * \param sub_model_type Sub model type
 * \param identifiers Pointer to information needed by the sub-model to idenify
 *                    the parameter requested. These must be in the format
 *                    used by the sub model
 * \return Index of the requested parameter, or -1 if it was not found that can
 *         be passed to sub_model_get_parameter_value() during solving
 */
int sub_model_get_parameter_id_sd(void *solver_data, int sub_model_type,
    void *identifiers)
{
  ModelData *model_data = &(((SolverData*)solver_data)->model_data);
  return sub_model_get_parameter_id(model_data, sub_model_type, identifiers);
}

/** \brief Get a pointer to a calcualted sub model parameter
 * \param model_data Pointer to the model data
 * \param type Sub model type
 * \param identifiers Pointer to information needed by the sub-model to idenify
 *                    the parameter requested. These must be in the format
 *                    used by the sub model
 * \return Index of the requested parameter, or -1 if it was not found that can
 *         be passed to sub_model_get_parameter_value() during solving
 */
int sub_model_get_parameter_id(ModelData *model_data, int type,
          void *identifiers)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Initialize the parameter id
  int parameter_id = -1;

  // Loop through the sub models to find the requested parameter.
  // The first model found to provide the parameter will be selected.
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Check if this is the requested type of sub model
    if (type==sub_model_type) {
      switch (sub_model_type) {
        case SUB_MODEL_PDFITE :
          sub_model_PDFiTE_get_parameter_id(
              sub_model_int_data, sub_model_float_data, identifiers, &parameter_id);
          break;
        case SUB_MODEL_UNIFAC :
          sub_model_UNIFAC_get_parameter_id(
              sub_model_int_data, sub_model_float_data, identifiers, &parameter_id);
          break;
        case SUB_MODEL_ZSR_AEROSOL_WATER :
          sub_model_ZSR_aerosol_water_get_parameter_id(
              sub_model_int_data, sub_model_float_data, identifiers, &parameter_id);
          break;
      }
    }

    // Check if the parameter was found
    if (parameter_id>=0) {
      parameter_id += (sub_model_float_data - model_data->sub_model_float_data);
      return parameter_id;
    }

  }
  return -1;
}

/** \brief Return a parameter by its index in the sub model data block
 * \param solver_data Pointer to the solver data
 * \param parameter_id Index of the parameter in the data block
 * \return The parameter value
 */
double sub_model_get_parameter_value_sd(void *solver_data, int parameter_id)
{
  ModelData *model_data = &(((SolverData*)solver_data)->model_data);
  return (double) sub_model_get_parameter_value(model_data, parameter_id);
}

/** \brief Return a parameter by its index in the sub model data block
 * \param model_data Pointer to the model data
 * \param parameter_id Index of the parameter in the data block
 * \return The parameter value
 */
double sub_model_get_parameter_value(ModelData *model_data, int parameter_id)
{
  double *param = model_data->sub_model_float_data;
  param += parameter_id;
  return *param;
}

/** \brief Update sub model data for a new environmental state
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
void sub_model_update_env_state(ModelData *model_data, double *env)
{

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to update the environmental conditions
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_update_env_state(
                  sub_model_int_data, sub_model_float_data, env);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_update_env_state(
                  sub_model_int_data, sub_model_float_data, env);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_update_env_state(
                  sub_model_int_data, sub_model_float_data, env);
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

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_calculate(
                  sub_model_int_data, sub_model_float_data, model_data);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_calculate(
                  sub_model_int_data, sub_model_float_data, model_data);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_calculate(
                  sub_model_int_data, sub_model_float_data, model_data);
        break;
    }
  }
}

/** \brief Calculate the Jacobian constributions from sub model calculations
 *
 * \param model_data Pointer to the model data
 * \param J Sub-model Jacobian
 * \param time_step Current time step [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_get_jac_contrib(ModelData *model_data, SUNMatrix J,
    realtype time_step)
{

  // Get a pointer to the Jacobian data
  realtype *J_data = SM_DATA_S(J);

  // Get the number of sub models
  int n_sub_model = model_data->sub_model_int_data[0];

  // Loop through the sub models to trigger their Jacobian calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    int *sub_model_int_data = model_data->sub_model_int_ptrs[i_sub_model];
    double *sub_model_float_data = model_data->sub_model_float_ptrs[i_sub_model];

    // Get the sub model type
    int sub_model_type = *(sub_model_int_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        sub_model_PDFiTE_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, model_data, J_data,
                  (double) time_step);
        break;
      case SUB_MODEL_UNIFAC :
        sub_model_UNIFAC_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, model_data, J_data,
                  (double) time_step);
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        sub_model_ZSR_aerosol_water_get_jac_contrib(
                  sub_model_int_data, sub_model_float_data, model_data, J_data,
                  (double) time_step);
        break;
    }
  }
}
#endif

/** \brief Add condensed data to the condensed data block for sub models
 *
 * \param sub_model_type Sub model type
 * \param n_int_param Number of integer parameters
 * \param n_float_param Number of floating-point parameters
 * \param int_param Pointer to integer parameter array
 * \param float_param Pointer to floating-point parameter array
 * \param solver_data Pointer to solver data
 */
void sub_model_add_condensed_data(int sub_model_type, int n_int_param,
          int n_float_param, int *int_param, double *float_param,
          void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *sub_model_int_data      = model_data->nxt_sub_model_int;
  double *sub_model_float_data = model_data->nxt_sub_model_float;

  // Save the pointers to this sub model's data
  model_data->sub_model_int_ptrs[model_data->n_added_sub_models] = sub_model_int_data;
  model_data->sub_model_float_ptrs[model_data->n_added_sub_models] = sub_model_float_data;
  ++(model_data->n_added_sub_models);

  // Add the sub model type
  *(sub_model_int_data++) = sub_model_type;

  // Add integer parameters
  for (; n_int_param>0; --n_int_param) *(sub_model_int_data++) = *(int_param++);

  // Add floating-point parameters
  for (; n_float_param>0; --n_float_param)
          *(sub_model_float_data++) = *(float_param++);

  // Set the pointers for the next free space in the sub model data arrays
  model_data->nxt_sub_model_int   = sub_model_int_data;
  model_data->nxt_sub_model_float = sub_model_float_data;
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

    // ... otherwise, call the update function for sub-model types that have
    // then
    switch (sub_model_type) {
      case SUB_MODEL_PDFITE :
        break;
      case SUB_MODEL_UNIFAC :
        break;
      case SUB_MODEL_ZSR_AEROSOL_WATER :
        break;
    }
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

