/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Sub model-specific functions for use by the solver
 */
/** \file
 * \brief Sub model solver functions
 */
#include "phlex_solver.h"
#include "sub_model_solver.h"

#define NUM_ENV_VAR 2

// Sub model types (Must match parameters in pmc_sub_model_factory)
#define SUB_MODEL_UNIFAC 1


/** \brief Update the time derivative and Jacobian array ids
 *
 * \param model_data Model data
 * \param deriv_size Number of elements per state on the derivative array
 * \param jac_size Number of elements per state on the Jacobian array
 * \param deriv_ids Ids for state variables on the time derivative array
 * \param jac_ids Ids for state variables on the Jacobian array
 */
void sub_model_update_ids(ModelData model_data, int deriv_size, int jac_size, 
            int *deriv_ids, int **jac_ids)
{

  // Get a pointer to the sub model data
  int *sub_model_data = (int*) (model_data.sub_model_data);

  // Get the number of sub models
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub_model type
    int sub_model_type = *(sub_model_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_UNIFAC :
        sub_model_data = (int*) sub_model_UNIFAC_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) sub_model_data);
        break;

    }
  }
}

/** \brief Get a pointer to a calcualted sub model parameter
 * \param solver_data Pointer to the solver data
 * \param state_id Index of the unique state
 * \param sub_model_type Sub model type
 * \param identifiers Pointer to information needed by the sub-model to idenify
 *                    the parameter requested. These must be in the format
 *                    used by the sub model
 * \return Index of the requested parameter, or -1 if it was not found that can
 *         be passed to sub_model_get_parameter_value() during solving
 */
int sub_model_get_parameter_id_sd(void *solver_data, int state_id, 
          int sub_model_type, void *identifiers)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *model_data = &( sd->model_data[ state_id ] );
  return sub_model_get_parameter_id(*model_data, sub_model_type, identifiers);
}

/** \brief Get a pointer to a calcualted sub model parameter
 * \param model_data Model data
 * \param type Sub model type
 * \param identifiers Pointer to information needed by the sub-model to idenify
 *                    the parameter requested. These must be in the format
 *                    used by the sub model
 * \return Index of the requested parameter, or -1 if it was not found that can
 *         be passed to sub_model_get_parameter_value() during solving
 */
int sub_model_get_parameter_id(ModelData model_data, int type, 
          void *identifiers)
{
 
  // Get the number of sub models
  int *sub_model_data = (int*) (model_data.sub_model_data);
  int n_sub_model = *(sub_model_data++);

  // Initialize the parameter id
  int parameter_id = -1;
  int curr_id = 1;

  // Loop through the sub models to find the requested parameter advancing the
  // sub_model_data pointer each time. The first model found to provide the
  // parameter will be selected.
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);
    curr_id++;

    // Check if this is the requested type of sub model
    if (type==sub_model_type) {
      switch (sub_model_type) {
        case SUB_MODEL_UNIFAC :
          sub_model_data = (int*) sub_model_UNIFAC_get_parameter_id(
              (void*) sub_model_data, identifiers, &parameter_id);
          break;
      }
      // ... otherwise skip past it
    } else {
      switch (sub_model_type) {
        case SUB_MODEL_UNIFAC :
          sub_model_data = (int*) sub_model_UNIFAC_skip((void*) sub_model_data);
          break;
      }
    }

    // Check if the parameter was found
    if (parameter_id>=0) {
      parameter_id += curr_id;
      return parameter_id;
    // ... if not, advance the index to the next sub model
    } else {
      curr_id += (int) (sub_model_data - (int*) (model_data.sub_model_data));
    }
  }
  return -1;
}

/** \brief Return a parameter by its index in the sub model data block
 * \param solver_data Pointer to the solver data
 * \param state_id Index of unique state to get value for
 * \param parameter_id Index of the parameter in the data block
 * \return The parameter value
 */
PMC_C_FLOAT sub_model_get_parameter_value_sd(void *solver_data, int state_id,
          int parameter_id)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *model_data = &( sd->model_data[ state_id ] );
  return (PMC_C_FLOAT) sub_model_get_parameter_value( *model_data,
            parameter_id );
}

/** \brief Return a parameter by its index in the sub model data block
 * \param model_data Model data
 * \param parameter_id Index of the parameter in the data block
 * \return The parameter value
 */
PMC_C_FLOAT sub_model_get_parameter_value(ModelData model_data,
          int parameter_id)
{
  int *sub_model_data = (int*) (model_data.sub_model_data);
  sub_model_data += parameter_id;
  return *((PMC_C_FLOAT*) sub_model_data);
}

/** \brief Update sub model data for a new environmental state
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
void sub_model_update_env_state(ModelData model_data)
{
  PMC_C_FLOAT *env = model_data.env;
  int *sub_model_data = (int*) (model_data.sub_model_data);
  
  // Get the number of sub models
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models to update the environmental conditions
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_UNIFAC :
        sub_model_data = (int*) sub_model_UNIFAC_update_env_state(
                  (void*) sub_model_data, env);
        break;
    }
  }
}

/** \brief Perform the sub model calculations for the current model state
 * \param model_data Model data
 */
void sub_model_calculate(ModelData model_data)
{

  int *sub_model_data = (int*) (model_data.sub_model_data);
  
  // Get the number of sub models
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models to trigger their calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_UNIFAC :
        sub_model_data = (int*) sub_model_UNIFAC_calculate(
                  (void*) sub_model_data, model_data);
        break;
    }
  }
}

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
          int n_float_param, int *int_param, PMC_C_FLOAT *float_param,
          void *solver_data)
{

  SolverData * sd = (SolverData*) solver_data;

  // Add condensed data for each state's ModelData object
  for( int i_state = 0; i_state < sd->n_states; i_state++ ) {

    ModelData *md = &( sd->model_data[ i_state ] );
    int *sub_model_data = (int*) (md->nxt_sub_model);

    // Add the sub model type
    *(sub_model_data++) = sub_model_type;

    // Add integer parameters
    for (int i=0; i<n_int_param; i++) *(sub_model_data++) = int_param[i];

    // Add floating-point parameters
    PMC_C_FLOAT *flt_ptr = (PMC_C_FLOAT*) sub_model_data;
    for (int i=0; i<n_float_param; i++)
            *(flt_ptr++) = (PMC_C_FLOAT) float_param[i];

    // Set the pointer for the next free space in sub_model_data
    md->nxt_sub_model = (void*) flt_ptr;

  }
}

/** \brief Update sub-model data
 *
 * Update data for one or more sub-models. Sub-models of a certain type are
 * passed a void pointer to updated data that must be in the format specified
 * by the sub-model type. This data could be used to find specific sub-models
 * of the specified type and change some model parameter(s).
 *
 * \param state_id Index of unique state to update
 * \param update_sub_model_type Type of the sub-model
 * \param update_data Pointer to updated data to pass to the sub-model
 * \param solver_data Pointer to solver data
 */
void sub_model_update_data(int state_id, int update_sub_model_type,
          void *update_data, void *solver_data)
{
  ModelData *model_data = 
          (ModelData*) &(((SolverData*)solver_data)->model_data[state_id]);

  // Get the number of sub models
  int *sub_model_data = (int*) (model_data->sub_model_data);
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);

    // Skip sub-models of other types
    if (sub_model_type!=update_sub_model_type) {
      switch (sub_model_type) {
        case SUB_MODEL_UNIFAC :
          sub_model_data = (int*) sub_model_UNIFAC_skip((void*) sub_model_data);
          break;
      }

    // ... otherwise, call the update function for sub-model types that have 
    // then
    } else {
      switch (sub_model_type) {
        case SUB_MODEL_UNIFAC :
          sub_model_data = (int*) sub_model_UNIFAC_skip((void*) sub_model_data);
          break;
      }
    }
  }
}

/** \brief Print the sub model data
 * \param model_data Pointer to the model data
 */
void sub_model_print_data(void *solver_data)
{

  SolverData *sd = (SolverData*) solver_data;

  // Loop through the unique states to solve
  for (int i_state = 0; i_state < sd->n_states; i_state++) {
  
    ModelData *model_data = &( sd->model_data[ i_state ] );
    int *sub_model_data = (int*) (model_data->sub_model_data);
  
    // Get the number of sub models
    int n_sub_model = *(sub_model_data++);

    printf("\nSub-Models in state %d\n", i_state);

    // Loop through the sub models to print their data
    // advancing the sub_model_data pointer each time
    for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

      // Get the sub model type
      int sub_model_type = *(sub_model_data++);

      // Call the appropriate function
      switch (sub_model_type) {
        case SUB_MODEL_UNIFAC :
          sub_model_data = (int*) sub_model_UNIFAC_print((void*) sub_model_data);
          break;
      }
    }
  }
}

