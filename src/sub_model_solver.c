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

// Sub model types (Must match parameters in pmc_sub_model_factory)
#define SUB_MODEL_UNIFAC 1

#ifdef PMC_USE_SUNDIALS

/** \brief Update sub model data for a new environmental state
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
void sub_model_update_env_state(ModelData *model_data, double *env)
{

  // Get the number of sub models
  int *sub_model_data = (int*) (model_data->sub_model_data);
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models to update the environmental conditions
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_UNIFAC :
        sub_model_data = (int*) sub_model_UNIFAC_update_env_state((void*) sub_model_data, env);
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
  int *sub_model_data = (int*) (model_data->sub_model_data);
  int n_sub_model = *(sub_model_data++);

  // Loop through the sub models to trigger their calculation
  // advancing the sub_model_data pointer each time
  for (int i_sub_model=0; i_sub_model<n_sub_model; i_sub_model++) {

    // Get the sub model type
    int sub_model_type = *(sub_model_data++);

    // Call the appropriate function
    switch (sub_model_type) {
      case SUB_MODEL_UNIFAC :
        sub_model_data = (int*) sub_model_UNIFAC_calculate((void*) sub_model_data, model_data);
        break;
    }
  }
}

/** \brief Print the sub model data
 * \param model_data Pointer to the model data
 */
void sub_model_print(ModelData *model_data)
{

  // Get the number of sub models
  int *sub_model_data = (int*) (model_data->sub_model_data);
  int n_sub_model = *(sub_model_data++);

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

#endif
