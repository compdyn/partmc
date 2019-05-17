/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Sub model-specific functions for use by the solver
 */
/** \file
 * \brief Sub model solver functions
 */
extern "C" {
#include "sub_model_solver_gpu.h"
//#include "sub_models_gpu.h"

// Sub model types (Must match parameters in pmc_sub_model_gpu_factory)
#define SUB_MODEL_UNIFAC 1

/** \brief Return a parameter by its index in the sub model data block
 * \param model_data Pointer to the model data
 * \param parameter_id Index of the parameter in the data block
 * \return The parameter value
 */
__device__ double sub_model_gpu_get_parameter_value(ModelDatagpu *model_data, int parameter_id) {
  int *sub_model_gpu_data = (int *) (model_data->sub_model_gpu_data);
  sub_model_gpu_data += parameter_id;
  return *((double *) sub_model_gpu_data);
}

}