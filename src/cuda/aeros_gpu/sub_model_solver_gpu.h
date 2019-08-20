/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for sub_model_gpu_solver.c
 *
 */
/** \file
 * \brief Header file for abstract sub model functions
 */
#ifndef SUB_MODEL_SOLVER_H
#define SUB_MODEL_SOLVER_H
#include "../phlex_gpu_solver.h"

__device__ double sub_model_gpu_get_parameter_value(ModelData *model_data, int parameter_id);

#endif
