/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for sub model calculations
 */
/** \file
 * \brief Header file for sub model functions
 */
#ifndef SUB_MODEL_SOLVER_H_
#define SUB_MODEL_SOLVER_H_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "phlex_solver.h"

// UNIFAC
void * sub_model_UNIFAC_get_used_jac_elem(
          void *sub_model_data, pmc_bool *jac_row);
void * sub_model_UNIFAC_update_ids(
          ModelData model_data, int *deriv_ids, int **jac_ids,
          void *sub_model_data);
void * sub_model_UNIFAC_get_parameter_id(
          void *sub_model_data, void* identifiers, int *parameter_id);
void * sub_model_UNIFAC_update_env_state(
          void *sub_model_data, PMC_C_FLOAT *env_data);
void * sub_model_UNIFAC_calculate(
          void *sub_model_data, ModelData model_data);
void * sub_model_UNIFAC_add_jac_contrib(
          void *sub_model_data, PMC_C_FLOAT base_val, PMC_C_FLOAT *jac_row);
void * sub_model_UNIFAC_skip(
          void *sub_model_data);
void * sub_model_UNIFAC_print(
          void *sub_model_data);

#endif
