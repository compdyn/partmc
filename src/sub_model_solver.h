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
#include "phlex_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef PMC_USE_SUNDIALS

// UNIFAC
void * sub_model_UNIFAC_get_used_jac_elem(void *sub_model_data, bool *jac_row);
void * sub_model_UNIFAC_update_ids(void *sub_model_data, int *jac_row);
void * sub_model_UNIFAC_get_parameter_id(void *sub_model_data, void* identifiers, int *parameter_id);
void * sub_model_UNIFAC_update_env_state(void *sub_model_data, realtype *env_data);
void * sub_model_UNIFAC_calculate(void *sub_model_data, ModelData *model_data);
void * sub_model_UNIFAC_add_jac_contrib(void *sub_model_data, realtype base_val, realtype *jac_row);
void * sub_model_UNIFAC_skip(void *sub_model_data);
void * sub_model_UNIFAC_print(void *sub_model_data);

#endif
#endif
