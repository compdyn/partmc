/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for sub_model_solver.c
 *
 */
/** \file
 * \brief Header file for abstract sub model functions
 */
#ifndef SUB_MODEL_SOLVER_H
#define SUB_MODEL_SOLVER_H
#include "phlex_common.h"

/** Public sub model functions **/

/* Solver functions */
void sub_model_get_used_jac_elem(ModelData *model_data, bool **jac_struct);
void sub_model_update_ids(ModelData *model_data, int *deriv_ids, int **jac_ids);
void sub_model_update_env_state(ModelData *model_data, double *env);
int sub_model_get_parameter_id(ModelData *model_data, int type,
          void *identifiers);
double sub_model_get_parameter_value(ModelData *model_data, int parameter_id);
void sub_model_calculate(ModelData *model_data);
void sub_model_print_data(void *solver_data);

/* Setup functions */
void sub_model_add_condensed_data(int sub_model_type, int n_int_param,
	  int n_float_param, int *int_param, double *float_param,
          void *solver_data);
void sub_model_update_data(int update_sub_model_type, void *update_data,
          void *solver_data);

/* Update data functions */
int sub_model_get_parameter_id_sd(void *solver_data, int sub_model_type,
          void *identifiers);
double sub_model_get_parameter_value_sd(void *solver_data, int parameter_id);

#endif
