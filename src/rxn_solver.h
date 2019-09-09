/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for rxn_solver.c
 *
 */
/** \file
 * \brief Header file for abstract reaction functions
 */
#ifndef RXN_SOLVER_H
#define RXN_SOLVER_H
#include "camp_common.h"

/** Public reaction functions **/

/* Solver functions */
void rxn_get_used_jac_elem(ModelData *model_data, bool **jac_struct);
void rxn_update_ids(ModelData *model_data, int *deriv_ids, int **jac_ids);
void rxn_update_env_state(ModelData *model_data);
void rxn_reset_state_adjustments(ModelData *model_data);
void rxn_adjust_state(ModelData *model_data);
void rxn_print_data(void *solver_data);
#ifdef PMC_USE_SUNDIALS
void rxn_calc_deriv(ModelData *model_data, double *deriv_data, double time_step);
void rxn_calc_jac(ModelData *model_data, double *J_data, double time_step);
#endif

/* Setup functions */
void rxn_add_condensed_data(int rxn_type, int n_int_param,
	  int n_float_param, int n_env_param, int *int_param,
          double *float_param, void *solver_data);

/* Update data functions */
void rxn_update_data(int cell_id, int update_rxn_type, void *update_data,
          void *solver_data);
void rxn_free_update_data(void *update_data);

#endif
