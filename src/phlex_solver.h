/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for solver functions
 *
*/
/** \file
 * \brief Header file for solver functions
*/
#ifndef PHLEX_SOLVER_H_
#define PHLEX_SOLVER_H_
#include "phlex_common.h"

/* Functions called by phlex-chem */
void * solver_new(int n_state_var, int *var_type, int n_rxn,
          int n_rxn_int_param, int n_rxn_float_param, int n_aero_phase,
          int n_aero_phase_int_param, int n_aero_phase_float_param,
          int n_aero_rep, int n_aero_rep_int_param, int n_aero_rep_float_param,
          int n_sub_model, int n_sub_model_int_param,
          int n_sub_model_float_param);
void solver_initialize(void *solver_data, double *abs_tol, double rel_tol,
          int max_steps, int max_conv_fails);
int solver_run(void *solver_data, double *state, double *env, double t_initial,
	  double t_final);
void solver_free(void *solver_data);
void model_free(ModelData model_data);

#ifdef PMC_USE_SUNDIALS
/* Functions called by the solver */
int f(realtype t, N_Vector y, N_Vector deriv, void *model_data);
int Jac(realtype t, N_Vector y, N_Vector deriv, SUNMatrix J, void *model_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* SUNDIALS support functions */
int phlex_solver_update_model_state(N_Vector solver_state,
          ModelData *model_data);
int phlex_solver_update_solver_state(N_Vector solver_state,
          ModelData *model_data);
SUNMatrix get_jac_init(SolverData *solver_data);
int check_flag(void *flag_value, char *func_name, int opt);
void check_flag_fail(void *flag_value, char *func_name, int opt);
static void solver_print_stats(void *cvode_mem);
#endif

#endif
