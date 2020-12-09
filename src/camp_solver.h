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
#ifndef CAMP_SOLVER_H_
#define CAMP_SOLVER_H_
#include "camp_common.h"

/* Functions called by camp-chem */
void *solver_new(int n_state_var, int n_cells, int *var_type, int n_rxn,
                 int n_rxn_int_param, int n_rxn_float_param,
                 int n_rxn_env_param, int n_aero_phase,
                 int n_aero_phase_int_param, int n_aero_phase_float_param,
                 int n_aero_rep, int n_aero_rep_int_param,
                 int n_aero_rep_float_param, int n_aero_rep_env_param,
                 int n_sub_model, int n_sub_model_int_param,
                 int n_sub_model_float_param, int n_sub_model_env_param);
void solver_set_spec_name(void *solver_data, char *spec_name,int size_spec_name, int i);
void solver_initialize(void *solver_data, double *abs_tol, double rel_tol,
                       int max_steps, int max_conv_fails);
#ifdef PMC_DEBUG
int solver_set_debug_out(void *solver_data, bool do_output);
int solver_set_eval_jac(void *solver_data, bool eval_Jac);
#endif
int solver_run(void *solver_data, double *state, double *env, double t_initial,
               double t_final, int n_cells);
void solver_get_statistics(void *solver_data, int *solver_flag, int *num_steps,
                           int *RHS_evals, int *LS_setups,
                           int *error_test_fails, int *NLS_iters,
                           int *NLS_convergence_fails, int *DLS_Jac_evals,
                           int *DLS_RHS_evals, double *last_time_step__s,
                           double *next_time_step__s, int *Jac_eval_fails,
                           int *RHS_evals_total, int *Jac_evals_total,
                           double *RHS_time__s, double *Jac_time__s,
                           double *timeCVode, double *timeCVodeTotal,
                           double *max_loss_precision);
void export_camp_input(void *solver_data, double *init_state, char *path);
void solver_free(void *solver_data);
void model_free(ModelData model_data);

#ifdef PMC_USE_SUNDIALS
/* Functions called by the solver */
int f(realtype t, N_Vector y, N_Vector deriv, void *model_data);
int Jac(realtype t, N_Vector y, N_Vector deriv, SUNMatrix J, void *model_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int guess_helper(const realtype t_n, const realtype h_n, N_Vector y_n,
                 N_Vector y_n1, N_Vector hf, void *solver_data, N_Vector tmp1,
                 N_Vector corr);
void error_handler(int error_code, const char *module, const char *function,
                   char *msg, void *sd);

/* SUNDIALS support functions */
int camp_solver_update_model_state(N_Vector solver_state, ModelData *model_data,
                                   realtype threshhold,
                                   realtype replacement_value);
SUNMatrix get_jac_init(SolverData *solver_data);
bool check_Jac(realtype t, N_Vector y, SUNMatrix J, N_Vector deriv,
               N_Vector tmp, N_Vector tmp1, void *solver_data);
int check_flag(void *flag_value, char *func_name, int opt);
void check_flag_fail(void *flag_value, char *func_name, int opt);
void solver_reset_timers(void *solver_data);
static void solver_print_stats(void *cvode_mem);
static void print_data_sizes(ModelData *md);
static void print_jacobian(SUNMatrix M);
static void print_derivative(N_Vector deriv);
bool is_anything_going_on_here(SolverData *sd, realtype t_initial,
                               realtype t_final);
#ifdef PMC_USE_GSL
double gsl_f(double x, void *param);
typedef struct {
  int ind_var;              // independent variable index
  int dep_var;              // dependent variable index
  realtype t;               // current model time (s)
  N_Vector y;               // dependent variable array
  N_Vector deriv;           // time derivative vector f(t,y)
  SolverData *solver_data;  // solver data
} GSLParam;
#endif
#endif

#endif
