/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * This is the c ODE solver for the chemistry module
 * It is currently set up to use the SUNDIALS BDF method, Newton
 * iteration with the KLU sparse linear solver.
 *
 * It uses a scalar relative tolerance and a vector absolute tolerance.
 *
*/
/** \file
 * \brief Interface to c solvers for chemistry
*/
#include <stdio.h>
#include <stdlib.h>
#include "phlex_solver.h"
#include "aero_rep_solver.h"
#include "rxn_solver.h"
#include "sub_model_solver.h"

#ifdef PMC_DEBUG
#define PMC_DEBUG_SPEC_ 0
#define PMC_DEBUG_PRINT(x) pmc_debug_print(sd->cvode_mem, x, false, 0, __LINE__, __func__)
#define PMC_DEBUG_PRINT_INT(x,y) pmc_debug_print(sd->cvode_mem, x, false, y, __LINE__, __func__)
#define PMC_DEBUG_PRINT_FULL(x) pmc_debug_print(sd->cvode_mem, x, true, 0, __LINE__, __func__)
#define PMC_DEBUG_JAC_STRUCT(x) pmc_debug_print_jac_struct((void*)sd, x)
void pmc_debug_print(void *cvode_mem, const char *message, bool do_full,
    const int int_val, const int line, const char *func)
{
  CVodeMem cv_mem = (CVodeMem) cvode_mem;
  if( !(cv_mem->cv_debug_out) ) return;
  printf("\n[DEBUG] line %4d in %-20s(): %-25s %-4.0d t_n = %le h = %le q = %d"
         "hin = %le species %d(zn[0] = %le zn[1] = %le tempv = %le tempv2 = %le)",
         line, func, message, int_val, cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_q,
         cv_mem->cv_hin, PMC_DEBUG_SPEC_,
         NV_DATA_S(cv_mem->cv_zn[0])[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_zn[1])[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_tempv)[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_tempv1)[PMC_DEBUG_SPEC_]);
  if (do_full) {
    for (int i=0; i<NV_LENGTH_S(cv_mem->cv_y); i++) {
      printf("\n  zn[0][%3d] = % -le zn[1][%3d] = % -le tempv[%3d] = % -le tempv2[%3d] = % -le",
         i, NV_DATA_S(cv_mem->cv_zn[0])[i],
         i, NV_DATA_S(cv_mem->cv_zn[1])[i],
         i, NV_DATA_S(cv_mem->cv_tempv)[i],
         i, NV_DATA_S(cv_mem->cv_tempv1)[i]);
    }
  }
}
void pmc_debug_print_jac_struct(void *solver_data, const char *message)
{
  SolverData *sd = (SolverData*) solver_data;

  if( !(sd->debug_out) ) return;
  int n_state_var = NV_LENGTH_S(sd->deriv);
  int i_elem = 0;
  int next_col = 0;
  printf("\n\n   Jacobian structure - %s\n     ", message);
  for (int i_dep=0; i_dep < n_state_var; i_dep++)
    printf("[%3d]", i_dep);
  for (int i_ind=0; i_ind < n_state_var; i_ind++) {
    printf("\n[%3d]", i_ind);
    next_col = SM_INDEXPTRS_S(sd->model_data.J_init)[i_ind+1];
    for (int i_dep=0; i_dep < n_state_var; i_dep++) {
      if (i_dep == SM_INDEXVALS_S(sd->model_data.J_init)[i_elem] &&
          i_elem < next_col) {
        printf(" %3d ", i_elem++);
      } else {
        printf("  -  ");
      }
    }
  }
}
#else
#define PMC_DEBUG_PRINT(x)
#define PMC_DEBUG_PRINT_INT(x,y)
#define PMC_DEBUG_PRINT_FULL(x)
#define PMC_DEBUG_JAC_STRUCT(x)
#endif

// Default solver initial time step relative to total integration time
#define DEFAULT_TIME_STEP 1.0
// Small number used in filtering
#define SMALL_NUMBER 1.1E-30
// Set MAX_TIMESTEP_WARNINGS to a negative number to prevent output
#define MAX_TIMESTEP_WARNINGS -1
// Relative time step size threshhold for trying to improve guess of yn
#define GUESS_THRESHHOLD 1.0

// Status codes for calls to phlex_solver functions
#define PHLEX_SOLVER_SUCCESS 0
#define PHLEX_SOLVER_FAIL 1

// State variable types (Must match parameters defined in pmc_chem_spec_data module)
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

/** \brief Get a new solver object
 *
 * Return a pointer to a new SolverData object
 *
 * \param n_state_var Number of variables on the state array
 * \param var_type Pointer to array of state variable types (solver, constant,
 *                 PSSA)
 * \param n_rxn Number of reactions to include
 * \param n_rxn_int_param Total number of integer reaction parameters
 * \param n_rxn_float_param Total number of floating-point reaction parameters
 * \param n_aero_phase Number of aerosol phases
 * \param n_aero_phase_int_param Total number of integer aerosol phase
 *                               parameters
 * \param n_aero_phase_float_param Total number of floating-point aerosol phase
 *                                 parameters
 * \param n_aero_rep Number of aerosol representations
 * \param n_aero_rep_int_param Total number of integer aerosol representation
 *                             parameters
 * \param n_aero_rep_float_param Total number of floating-point aerosol
 *                               representation parameters
 * \param n_sub_model Number of sub models
 * \param n_sub_model_int_param Total number of integer sub model parameters
 * \param n_sub_model_float_param Total number of floating-point sub model
 *                                parameters
 * \return Pointer to the new SolverData object
 */
void * solver_new(int n_state_var, int *var_type, int n_rxn,
          int n_rxn_int_param, int n_rxn_float_param, int n_aero_phase,
          int n_aero_phase_int_param, int n_aero_phase_float_param,
          int n_aero_rep, int n_aero_rep_int_param, int n_aero_rep_float_param,
          int n_sub_model, int n_sub_model_int_param, int n_sub_model_float_param)
{
  // Create the SolverData object
  SolverData *sd = (SolverData*) malloc(sizeof(SolverData));
  if (sd==NULL) {
    printf("\n\nERROR allocating space for SolverData\n\n");
    exit(1);
  }

#ifdef PMC_DEBUG
  // Default to no debugging output
  sd->debug_out = SUNFALSE;
#endif

  // Save the number of state variables
  sd->model_data.n_state_var = n_state_var;

  // Add the variable types to the solver data
  sd->model_data.var_type = (int*) malloc(n_state_var * sizeof(int));
  if (sd->model_data.var_type==NULL) {
    printf("\n\nERROR allocating space for variable types\n\n");
    exit(1);
  }
  for (int i=0; i<n_state_var; i++)
    sd->model_data.var_type[i] = var_type[i];

  // Create arrays for adjustments to the state array from fast rxns applied
  // during calculations of derivatives and Jacobian
  sd->model_data.state_adj     = (double*) malloc(n_state_var * sizeof(double));
  sd->model_data.rel_flux      = (double*) malloc(n_state_var * sizeof(double));
  if (sd->model_data.state_adj == NULL ||
      sd->model_data.rel_flux  == NULL) {
    printf("\n\nERROR allocating space for state adjustment arrays\n\n");
    exit(1);
  }

  // Get the number of solver variables
  int n_dep_var = 0;
  for (int i=0; i<n_state_var; i++)
    if (var_type[i]==CHEM_SPEC_VARIABLE) n_dep_var++;

#ifdef PMC_USE_SUNDIALS
  // Set up the solver variable array and helper derivative array
  sd->y     = N_VNew_Serial(n_dep_var);
  sd->deriv = N_VNew_Serial(n_dep_var);
#endif

  // Allocate space for the reaction data and set the number
  // of reactions (including one int for the number of reactions
  // and one int per reaction to store the reaction type)
  sd->model_data.rxn_data = (void*) malloc(
		  (n_rxn_int_param + 1 + n_rxn) * sizeof(int)
		  + n_rxn_float_param * sizeof(double));
  if (sd->model_data.rxn_data==NULL) {
    printf("\n\nERROR allocating space for reaction data\n\n");
    exit(1);
  }
  int *ptr = sd->model_data.rxn_data;
  ptr[0] = n_rxn;
  sd->model_data.nxt_rxn = (void*) &(ptr[1]);

  // If there are no reactions, flag the solver not to run
  sd->no_solve = (n_rxn==0);

  // Allocate space for the aerosol phase data and st the number
  // of aerosol phases (including one int for the number of
  // phases)
  sd->model_data.aero_phase_data = (void*) malloc(
                  (n_aero_phase_int_param + 1) * sizeof(int)
                  + n_aero_phase_float_param * sizeof(double));
  if (sd->model_data.aero_phase_data==NULL) {
    printf("\n\nERROR allocating space for aerosol phase data\n\n");
    exit(1);
  }
  ptr = sd->model_data.aero_phase_data;
  ptr[0] = n_aero_phase;
  sd->model_data.nxt_aero_phase = (void*) &(ptr[1]);

  // Allocate space for the aerosol representation data and set
  // the number of aerosol representations (including one int
  // for the number of aerosol representations and one int per
  // aerosol representation to store the aerosol representation
  // type)
  sd->model_data.aero_rep_data = (void*) malloc(
		  (n_aero_rep_int_param + 1 + n_aero_rep) * sizeof(int)
		  + n_aero_rep_float_param * sizeof(double));
  if (sd->model_data.aero_rep_data==NULL) {
    printf("\n\nERROR allocating space for aerosol representation data\n\n");
    exit(1);
  }
  ptr = sd->model_data.aero_rep_data;
  ptr[0] = n_aero_rep;
  sd->model_data.nxt_aero_rep = (void*) &(ptr[1]);

  // Allocate space for the sub model data and set the number of sub models
  // (including one int for the number of sub models and one int per sub
  // model to store the sub model type)
  sd->model_data.sub_model_data = (void*) malloc(
                  (n_sub_model_int_param + 1 + n_sub_model) * sizeof(int)
                  + n_sub_model_float_param * sizeof(double));
  if (sd->model_data.sub_model_data==NULL) {
    printf("\n\nERROR allocating space for sub model data\n\n");
    exit(1);
  }
  ptr = sd->model_data.sub_model_data;
  ptr[0] = n_sub_model;
  sd->model_data.nxt_sub_model = (void*) &(ptr[1]);

  // Return a pointer to the new SolverData object
  return (void*) sd;
}

/** \brief Solver initialization
 *
 * Allocate and initialize solver objects
 *
 * \param solver_data Pointer to a SolverData object
 * \param abs_tol Pointer to array of absolute tolerances
 * \param rel_tol Relative integration tolerance
 * \param max_steps Maximum number of internal integration steps
 * \param max_conv_fails Maximum number of convergence failures
 * \return Pointer to an initialized SolverData object
 */
void solver_initialize(void *solver_data, double *abs_tol, double rel_tol,
	       int max_steps, int max_conv_fails)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd;		// SolverData object
  int flag;			// return code from SUNDIALS functions
  int n_dep_var;		// number of dependent variables
  int i_dep_var; 		// index of dependent variables in loops
  int n_state_var; 		// number of variables on the state array
  int *var_type;		// state variable types

  // Get a pointer to the SolverData
  sd = (SolverData*) solver_data;

  // Create a new solver object
  sd->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  check_flag_fail((void *)sd->cvode_mem, "CVodeCreate", 0);

  // Get the number of total and dependent variables on the state array,
  // and the type of each state variable
  n_state_var = sd->model_data.n_state_var;
  n_dep_var = NV_LENGTH_S(sd->y);
  var_type = sd->model_data.var_type;

  // Set the solver data
  flag = CVodeSetUserData(sd->cvode_mem, sd);
  check_flag_fail(&flag, "CVodeSetUserData", 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * right-hand side function in y'=f(t,y), the initial time t0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(sd->cvode_mem, f, (realtype) 0.0, sd->y);
  check_flag_fail(&flag, "CVodeInit", 1);

  // Set the relative and absolute tolerances
  sd->abs_tol_nv = N_VNew_Serial(n_dep_var);
  i_dep_var = 0;
  for (int i=0; i<n_state_var; i++)
    if (var_type[i]==CHEM_SPEC_VARIABLE)
            NV_Ith_S(sd->abs_tol_nv, i_dep_var++) = (realtype) abs_tol[i];
  flag = CVodeSVtolerances(sd->cvode_mem, (realtype) rel_tol, sd->abs_tol_nv);
  check_flag_fail(&flag, "CVodeSVtolerances", 1);

  // Add a pointer in the model data to the absolute tolerances for use during
  // solving. TODO find a better way to do this
  sd->model_data.abs_tol = abs_tol;

  // Set the maximum number of iterations
  flag = CVodeSetMaxNumSteps(sd->cvode_mem, max_steps);
  check_flag_fail(&flag, "CVodeSetMaxNumSteps", 1);

  // Set the maximum number of convergence failures
  flag = CVodeSetMaxConvFails(sd->cvode_mem, max_conv_fails);
  check_flag_fail(&flag, "CVodeSetMaxConvFails", 1);

  // Set the maximum number of error test failures (TODO make separate input?)
  flag = CVodeSetMaxErrTestFails(sd->cvode_mem, max_conv_fails);
  check_flag_fail(&flag, "CVodeSetMaxErrTestFails", 1);

  // Set the maximum number of warnings about a too-small time step
  flag = CVodeSetMaxHnilWarns(sd->cvode_mem, MAX_TIMESTEP_WARNINGS);
  check_flag_fail(&flag, "CVodeSetMaxHnilWarns", 1);

  // Get the structure of the Jacobian matrix
  sd->J = get_jac_init(sd);
  sd->model_data.J_init = SUNMatClone(sd->J);
  SUNMatCopy(sd->J, sd->model_data.J_init);

  // Create a Jacobian matrix for correcting negative predicted concentrations
  // during solving
  sd->J_guess = SUNMatClone(sd->J);
  SUNMatCopy(sd->J, sd->J_guess);

  // Create a KLU SUNLinearSolver
  sd->ls = SUNKLU(sd->y, sd->J);
  check_flag_fail((void*) sd->ls, "SUNKLU", 0);

  // Attach the linear solver and Jacobian to the CVodeMem object
  flag = CVDlsSetLinearSolver(sd->cvode_mem, sd->ls, sd->J);
  check_flag_fail(&flag, "CVDlsSetLinearSolver", 1);

  // Set the Jacobian function to Jac
  flag = CVDlsSetJacFn(sd->cvode_mem, Jac);
  check_flag_fail(&flag, "CVDlsSetJacFn", 1);

  // Set a function to improve guesses for y sent to the linear solver
  flag = CVodeSetDlsGuessHelper(sd->cvode_mem, guess_helper);
  check_flag_fail(&flag, "CVodeSetDlsGuessHelper", 1);

#ifndef FAILURE_DETAIL
  // Set a custom error handling function
  flag = CVodeSetErrHandlerFn(sd->cvode_mem, error_handler, (void*) sd );
  check_flag_fail(&flag, "CVodeSetErrHandlerFn", 0);
#endif

#endif
}

/** \brief Set the flag indicating whether to output debugging information
 *
 * \param solver_data A pointer to the solver data
 * \param do_output Whether to output debugging information during solving
 */
#ifdef PMC_DEBUG
int solver_set_debug_out(void *solver_data, bool do_output)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;

  if (do_output == true) {
    sd->debug_out = SUNTRUE;
  } else {
    sd->debug_out = SUNFALSE;
  }
  return PHLEX_SOLVER_SUCCESS;
#endif
}
#endif

/** \brief Solve for a given timestep
 *
 * \param solver_data A pointer to the initialized solver data
 * \param state A pointer to the state array
 * \param env A pointer to the array of environmental conditions
 * \param t_initial Initial time (s)
 * \param t_final (s)
 * \return Flag indicating PHLEX_SOLVER_SUCCESS or PHLEX_SOLVER_FAIL
 */
int solver_run(void *solver_data, double *state, double *env, double t_initial,
		double t_final)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;
  int flag;

  // Update the dependent variables
  for (int i_spec=0, i_dep_var=0; i_spec<sd->model_data.n_state_var; i_spec++)
    if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE)
      NV_Ith_S(sd->y,i_dep_var++) = (realtype) state[i_spec];

  // Update model data pointers
  sd->model_data.state = state;
  sd->model_data.env = env;

#ifdef PMC_DEBUG
  // Update the debug output flag in CVODES and the linear solver
  flag = CVodeSetDebugOut(sd->cvode_mem, sd->debug_out);
  check_flag_fail(&flag, "CVodeSetDebugOut", 1);
  flag = SUNKLUSetDebugOut(sd->ls, sd->debug_out);
  check_flag_fail(&flag, "SUNKLUSetDebugOut", 1);
#endif

  // Update data for new environmental state
  // (This is set up to assume the environmental variables do not change during
  //  solving. This can be changed in the future if necessary.)
  aero_rep_update_env_state(&(sd->model_data), env);
  sub_model_update_env_state(&(sd->model_data), env);
  rxn_update_env_state(&(sd->model_data), env);

  // Reset the state adjustment arrays
  sd->model_data.use_adj = true;
  rxn_reset_state_adjustments(&(sd->model_data));

  PMC_DEBUG_JAC_STRUCT("Begin solving");

  // Reset the flag indicating a current J_guess
  sd->curr_J_guess = false;

  // Set the initial time step
  sd->init_time_step = (t_final - t_initial) * DEFAULT_TIME_STEP;

  // Check whether there is anything to solve (filters empty air masses with no
  // emissions)
  if( is_anything_going_on_here( sd, t_initial, t_final ) == false )
    return PHLEX_SOLVER_SUCCESS;

  // Reinitialize the solver
  flag = CVodeReInit(sd->cvode_mem, t_initial, sd->y);
  check_flag_fail(&flag, "CVodeReInit", 1);

  // Reinitialize the linear solver
  flag = SUNKLUReInit(sd->ls, sd->J, SM_NNZ_S(sd->J), SUNKLU_REINIT_PARTIAL);
  check_flag_fail(&flag, "SUNKLUReInit", 1);

  // Set the inital time step
  flag = CVodeSetInitStep(sd->cvode_mem, sd->init_time_step);
  check_flag_fail(&flag, "CVodeSetInitStep", 1);

  // Run the solver
  realtype t_rt = (realtype) t_initial;
  if (!sd->no_solve) {
    flag = CVode(sd->cvode_mem, (realtype) t_final, sd->y, &t_rt, CV_NORMAL);
#ifndef FAILURE_DETAIL
    if (flag < 0 ) {
#else
    if (check_flag(&flag, "CVode", 1)==PHLEX_SOLVER_FAIL) {
      N_Vector deriv = N_VClone(sd->y);
      flag = f(t_initial, sd->y, deriv, sd);
      if (flag!=0) printf("\nCall to f() at failed state failed with flag %d\n", flag);
      printf("temp = %le pressure = %le\n", env[0], env[1]);
      for (int i_spec=0, i_dep_var=0; i_spec<sd->model_data.n_state_var; i_spec++)
        if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE) {
          printf("spec %d = %le deriv = %le\n", i_spec, NV_Ith_S(sd->y,i_dep_var), NV_Ith_S(deriv, i_dep_var));
          i_dep_var++;
        } else {
          printf("spec %d = %le\n", i_spec, state[i_spec]);
        }
      solver_print_stats(sd->cvode_mem);
#endif
      return PHLEX_SOLVER_FAIL;
    }
  }

  // Update the species concentrations on the state array
  for (int i_spec=0, i_dep_var=0; i_spec<sd->model_data.n_state_var; i_spec++) {
    if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      state[i_spec] = (double) ( NV_Ith_S(sd->y, i_dep_var) > 0.0 ?
                                 NV_Ith_S(sd->y, i_dep_var) : 0.0 );
      i_dep_var++;
    }
  }

  // Re-run the pre-derivative calculations to update equilibrium species
  // and apply adjustments to final state
  sub_model_calculate(&(sd->model_data));
  rxn_pre_calc(&(sd->model_data), 0.0);

  return PHLEX_SOLVER_SUCCESS;
#else
  return PHLEX_SOLVER_FAIL;
#endif
}

/** \brief Get solver statistics after an integration attempt
 *
 * \param solver_data           Pointer to the solver data
 * \param num_steps             Pointer to set to the number of integration
 *                              steps
 * \param RHS_evals             Pointer to set to the number of right-hand side
 *                              evaluations
 * \param LS_setups             Pointer to set to the number of linear solver
 *                              setups
 * \param error_test_fails      Pointer to set to the number of error test
 *                              failures
 * \param NLS_iters             Pointer to set to the non-linear solver
 *                              iterations
 * \param NLS_convergence_fails Pointer to set to the non-linear solver
 *                              convergence failures
 * \param DLS_Jac_evals         Pointer to set to the direct linear solver
 *                              Jacobian evaluations
 * \param DLS_RHS_evals         Pointer to set to the direct linear solver
 *                              right-hand side evaluations
 * \param last_time_step__s     Pointer to set to the last time step size [s]
 * \param next_time_step__s     Pointer to set to the next time step size [s]
 */
void solver_get_statistics( void *solver_data, int *num_steps, int *RHS_evals,
                    int *LS_setups, int *error_test_fails,
                    int *NLS_iters, int *NLS_convergence_fails,
                    int *DLS_Jac_evals, int *DLS_RHS_evals,
                    double *last_time_step__s, double *next_time_step__s ) {

#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  realtype last_h, curr_h;
  int flag;

  flag = CVodeGetNumSteps(sd->cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)==PHLEX_SOLVER_FAIL)
          return;
  *num_steps = (int) nst;
  flag = CVodeGetNumRhsEvals(sd->cvode_mem, &nfe);
  if (check_flag(&flag, "CVodeGetNumRhsEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  *RHS_evals = (int) nfe;
  flag = CVodeGetNumLinSolvSetups(sd->cvode_mem, &nsetups);
  if (check_flag(&flag, "CVodeGetNumLinSolveSetups", 1)==PHLEX_SOLVER_FAIL)
          return;
  *LS_setups = (int) nsetups;
  flag = CVodeGetNumErrTestFails(sd->cvode_mem, &netf);
  if (check_flag(&flag, "CVodeGetNumErrTestFails", 1)==PHLEX_SOLVER_FAIL)
          return;
  *error_test_fails = (int) netf;
  flag = CVodeGetNumNonlinSolvIters(sd->cvode_mem, &nni);
  if (check_flag(&flag, "CVodeGetNonlinSolvIters", 1)==PHLEX_SOLVER_FAIL)
          return;
  *NLS_iters = (int) nni;
  flag = CVodeGetNumNonlinSolvConvFails(sd->cvode_mem, &ncfn);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1)==PHLEX_SOLVER_FAIL)
          return;
  *NLS_convergence_fails = ncfn;
  flag = CVDlsGetNumJacEvals(sd->cvode_mem, &nje);
  if (check_flag(&flag, "CVDlsGetNumJacEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  *DLS_Jac_evals = (int) nje;
  flag = CVDlsGetNumRhsEvals(sd->cvode_mem, &nfeLS);
  if (check_flag(&flag, "CVDlsGetNumRhsEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  *DLS_RHS_evals = (int) nfeLS;
  flag = CVodeGetLastStep(sd->cvode_mem, &last_h);
  if (check_flag(&flag, "CVodeGetLastStep", 1)==PHLEX_SOLVER_FAIL)
          return;
  *last_time_step__s = (double) last_h;
  flag = CVodeGetCurrentStep(sd->cvode_mem, &curr_h);
  if (check_flag(&flag, "CVodeGetCurrentStep", 1)==PHLEX_SOLVER_FAIL)
          return;
  *next_time_step__s = (double) curr_h;

#endif
}

#ifdef PMC_USE_SUNDIALS
/** \brief Update the model state from the current solver state
 *
 * \param solver_state Solver state vector
 * \param model_data Pointer to the model data (including the state array)
 * \param threshhold A lower limit for model concentrations below which the
 *                   solver value is replaced with a replacement value
 * \param replacement_value Replacement value for low concentrations
 * \return PHLEX_SOLVER_SUCCESS for successful update or
 *         PHLEX_SOLVER_FAIL for negative concentration
 */
int phlex_solver_update_model_state(N_Vector solver_state,
          ModelData *model_data, realtype threshhold,
          realtype replacement_value)
{
  for (int i_spec=0, i_dep_var=0; i_spec<model_data->n_state_var; i_spec++) {
    if (model_data->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      if (NV_DATA_S(solver_state)[i_dep_var] < -SMALL_NUMBER) {
#ifdef FAILURE_DETAIL
        printf("\nFailed model state update: [spec %d] = %le", i_spec,
            NV_DATA_S(solver_state)[i_dep_var]);
#endif
        return PHLEX_SOLVER_FAIL;
      }
      model_data->state[i_spec] =
        NV_DATA_S(solver_state)[i_dep_var] > threshhold ?
        NV_DATA_S(solver_state)[i_dep_var] : replacement_value;
      i_dep_var++;
    }
  }
  return PHLEX_SOLVER_SUCCESS;
}

/** \brief Compute the time derivative f(t,y)
 *
 * \param t Current model time (s)
 * \param y Dependent variable array
 * \param deriv Time derivative vector f(t,y) to calculate
 * \param solver_data Pointer to the solver data
 * \return Status code
 */
int f(realtype t, N_Vector y, N_Vector deriv, void *solver_data)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  realtype time_step;

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);

  // On the first call to f(), the time step hasn't been set yet, so use the
  // default value
  time_step = time_step > ZERO ? time_step : sd->init_time_step;

  // Update the state array with the current dependent variable values.
  // Signal a recoverable error (positive return value) for negative
  // concentrations.
  if (phlex_solver_update_model_state(y, md, SMALL, ZERO)
      != PHLEX_SOLVER_SUCCESS) return 1;

  // Reset the derivative vector
  N_VConst(ZERO, deriv);

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Run pre-derivative calculations
  rxn_pre_calc(md, (double) time_step);

  // Calculate the time derivative f(t,y)
  rxn_calc_deriv(md, deriv, (double) time_step);

  return (0);

}

/** \brief Compute the Jacobian
 *
 * \param t Current model time (s)
 * \param y Dependent variable array
 * \param deriv Time derivative vector f(t,y)
 * \param J Jacobian to calculate
 * \param solver_data Pointer to the solver data
 * \param tmp1 Unused vector
 * \param tmp2 Unused vector
 * \param tmp3 Unused vector
 * \return Status code
 */
int Jac(realtype t, N_Vector y, N_Vector deriv, SUNMatrix J, void *solver_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  realtype time_step;

  // Update the state array with the current dependent variable values
  // Signal a recoverable error (positive return value) for negative
  // concentrations.
  if (phlex_solver_update_model_state(y, md, TINY, TINY)
      != PHLEX_SOLVER_SUCCESS) return 1;

  // Advance the state by a small amount to get more accurate Jac values
  // for species that are currently at zero concentration
  for (int i_spec=0, i_dep_var=0; i_spec<md->n_state_var; i_spec++) {
    if (md->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      md->state[i_spec] += NV_DATA_S(deriv)[i_dep_var] * SMALL_NUMBER;
      i_dep_var++;
    }
  }

  // Reset the Jacobian
  SM_NNZ_S(J) = SM_NNZ_S(md->J_init);
  for (int i=0; i<SM_NNZ_S(J); i++) {
    (SM_DATA_S(J))[i] = (realtype)0.0;
    (SM_INDEXVALS_S(J))[i] = (SM_INDEXVALS_S(md->J_init))[i];
  }
  for (int i=0; i<=SM_NP_S(J); i++) {
    (SM_INDEXPTRS_S(J))[i] = (SM_INDEXPTRS_S(md->J_init))[i];
  }


  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);

  // Run pre-Jacobian calculations
  rxn_pre_calc(md, (double) time_step);

  // Calculate the Jacobian
  rxn_calc_jac(md, J, time_step);

  return (0);

}

/** \brief Try to improve guesses of y sent to the linear solver
 *
 * \param t_n Current time [s]
 * \param h_n Current time step size [s]
 * \param y_n Current guess for y_n
 * \param hf Current guess for change in y from t_n-1 to t_n
 * \param solver_data Solver data
 * \param tmp1 Temporary vector for calculations
 * \param corr Vector of calculated adjustments to y [output]
 * \return 1 if corrections were calculated, 0 if not
 */
int guess_helper(const realtype t_n, const realtype h_n, N_Vector y_n,
                 N_Vector hf, void *solver_data, N_Vector tmp1, N_Vector corr)
{
  SolverData *sd  = (SolverData*) solver_data;
  realtype *ay_n  = NV_DATA_S(y_n);
  realtype *atmp1 = NV_DATA_S(tmp1);
  realtype *acorr = NV_DATA_S(corr);
  int n_elem      = NV_LENGTH_S(y_n);

  // TODO Figure out if/how to make this work with q > 1
  if (((CVodeMem)sd->cvode_mem)->cv_q > 1) return 0;

  // Only perform this kind-of expensive guess improvement if the step size is
  // getting really small
  if (h_n == ZERO || h_n > sd->init_time_step * GUESS_THRESHHOLD) return 0;

  // Only try improvements when negative concentrations are predicted
  if (N_VMin(y_n) > -SMALL_NUMBER) return 0;

  PMC_DEBUG_PRINT_FULL("Trying to improve guess");

  // Determine if J_guess is current
  if (sd->curr_J_guess)
    if ( ((t_n - h_n) - sd->J_guess_t) > SMALL_NUMBER )
      sd->curr_J_guess = false;

  // Calculate the Jacobian
  N_VLinearSum(ONE, y_n, -ONE, hf, tmp1);
  PMC_DEBUG_PRINT_FULL("Got y0");
  N_VScale(ONE/h_n, hf, corr);
  PMC_DEBUG_PRINT_FULL("Got f0");
  if (!(sd->curr_J_guess)) {
    if (Jac(t_n-h_n, tmp1, corr, sd->J_guess, solver_data,
            NULL, NULL, NULL) != 0) return 0;
    sd->curr_J_guess = true;
    sd->J_guess_t    = (t_n - h_n); // J is calculatd for t_(n-1)
  }
  PMC_DEBUG_PRINT_FULL("Calculated Jacobian");

  // Estimate change in y0 needed to keep values in yn positive
  // assuming that f0(y) ~ y * f0' where df0'/dy = 0
  for (int i = 0; i < n_elem; i++) {
    if (ay_n[i] >= ZERO) {
      atmp1[i] = ZERO;
    } else {
      atmp1[i] = atmp1[i] * (ay_n[i] / (atmp1[i] - ay_n[i]));
    }
  }
  PMC_DEBUG_PRINT_FULL("Estimated primary adjustments in y0");

  // Get changes in hf from adjustments to y0
  SUNMatScaleAddI(ONE, sd->J_guess);
  SUNMatMatvec(sd->J_guess, tmp1, corr);
  N_VScale(h_n, corr, corr);
  PMC_DEBUG_PRINT_FULL("Applied Jacobian to adjustments");

  // Recalculate adjustments to y0
  for (int i = 0; i < n_elem; i++)
    if (ay_n[i] < ZERO && acorr[i] > ZERO)
      atmp1[i] *= (-ay_n[i]/acorr[i]);
  SUNMatMatvec(sd->J_guess, tmp1, corr);
  N_VScale(h_n, corr, corr);
  PMC_DEBUG_PRINT_FULL("Applied Jacobian to recalculated adjustments");

  return 1;
}

/** \brief Create a sparse Jacobian matrix based on model data
 *
 * \param solver_data A pointer to the SolverData
 * \return Sparse Jacobian matrix with all possible non-zero elements intialize
 *         to 1.0
 */
SUNMatrix get_jac_init(SolverData *solver_data)
{
  int n_rxn;			/* number of reactions in the mechanism
  				 * (stored in first position in *rxn_data) */
  bool **jac_struct;		/* structure of Jacobian with flags to indicate
				 * elements that could be used. */
  sunindextype n_jac_elem; 	/* number of potentially non-zero Jacobian
                                   elements */

  // Number of variables on the state array (these are the ids the reactions
  // are initialized with)
  int n_state_var = solver_data->model_data.n_state_var;

  // Set up the 2D array of flags
  jac_struct = (bool**) malloc(sizeof(int*) * n_state_var);
  if (jac_struct==NULL) {
    printf("\n\nERROR allocating space for jacobian structure array\n\n");
    exit(1);
  }
  for (int i_spec=0; i_spec < n_state_var; i_spec++) {
    jac_struct[i_spec] = (bool*) malloc(sizeof(int) * n_state_var);
    if (jac_struct[i_spec]==NULL) {
      printf("\n\nERROR allocating space for jacobian structure array "
                "row %d\n\n", i_spec);
      exit(1);
    }
    // Add diagnonal elements by default
    for (int j_spec=0; j_spec < n_state_var; j_spec++)
            jac_struct[i_spec][j_spec] = i_spec==j_spec ? true : false;
  }

  // Fill in the 2D array of flags with Jacobian elements used by the
  // mechanism reactions
  rxn_get_used_jac_elem(&(solver_data->model_data), jac_struct);

  // Determine the number of non-zero Jacobian elements
  n_jac_elem = 0;
  for (int i=0; i<n_state_var; i++)
    for (int j=0; j<n_state_var; j++)
      if (jac_struct[i][j]==true &&
	  solver_data->model_data.var_type[i]==CHEM_SPEC_VARIABLE &&
	  solver_data->model_data.var_type[j]==CHEM_SPEC_VARIABLE) n_jac_elem++;

  // Initialize the sparse matrix
  int n_dep_var = NV_LENGTH_S(solver_data->y);
  SUNMatrix M = SUNSparseMatrix(n_dep_var, n_dep_var, n_jac_elem, CSC_MAT);

  // Set the column and row indices
  int i_col=0, i_elem=0;
  for (int i=0; i<n_state_var; i++) {
    if (solver_data->model_data.var_type[i]!=CHEM_SPEC_VARIABLE) continue;
    (SM_INDEXPTRS_S(M))[i_col] = i_elem;
    for (int j=0, i_row=0; j<n_state_var; j++) {
      if (solver_data->model_data.var_type[j]!=CHEM_SPEC_VARIABLE) continue;
      if (jac_struct[j][i]==true) {
	(SM_DATA_S(M))[i_elem] = (realtype) 1.0;
	(SM_INDEXVALS_S(M))[i_elem++] = i_row;
      }
      i_row++;
    }
    i_col++;
  }
  (SM_INDEXPTRS_S(M))[i_col] = i_elem;

  // Build the set of time derivative ids
  int *deriv_ids = (int*) malloc(sizeof(int) * n_state_var);
  if (deriv_ids==NULL) {
    printf("\n\nERROR allocating space for derivative ids\n\n");
    exit(1);
  }
  for (int i_spec=0, i_dep_var=0; i_spec < n_state_var; i_spec++)
    if (solver_data->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      deriv_ids[i_spec] = i_dep_var++;
    } else {
      deriv_ids[i_spec] = -1;
    }

  // Build the set of Jacobian ids
  int **jac_ids = (int**) jac_struct;
  for (int i_ind=0, i_jac_elem=0; i_ind < n_state_var; i_ind++)
    for (int i_dep=0; i_dep < n_state_var; i_dep++)
      if (solver_data->model_data.var_type[i_dep]==CHEM_SPEC_VARIABLE &&
          solver_data->model_data.var_type[i_ind]==CHEM_SPEC_VARIABLE &&
	  jac_struct[i_dep][i_ind]==true) {
	jac_ids[i_dep][i_ind] = i_jac_elem++;
      } else {
	jac_ids[i_dep][i_ind] = -1;
      }

  // Update the ids in the reaction data
  rxn_update_ids(&(solver_data->model_data), deriv_ids, jac_ids);

  // Free the memory used
  for (int i_spec=0; i_spec<n_state_var; i_spec++) free(jac_struct[i_spec]);
  free(jac_struct);
  free(deriv_ids);

  return M;

}

/** \brief Check the return value of a SUNDIALS function
 *
 * \param flag_value A pointer to check (either for NULL, or as an int pointer
 *                   giving the flag value
 * \param func_name A string giving the function name returning this result code
 * \param opt A flag indicating the type of check to perform (0 for NULL
 *            pointer check; 1 for integer flag check)
 * \return Flag indicating PHLEX_SOLVER_SUCCESS or PHLEX_SOLVER_FAIL
 */
int check_flag(void *flag_value, char *func_name, int opt)
{
  int *err_flag;

  /* Check for a NULL pointer */
  if (opt==0 && flag_value == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
		    func_name);
    return PHLEX_SOLVER_FAIL;
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    err_flag = (int *) flag_value;
    if (*err_flag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
		      func_name, *err_flag);
      return PHLEX_SOLVER_FAIL;
    }
  }
  return PHLEX_SOLVER_SUCCESS;
}

/** \brief Check the return value of a SUNDIALS function and exit on failure
 *
 * \param flag_value A pointer to check (either for NULL, or as an int pointer
 *                   giving the flag value
 * \param func_name A string giving the function name returning this result code
 * \param opt A flag indicating the type of check to perform (0 for NULL
 *            pointer check; 1 for integer flag check)
 */
void check_flag_fail(void *flag_value, char *func_name, int opt)
{
  if (check_flag(flag_value, func_name, opt)==PHLEX_SOLVER_FAIL) {
    exit(EXIT_FAILURE);
  }
}

/** \brief Print solver statistics
 *
 * \param cvode_mem Solver object
 */
static void solver_print_stats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  realtype last_h, curr_h;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if (check_flag(&flag, "CVodeGetNumRhsEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if (check_flag(&flag, "CVodeGetNumLinSolveSetups", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if (check_flag(&flag, "CVodeGetNumErrTestFails", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if (check_flag(&flag, "CVodeGetNonlinSolvIters", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  if (check_flag(&flag, "CVDlsGetNumJacEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  if (check_flag(&flag, "CVDlsGetNumRhsEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  if (check_flag(&flag, "CVodeGetNumGEvals", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetLastStep(cvode_mem, &last_h);
  if (check_flag(&flag, "CVodeGetLastStep", 1)==PHLEX_SOLVER_FAIL)
          return;
  flag = CVodeGetCurrentStep(cvode_mem, &curr_h);
  if (check_flag(&flag, "CVodeGetCurrentStep", 1)==PHLEX_SOLVER_FAIL)
          return;

  printf("\nSUNDIALS Solver Statistics:\n");
  printf("number of steps = %-6ld RHS evals = %-6ld LS setups = %-6ld\n", nst,
            nfe, nsetups);
  printf("error test fails = %-6ld LS iters = %-6ld NLS iters = %-6ld\n", netf,
            nni, ncfn);
  printf("NL conv fails = %-6ld Dls Jac evals = %-6ld Dls RHS evals = %-6ld G evals ="
            " %-6ld\n", ncfn, nje, nfeLS, nge);
  printf("Last time step = %le Next time step = %le\n", last_h, curr_h);
}

#endif

/** \brief Free a SolverData object
 *
 * \param solver_data Pointer to the SolverData object to free
 */
void solver_free(void *solver_data)
{
  SolverData *sd = (SolverData*) solver_data;

#ifdef PMC_USE_SUNDIALS
  // free the SUNDIALS solver
  CVodeFree(&(sd->cvode_mem));

  // free the absolute tolerance vector
  N_VDestroy(sd->abs_tol_nv);

  // free the derivative vector
  N_VDestroy(sd->y);

  // destroy the Jacobian marix
  SUNMatDestroy(sd->J);

  // destroy Jacobian matrix for guessing state
  SUNMatDestroy(sd->J_guess);

  // free the linear solver
  SUNLinSolFree(sd->ls);
#endif

  // Free the allocated ModelData
  model_free(sd->model_data);

  // free the SolverData object
  free(sd);

}

#ifdef PMC_USE_SUNDIALS
/** \brief Determine if there is anything to solve
 *
 * If the solver state concentrations and the derivative vector are very small,
 * there is no point running the solver
 */
bool is_anything_going_on_here(SolverData *sd, realtype t_initial, realtype t_final) {

  if( f(t_initial, sd->y, sd->deriv, sd) ) {
    for (int i_spec=0, i_dep_var=0; i_spec<sd->model_data.n_state_var; i_spec++) {
      if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE) {
        if( NV_Ith_S(sd->y, i_dep_var) >
            NV_Ith_S(sd->abs_tol_nv, i_dep_var) * 1.0e-10 ) return true;
        if( NV_Ith_S(sd->deriv, i_dep_var) * ( t_final - t_initial ) >
            NV_Ith_S(sd->abs_tol_nv, i_dep_var) * 1.0e-10 ) return true;
        i_dep_var++;
      }
    }
    return false;
  }

  return true;
}
#endif

/** \brief Custom error handling function
 *
 * This is used for quiet operation. Solver failures are returned with a flag
 * from the solver_run() function.
 */
void error_handler(int error_code, const char *module,
      const char *function, char *msg, void *sd) {
  // Do nothing
}

/** \brief Free a ModelData object
 *
 * \param model_data Pointer to the ModelData object to free
 */
void model_free(ModelData model_data)
{

#ifdef PMC_USE_SUNDIALS
  // Destroy the initialized Jacbobian matrix
  SUNMatDestroy(model_data.J_init);
#endif
  free(model_data.var_type);
  free(model_data.rxn_data);
  free(model_data.aero_phase_data);
  free(model_data.aero_rep_data);
  free(model_data.sub_model_data);

}

/** \brief Free update data
 *
 * \param update_data Object to free
 */
void solver_free_update_data(void *update_data)
{
  free(update_data);
}

