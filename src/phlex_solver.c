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
#include <time.h>
#include <math.h>
#include "phlex_solver.h"
#ifdef PMC_USE_GPU
#include "cuda/phlex_gpu_solver.h"
#endif
#include "aero_rep_solver.h"
#include "rxn_solver.h"
#include "sub_model_solver.h"

// FIXME are these necessary? CVODE already has counters for these calls
// CVODE counters reset each time we call solver_run, so we need to accumulate them
// adding a extra counter to know the total iterations (since the iters can vary between solver_run iters)
unsigned int counterDeriv=0;
unsigned int counterJac=0;
double timeDerivgpu=0;
double timeDeriv=0;
double timeJac=0;

// Define PMC_DEBUG to turn on output of debug info

// Define FAILURE_DETAIL to print out conditions before and after solver failures

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
  printf("\n[DEBUG] line %4d in %-20s(): %-25s %-4.0d t_n = %le h = %le q = %d "
         "hin = %le species %d(zn[0] = %le zn[1] = %le tempv = %le tempv1 = %le "
         "tempv2 = %le acor_init = %le last_yn = %le",
         line, func, message, int_val, cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_q,
         cv_mem->cv_hin, PMC_DEBUG_SPEC_,
         NV_DATA_S(cv_mem->cv_zn[0])[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_zn[1])[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_tempv)[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_tempv1)[PMC_DEBUG_SPEC_],
         //NV_DATA_S(cv_mem->cv_tempv2)[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_acor_init)[PMC_DEBUG_SPEC_],
         NV_DATA_S(cv_mem->cv_last_yn)[PMC_DEBUG_SPEC_]);
  if (do_full) {
    for (int i=0; i<NV_LENGTH_S(cv_mem->cv_y); i++) {
      printf("\n  zn[0][%3d] = % -le zn[1][%3d] = % -le tempv[%3d] = % -le "
             "tempv1[%3d] = % -le tempv2[%3d] = % -le acor_init[%3d] = % -le "
             "last_yn[%3d] = % -le",
         i, NV_DATA_S(cv_mem->cv_zn[0])[i],
         i, NV_DATA_S(cv_mem->cv_zn[1])[i],
         i, NV_DATA_S(cv_mem->cv_tempv)[i],
         i, NV_DATA_S(cv_mem->cv_tempv1)[i],
         //i, NV_DATA_S(cv_mem->cv_tempv2)[i],
         i, NV_DATA_S(cv_mem->cv_acor_init)[i],
         i, NV_DATA_S(cv_mem->cv_last_yn)[i]);
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
// State advancement factor for Jacobian element evaluation
#define JAC_CHECK_ADV 1.0E-8
// Tolerance for Jacobian element evaluation
#define JAC_CHECK_TOL 1.0E-6
// Set MAX_TIMESTEP_WARNINGS to a negative number to prevent output
#define MAX_TIMESTEP_WARNINGS -1
// Maximum number of steps in discreet addition guess helper
#define GUESS_MAX_ITER 50
// Guess advance scaling factor
#define GUESS_ADV_SCALE 0.618

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
 * \param n_state_var Number of variables on the state array per grid cell
 * \param n_cells Number of grid cells to solve simultaneously
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
void * solver_new(int n_state_var, int n_cells, int *var_type, int n_rxn,
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

  // Save the number of state variables per grid cell
  sd->model_data.n_state_var = n_state_var;

  // Set number of cells to compute simultaneously
  sd->model_data.n_cells = n_cells;

  // Add the variable types to the solver data
  sd->model_data.var_type = (int*) malloc(n_state_var * sizeof(int));
  if (sd->model_data.var_type==NULL) {
    printf("\n\nERROR allocating space for variable types\n\n");
    exit(1);
  }
  for (int i=0; i<n_state_var; i++)
    sd->model_data.var_type[i] = var_type[i];

  // Get the number of solver variables per grid cell
  int n_dep_var = 0;
  for (int i=0; i<n_state_var; i++)
    if (var_type[i]==CHEM_SPEC_VARIABLE) n_dep_var++;

  // Save the number of solver variables per grid cell
  sd->model_data.n_dep_var = n_dep_var;

#ifdef PMC_USE_SUNDIALS
  // Set up the solver variable array and helper derivative array
  sd->y     = N_VNew_Serial(n_dep_var*n_cells);
  sd->deriv = N_VNew_Serial(n_dep_var*n_cells);
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
  //TODO: Maybe move this rxn_data away to model_data, avoiding hiding this value
  int *ptr = sd->model_data.rxn_data;
  ptr[0] = n_rxn;
  sd->model_data.nxt_rxn = (void*) &(ptr[1]);

  // If there are no reactions, flag the solver not to run
  sd->no_solve = (n_rxn==0);

  // TODO Move this to rxn_data
  sd->model_data.rate_constants = (void*) malloc(
          (n_rxn*n_cells)*sizeof(double));

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

  //GPU
#ifdef PMC_USE_GPU
  solver_new_gpu_cu(n_dep_var,
  n_state_var, n_rxn,
  n_rxn_int_param, n_rxn_float_param, n_cells);
#endif

#ifdef PMC_DEBUG
  print_data_sizes(sd->model_data);
#endif

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
  int n_cells;                  // number of cells to solve simultaneously
  int *var_type;		// state variable types

  // Get a pointer to the SolverData
  sd = (SolverData*) solver_data;

  // Create a new solver object
  sd->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  check_flag_fail((void *)sd->cvode_mem, "CVodeCreate", 0);

  // Get the number of total and dependent variables on the state array,
  // and the type of each state variable. All values are per-grid-cell.
  n_state_var = sd->model_data.n_state_var;
  n_dep_var = sd->model_data.n_dep_var;
  var_type = sd->model_data.var_type;
  n_cells = sd->model_data.n_cells;

  // Set the solver data
  flag = CVodeSetUserData(sd->cvode_mem, sd);
  check_flag_fail(&flag, "CVodeSetUserData", 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * right-hand side function in y'=f(t,y), the initial time t0, and
   * the initial dependent variable vector y. */
  //TODO: sd->y is deriv for cvode calculations(jac and f)
  flag = CVodeInit(sd->cvode_mem, f, (realtype) 0.0, sd->y);
  check_flag_fail(&flag, "CVodeInit", 1);

  // Set the relative and absolute tolerances
  sd->abs_tol_nv = N_VNew_Serial(n_dep_var*n_cells);
  i_dep_var = 0;
  for (int i_cell=0; i_cell<n_cells; ++i_cell)
    for (int i_var=0; i_var<n_state_var; ++i_var)
      if (var_type[i_var]==CHEM_SPEC_VARIABLE)
            NV_Ith_S(sd->abs_tol_nv, i_dep_var++) = (realtype) abs_tol[i_var];
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

  sd->debug_out = do_output == true ? SUNTRUE : SUNFALSE;
  return PHLEX_SOLVER_SUCCESS;
#endif
}
#endif

/** \brief Set the flag indicating whether to evalute the Jacobian during
 **        solving
 *
 * \param solver_data A pointer to the solver data
 * \param do_output Whether to evaluate the Jacobian during solving
 */
#ifdef PMC_DEBUG
int solver_set_eval_jac(void *solver_data, bool eval_Jac)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;

  sd->eval_Jac = eval_Jac == true ? SUNTRUE : SUNFALSE;
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
		double t_final) {

#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  int flag;

  // Update the dependent variables
  int i_dep_var = 0;
  int n_state_var = sd->model_data.n_state_var;
  int n_cells = sd->model_data.n_cells;
  for (int i_cell=0; i_cell<n_cells; ++i_cell)
    for (int i_spec=0; i_spec<n_state_var; ++i_spec)
      if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE)
        NV_Ith_S(sd->y,i_dep_var++) = (realtype) state[i_spec+i_cell*n_state_var];

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

  // Reset the counter of Jacobian evaluation failures
  sd->Jac_eval_fails = 0;

  // Update data for new environmental state
  // (This is set up to assume the environmental variables do not change during
  //  solving. This can be changed in the future if necessary.)
  aero_rep_update_env_state(&(sd->model_data), env);
  sub_model_update_env_state(&(sd->model_data), env);
  rxn_update_env_state(&(sd->model_data), env);

  PMC_DEBUG_JAC_STRUCT("Begin solving");

  // Reset the flag indicating a current J_guess
  sd->curr_J_guess = false;

  // Set the initial time step
  sd->init_time_step = (t_final - t_initial) * DEFAULT_TIME_STEP;

#ifdef PMC_USE_GPU
  //Set gpu rxn values only 1 time
  //TODO: this should be reordered by setting after initializations and before the run
  solver_set_data_gpu(&(sd->model_data));
#endif

  // Check whether there is anything to solve (filters empty air masses with no
  // emissions)
  if( is_anything_going_on_here( sd, t_initial, t_final ) == false )
    return PHLEX_SOLVER_SUCCESS;
  // CVODE makes 3+ calls to deriv before it decides there is nothing
  // to solve, so this could still be faster - we need to evaluate it
  // in the full MONARCH model

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
  i_dep_var = 0;
  for (int i_cell=0; i_cell<n_cells; ++i_cell) {
    for (int i_spec=0; i_spec<n_state_var; ++i_spec) {
      if (sd->model_data.var_type[i_spec]==CHEM_SPEC_VARIABLE) {
        state[i_spec+i_cell*n_state_var] =
            (double) ( NV_Ith_S(sd->y, i_dep_var) > 0.0 ?
                       NV_Ith_S(sd->y, i_dep_var) : 0.0 );
        ++i_dep_var;
      }
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
 * \param Jac_eval_fails        Number of Jacobian evaluation failures
 */
void solver_get_statistics(void *solver_data, int *num_steps, int *RHS_evals,
                           int *LS_setups, int *error_test_fails,
                           int *NLS_iters, int *NLS_convergence_fails,
                           int *DLS_Jac_evals, int *DLS_RHS_evals,
                           double *last_time_step__s, double *next_time_step__s,
                           int *Jac_eval_fails) {

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
  *Jac_eval_fails = sd->Jac_eval_fails;

#endif
}
//f and jac
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
  int i_dep_var=0;
  int n_state_var = model_data->n_state_var;
  int n_cells = model_data->n_cells;
  for (int i_cell=0; i_cell<n_cells; ++i_cell) {
    for (int i_spec=0; i_spec<n_state_var; ++i_spec) {
      if (model_data->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
        if (NV_DATA_S(solver_state)[i_dep_var] < -SMALL_NUMBER) {
#ifdef FAILURE_DETAIL
          printf("\nFailed model state update: [spec %d] = %le", i_spec,
              NV_DATA_S(solver_state)[i_dep_var]);
#endif
          return PHLEX_SOLVER_FAIL;
        }

        //Assign model state to solver_state
        model_data->state[i_spec+i_cell*n_state_var] =
          NV_DATA_S(solver_state)[i_dep_var] > threshhold ?
          NV_DATA_S(solver_state)[i_dep_var] : replacement_value;
        ++i_dep_var;
      }
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

  // On the first call to f(), the time step hasn't beRen set yet, so use the
  // default value
  time_step = time_step > ZERO ? time_step : sd->init_time_step;

  // Update the state array with the current dependent variable values.
  // Signal a recoverable error (positive return value) for negative
  // concentrations.
  if (phlex_solver_update_model_state(y, md, SMALL, ZERO)
      != PHLEX_SOLVER_SUCCESS) return 1;

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Reset the derivative vector
  N_VConst(ZERO, deriv);

  // Run pre-derivative calculations
  rxn_pre_calc(md, (double) time_step);

#ifndef PMC_USE_GPU

  #ifdef PMC_DEBUG_PRINT

    N_Vector derivtest = N_VNew_Serial(NV_LENGTH_S(deriv));
    N_VConst(ZERO, derivtest);

    clock_t start = clock();

    // Calculate the time derivative f(t,y)
    rxn_calc_deriv(md, deriv, (double) time_step);

    clock_t end = clock();
    timeDeriv+= ((double) (end - start));
    counterDeriv++;

  #else

    // Calculate the time derivative f(t,y)
    rxn_calc_deriv(md, deriv, (double) time_step);

  #endif

#else

  #ifdef PMC_DEBUG_PRINT
  clock_t start2 = clock();

  // Calculate the time derivative f(t,y)
  rxn_calc_deriv_gpu(md, deriv, (double) time_step);

  clock_t end2 = clock();
  timeDerivgpu+= ((double) (end2 - start2));
  counterDeriv++;
  #else

    // Calculate the time derivative f(t,y)
    rxn_calc_deriv(md, deriv, (double) time_step);

  #endif

#endif

#ifndef PMC_DEBUG_PRINT

  if(counterDeriv==10)
    print_derivative(deriv);

#endif

  //Return 0 if success
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

  // !!!! Do not use tmp2 - it is the same as y !!!! //
  // FIXME Find out why cvode is sending tmp2 as y

  // Update the state array with the current dependent variable values
  // Signal a recoverable error (positive return value) for negative
  // concentrations.
  if (phlex_solver_update_model_state(y, md, TINY, TINY)
      != PHLEX_SOLVER_SUCCESS) return 1;

  // Advance the state by a small amount to get more accurate Jac values
  // for species that are currently at zero concentration
  int i_dep_var = 0;
  int n_state_var = md->n_state_var;
  int n_cells = md->n_cells;
  for (int i_cell=0; i_cell<n_cells; ++i_cell) {
    for (int i_spec=0; i_spec<n_state_var; ++i_spec) {
      if (md->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
        md->state[i_spec+i_cell*n_state_var] +=
            NV_DATA_S(deriv)[i_dep_var] * SMALL_NUMBER;
        ++i_dep_var;
      }
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

#ifndef PMC_USE_GPU

  #ifdef PMC_DEBUG_PRINT

    clock_t start = clock();

    // Calculate the Jacobian
    rxn_calc_jac(md, J, time_step);

    clock_t end = clock();
    timeJac+= ((double) (end - start));
    counterJac++;

  #else

    // Calculate the Jacobian
    rxn_calc_jac(md, J, time_step);

  #endif

#else

  #ifdef PMC_DEBUG_PRINT

    clock_t start2 = clock();

    // Calculate the Jacobian
    rxn_calc_jac_gpu(md, J, time_step);

    clock_t end2 = clock();
    timeDerivgpu+= ((double) (end2 - start2));
    counterJac++;

  #else

    // Calculate the Jacobian
    rxn_calc_jac_gpu(md, J, time_step);

  #endif

#endif

#ifdef PMC_DEBUG
  // Evaluate the Jacobian if flagged to do so
  if (sd->eval_Jac==SUNTRUE) {
    if (!check_Jac(t, y, J, deriv, tmp1, tmp3, solver_data)) {
      ++(sd->Jac_eval_fails);
    }
  }

  if(counterJac==1)
    print_jacobian(sd->J);

#endif

  return (0);

}

/** \brief Check a Jacobian for accuracy
 *
 * This function compares Jacobian elements against differences in derivative
 * calculations for small changes to the state array:
 * \f[
 *   J_{ij}(x) = \frac{f_i(x+\sum_j e_j) - f_i(x)}{\epsilon}
 * \f]
 * where \f$\epsilon_j = 10^{-8} \left|x_j\right|\f$
 *
 * \param t Current time [s]
 * \param y Current state array
 * \param J Jacobian matrix to evaluate
 * \param deriv Current derivative \f$f(y)\f$
 * \param tmp Working array the size of \f$y\f$
 * \param tmp1 Working array the size of \f$y\f$
 * \param solver_data Solver data
 * \return True if Jacobian values are accurate, false otherwise
 */
bool check_Jac(realtype t, N_Vector y, SUNMatrix J, N_Vector deriv,
    N_Vector tmp, N_Vector tmp1, void *solver_data) {

  realtype * d_state     = NV_DATA_S(y);
  realtype * d_deriv     = NV_DATA_S(deriv);
  realtype * d_adj_state = NV_DATA_S(tmp);
  realtype * d_adj_deriv = NV_DATA_S(tmp1);
  bool retval = true;

  // Calculate the the derivative for the current state y
  if (f(t, y, deriv, solver_data)!=0) {
    printf("\n Derivative calculation failed.\n");
    return false;
  }

  // Copy the current state into an array that will be used to make small
  // adjustments to the state and recalculate the derivative
  N_VScale(ONE, y, tmp);

  // Loop through the independent variables, adding a small amount to the
  // state value for each, one at a time
  for (int i_ind = 0; i_ind < NV_LENGTH_S(y); ++i_ind) {
    d_adj_state[i_ind] += JAC_CHECK_ADV * fabs(d_state[i_ind]);

    // Recalculate the derivative for the adjusted state
    if (f(t, tmp, tmp1, solver_data)!=0) {
      printf("\n Derivative calculation failed.\n");
      return false;
    }

    // Loop through the Jacobian data and evaluate adjustments to f()
    // based on adjustments to y(i_ind)
    for (int i_elem = SM_INDEXPTRS_S(J)[i_ind];
         i_elem < SM_INDEXPTRS_S(J)[i_ind+1]; ++i_elem ) {
      int i_dep = SM_INDEXVALS_S(J)[i_elem];

      // Calculate f_adj(i_dep) = f(i_dep) + J(i_ind, i_dep) * dy(i_ind)
      realtype jac_adj = (d_deriv[i_dep] +
                          SM_DATA_S(J)[i_elem] * JAC_CHECK_ADV *
                              fabs(d_state[i_ind]));

      // Determine the difference between f_adj(i_dep) from Jacobian-based
      // calculations, and recalculations of f() from the adjusted state y_adj
      realtype jac_diff = d_adj_deriv[i_dep] - jac_adj;
      if (jac_diff != ZERO) {

        // Use the average of the two f_adj(i_dep) to normalize the difference
        realtype f_adj_avg = (fabs(d_adj_deriv[i_dep]) + fabs(jac_adj)) / 2.0;
        jac_diff /= f_adj_avg;

        // Evaluate the difference with a relative and absolute tolerance
        if (fabs(jac_diff) > JAC_CHECK_TOL && fabs(f_adj_avg) > SMALL_NUMBER) {
          printf("\nError in Jacobian[%d][%d]: factor of %le relative to f[%d]"
                 " = %le", i_ind, i_dep, jac_diff, i_dep,
                 (fabs(d_adj_deriv[i_dep]) + fabs(jac_adj)) / 2.0);
          printf("\n  f_J[%d] = %le f_adj[%d] = %le", i_dep, jac_adj, i_dep,
                 d_adj_deriv[i_dep]);
          printf("\n  df_J[%d] = %le df_adj[%d] = %le", i_dep,
                 jac_adj - d_deriv[i_dep], i_dep,
                 d_adj_deriv[i_dep] - d_deriv[i_dep]);
          printf("\n  state[%d] = %le state_adj[%d] = %le", i_ind,
                 d_state[i_ind], i_ind, d_adj_state[i_ind]);
          retval = false;
        }
      }
    }

    // Reset the adjusted state y_adj
    d_adj_state[i_ind] = d_state[i_ind];
  }
  return retval;
}

/** \brief Try to improve guesses of y sent to the linear solver
 *
 * This function checks if there are any negative guessed concentrations,
 * and if there are it calculates a set of initial corrections to the
 * guessed state using the state at time \f$t_{n-1}\f$ and the derivative
 * \f$f_{n-1}\f$ and advancing the state according to:
 * \f[
 *   y_n = y_{n-1} + \sum_{j=1}^m h_j * f_j
 * \f]
 * where \f$h_j\f$ is the largest timestep possible where
 * \f[
 *   y_{j-1} + h_j * f_j > 0
 * \f]
 * and
 * \f[
 *   t_n = t_{n-1} + \sum_{j=1}^m h_j
 * \f]
 *
 * \param t_n Current time [s]
 * \param h_n Current time step size [s] If this is set to zero, the change hf
 *            is assumed to be an adjustment where y_n = y_n1 + hf
 * \param y_n Current guess for \f$y(t_n)\f$
 * \param y_n1 \f$y(t_{n-1})\f$
 * \param hf Current guess for change in \f$y\f$ from \f$t_{n-1}\f$ to
 *            \f$t_n\f$ [input/output]
 * \param solver_data Solver data
 * \param tmp1 Temporary vector for calculations
 * \param corr Vector of calculated adjustments to \f$y(t_n)\f$ [output]
 * \return 1 if corrections were calculated, 0 if not
 */
int guess_helper(const realtype t_n, const realtype h_n, N_Vector y_n,
                 N_Vector y_n1, N_Vector hf, void *solver_data, N_Vector tmp1,
                 N_Vector corr)
{
  SolverData *sd  = (SolverData*) solver_data;
  realtype *ay_n  = NV_DATA_S(y_n);
  realtype *ay_n1 = NV_DATA_S(y_n1);
  realtype *atmp1 = NV_DATA_S(tmp1);
  realtype *acorr = NV_DATA_S(corr);
  realtype *ahf   = NV_DATA_S(hf);
  int n_elem      = NV_LENGTH_S(y_n);

  // Only try improvements when negative concentrations are predicted
  if (N_VMin(y_n) > -SMALL_NUMBER) return 0;

  PMC_DEBUG_PRINT_FULL("Trying to improve guess");

  // Copy \f$y(t_{n-1})\f$ to working array
  N_VScale(ONE, y_n1, tmp1);

  // Get  \f$f(t_{n-1})\f$
  if (h_n > ZERO) {
    N_VScale(ONE/h_n, hf, corr);
  } else {
    N_VScale(ONE, hf, corr);
  }
  PMC_DEBUG_PRINT("Got f0");

  // Advance state interatively
  realtype t_0 = h_n > ZERO ? t_n - h_n : t_n - ONE;
  realtype t_j = ZERO;
  int iter = 0;
  for (; iter < GUESS_MAX_ITER && t_0 + t_j < t_n; iter++) {

    // Calculate \f$h_j\f$
    realtype h_j = t_n - ( t_0 + t_j );
    int i_fast = -1;
    for (int i = 0; i < n_elem; i++) {
      realtype t_star = - atmp1[i] / acorr[i];
      if( ( t_star > ZERO || (t_star == ZERO && acorr[i] < ZERO) )
          && t_star < h_j ) {
        h_j = t_star;
        i_fast = i;
      }
    }

    // Scale h_j unless t_n can be reached
    if( i_fast >= 0 && h_n > ZERO ) h_j *= GUESS_ADV_SCALE;

    // Only make small changes to adjustment vectors used in Newton iteration
    if( h_n == ZERO && ONE - h_j > ((CVodeMem)sd->cvode_mem)->cv_reltol ) return 0;

    // Avoid advancing state past zero
    if( h_n > ZERO ) h_j = nextafter(h_j, ZERO);

    // Advance the state
    N_VLinearSum(ONE, tmp1, h_j, corr, tmp1);

    PMC_DEBUG_PRINT_FULL("Advanced state");

    // If just scaling an adjustment vector, exit the loop
    if( h_n == ZERO ) break;

    h_j = nextafter(h_j, HUGE_VAL);
    h_j = t_0 + t_j + h_j > t_n ? t_n - (t_0 + t_j) : h_j;

    // Advance t_j
    t_j += h_j;

    // Recalculate the time derivative \f$f(t_j)\f$
    if (f(t_0 + t_j, tmp1, corr, solver_data) != 0) {
      PMC_DEBUG_PRINT("Unexpected failure in guess helper!");
      return 0;
    }
    ((CVodeMem)sd->cvode_mem)->cv_nfe++;


#ifdef PMC_DEBUG
    if (iter == GUESS_MAX_ITER-1 && t_0 + t_j < t_n) {
      PMC_DEBUG_PRINT("Max guess iterations reached!");
    }
#endif
  }

  PMC_DEBUG_PRINT_INT("Guessed y_h in steps:", iter);

  // Set the correction vector
  N_VLinearSum(ONE, tmp1, -ONE, y_n, corr);

  // Update the hf vector
  N_VLinearSum(ONE, tmp1, -ONE, y_n1, hf);

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
                                   elements per grid cell*/
  sunindextype n_jac_elem_total;/* number of potentially non-zero Jacobian
                                   elements for whole domain */

  // Number of variables on the state array per grid cell
  // (these are the ids the reactions are initialized with)
  int n_state_var = solver_data->model_data.n_state_var;

  // Number of solver variables per grid cell
  int n_dep_var = solver_data->model_data.n_dep_var;

  // Number of grid cells
  int n_cells = solver_data->model_data.n_cells;

  // Set up the 2D array of flags for a single grid cell
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
  // mechanism reactions for a single grid cell
  rxn_get_used_jac_elem(&(solver_data->model_data), jac_struct);

  // Determine the number of non-zero Jacobian elements per grid cell
  n_jac_elem = 0;
  for (int i=0; i<n_state_var; i++)
    for (int j=0; j<n_state_var; j++)
      if (jac_struct[i][j]==true &&
      solver_data->model_data.var_type[i]==CHEM_SPEC_VARIABLE &&
      solver_data->model_data.var_type[j]==CHEM_SPEC_VARIABLE) n_jac_elem++;

  // Initialize the sparse matrix
  int deriv_length = NV_LENGTH_S(solver_data->y);
  solver_data->model_data.n_jac_elem = (int) n_jac_elem;
  n_jac_elem_total = n_jac_elem * n_cells;
  SUNMatrix M = SUNSparseMatrix(deriv_length, deriv_length,  n_jac_elem_total, CSC_MAT);

  // Set the column and row indices
  int i_col=0, i_elem=0;
  for (int i_cell=0; i_cell<n_cells; ++i_cell){
    for (int i=0; i<n_state_var; i++) {
      if (solver_data->model_data.var_type[i]!=CHEM_SPEC_VARIABLE) continue;
      (SM_INDEXPTRS_S(M))[i_col] = i_elem;
      for (int j=0, i_row=0; j<n_state_var; j++) {
        if (solver_data->model_data.var_type[j]!=CHEM_SPEC_VARIABLE) continue;
        if (jac_struct[j][i]==true) {
      (SM_DATA_S(M))[i_elem] = (realtype) 1.0;
      (SM_INDEXVALS_S(M))[i_elem++] = i_row+n_dep_var*i_cell;
        }
        i_row++;
      }
      i_col++;
    }
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
  for (int
  i_ind=0, i_jac_elem=0; i_ind < n_state_var; i_ind++)
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

/** \brief Print some phlex-chem data sizes
 *
 * \param md Pointer to the model data
 */
static void print_data_sizes(ModelData *md)
{

  int *ptr = md->rxn_data;
  int n_rxn = ptr[0];

  printf("n_rxn: %d " , n_rxn);
  printf("n_state_var: %d" ,md->n_state_var*md->n_cells);
  printf("n_dep_var: %d" ,md->n_dep_var*md->n_cells);

}

/** \brief Print Jacobian matrix in format KLU SPARSE
 *
 * \param M Jacobian matrix
 */
static void print_jacobian(SUNMatrix M)
{

  printf("\n NNZ JAC:\n");
  printf ("%d ", SM_NNZ_S(M));
  printf("INDEXVALS:\n");
  for (int i=0; i<SM_NNZ_S(M); i++) {
    printf ("% -le \n", (SM_DATA_S(M))[i]);
    printf ("%d \n", (SM_INDEXVALS_S(M))[i]);
  }
  printf("PTRS:\n");
  for (int i=0; i<=SM_NP_S(M); i++) {
    printf ("%d \n", (SM_INDEXPTRS_S(M))[i]);
  }

}

/** \brief Print derivative array
 *
 * \param deriv Derivative array
 */
static void print_derivative(N_Vector deriv)
{

  //printf(" deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i=0; i<NV_LENGTH_S(deriv); i++) {//NV_LENGTH_S(deriv)
    printf(" deriv: % -le", NV_DATA_S(deriv)[i]);
    printf(" index: %d \n", i);

  }

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

#ifndef PMC_DEBUG_PRINT
  printf ("Total Time Derivgpu= %f",timeDerivgpu / CLOCKS_PER_SEC);
  printf (", Total Time Deriv= %f",timeDeriv / CLOCKS_PER_SEC);
  printf (", Total Time Jac= %f\n",timeJac / CLOCKS_PER_SEC);
  printf ("counterDeriv: %d ", counterDeriv);
  printf ("counterJac: %d ", counterJac);
#endif

#ifdef PMC_USE_GPU
  free_gpu_cu();
#endif

#ifdef PMC_USE_SUNDIALS
  // Destroy the initialized Jacbobian matrix
  SUNMatDestroy(model_data.J_init);
#endif
  free(model_data.var_type);
  free(model_data.rxn_data);
  free(model_data.aero_phase_data);
  free(model_data.aero_rep_data);
  free(model_data.sub_model_data);
  free(model_data.rate_constants);

}

/** \brief Free update data
 *
 * \param update_data Object to free
 */
void solver_free_update_data(void *update_data)
{
  free(update_data);
}

