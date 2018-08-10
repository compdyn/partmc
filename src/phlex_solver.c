/* Copyright (C) 2015-2017 Matthew Dawson
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
#include "phlex_solver.h"

/* Header files with GPU structures and functions */
#ifdef PMC_USE_GPU
#include "phlex_gpu_solver.h"
#endif

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
 * \param n_state_var Number of variables per state on the state array
 * \param n_env_var Number of variables per state on the environmental state
 *                  array
 * \param n_states Number of states to solve simultaneously
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
void * solver_new(int n_state_var, int n_env_var, int n_states, int *var_type,
          int n_rxn, int n_rxn_int_param, int n_rxn_float_param, int n_aero_phase,
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

  // Save the number of state variables
  sd->model_data.n_state_var = n_state_var;

  // Save the number of environmental variables
  sd->model_data.n_env_var = n_env_var;

  // Save the number of states to solve simultaneously
  sd->model_data.n_states = n_states;

  // Add the variable types to the solver data
  sd->model_data.var_type = (int*) malloc(n_state_var * sizeof(int));
  if (sd->model_data.var_type==NULL) {
    printf("\n\nERROR allocating space for variable types\n\n");
    exit(1);
  }
  for (int i=0; i<n_state_var; i++)
    sd->model_data.var_type[i] = var_type[i];

  // Get the number of solver variables
  int n_dep_var = 0;
  for (int i=0; i<n_state_var; i++) 
    if (var_type[i]==CHEM_SPEC_VARIABLE) n_dep_var++;

#ifdef PMC_USE_SUNDIALS
  // Set up the solver variable array
  sd->y = N_VNew_Serial(n_dep_var * n_states);
#endif

#ifndef PMC_USE_DOUBLE_PRECISION
  // Set up a working derivative array
  sd->model_data.deriv = (PMC_SOLVER_C_FLOAT*) malloc(n_dep_var * 
                                        n_states *
                                        sizeof(PMC_SOLVER_C_FLOAT));
  if (sd->model_data.deriv==NULL) {
    printf("\n\nERROR allocating space for working derivative\n\n");
    exit(1);
  }

  // Save the size of the derivative array
  sd->model_data.deriv_size = n_dep_var * n_states;
#endif

  // Allocate space for the reaction data and set the number
  // of reactions (including one int for the number of reactions
  // and one int per reaction to store the reaction type)
  sd->model_data.rxn_data_size = 
	          (n_rxn_int_param + 1 + n_rxn) * sizeof(int) 
		  + n_rxn_float_param * sizeof(PMC_C_FLOAT);
  sd->model_data.rxn_data = (void*) malloc(
                  sd->model_data.rxn_data_size * n_states);
  if (sd->model_data.rxn_data==NULL) {
    printf("\n\nERROR allocating space for reaction data\n\n");
    exit(1);
  }
  int *ptr = sd->model_data.rxn_data;
  for (int i=0; i<n_states; i++) 
    ptr[i*sd->model_data.rxn_data_size/sizeof(int)] = n_rxn;
  sd->model_data.nxt_rxn = (void*) &(ptr[1]);

  // If there are no reactions, flag the solver not to run
  sd->no_solve = (n_rxn==0);

  // Allocate space for the aerosol phase data and st the number
  // of aerosol phases (including one int for the number of
  // phases)
  sd->model_data.aero_phase_data_size = 
                  (n_aero_phase_int_param + 1) * sizeof(int)
                  + n_aero_phase_float_param * sizeof(PMC_C_FLOAT);
  sd->model_data.aero_phase_data = (void*) malloc(
                  sd->model_data.aero_phase_data_size * n_states);
  if (sd->model_data.aero_phase_data==NULL) {
    printf("\n\nERROR allocating space for aerosol phase data\n\n");
    exit(1);
  }
  ptr = sd->model_data.aero_phase_data;
  for (int i=0; i<n_states; i++) 
    ptr[i*sd->model_data.aero_phase_data_size/sizeof(int)] = n_aero_phase;
  sd->model_data.nxt_aero_phase = (void*) &(ptr[1]);

  // Allocate space for the aerosol representation data and set
  // the number of aerosol representations (including one int
  // for the number of aerosol representations and one int per
  // aerosol representation to store the aerosol representation
  // type)
  sd->model_data.aero_rep_data_size = 
		  (n_aero_rep_int_param + 1 + n_aero_rep) * sizeof(int)
		  + n_aero_rep_float_param * sizeof(PMC_C_FLOAT);
  sd->model_data.aero_rep_data = (void*) malloc(
                  sd->model_data.aero_rep_data_size * n_states);
  if (sd->model_data.aero_rep_data==NULL) {
    printf("\n\nERROR allocating space for aerosol representation data\n\n");
    exit(1);
  }
  ptr = sd->model_data.aero_rep_data;
  for (int i=0; i<n_states; i++) 
    ptr[i*sd->model_data.aero_rep_data_size/sizeof(int)] = n_aero_rep;
  sd->model_data.nxt_aero_rep = (void*) &(ptr[1]);

  // Allocate space for the sub model data and set the number of sub models
  // (including one int for the number of sub models and one int per sub
  // model to store the sub model type)
  sd->model_data.sub_model_data_size = 
                  (n_sub_model_int_param + 1 + n_sub_model) * sizeof(int)
                  + n_sub_model_float_param * sizeof(PMC_C_FLOAT);
  sd->model_data.sub_model_data = (void*) malloc(
                  sd->model_data.sub_model_data_size * n_states);
  if (sd->model_data.sub_model_data==NULL) {
    printf("\n\nERROR allocating space for sub model data\n\n");
    exit(1);
  }
  ptr = sd->model_data.sub_model_data;
  for (int i=0; i<n_states; i++) 
    ptr[i*sd->model_data.sub_model_data_size/sizeof(int)] = n_sub_model;
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
void solver_initialize(void *solver_data, PMC_C_FLOAT *abs_tol, PMC_C_FLOAT rel_tol,
	       int max_steps, int max_conv_fails)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd;		// SolverData object
  int flag;			// return code from SUNDIALS functions
  int n_dep_var;		// number of dependent variables
  int i_dep_var; 		// index of dependent variables in loops
  int n_state_var; 		// number of variables on the state array
  int *var_type;		// state variable types
  int n_states;                 // number of states to solve at once

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
  n_states = sd->model_data.n_states;

  // Set the solver data
  flag = CVodeSetUserData(sd->cvode_mem, sd);
  check_flag_fail(&flag, "CVodeSetUserData", 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * right-hand side function in y'=f(t,y), the initial time t0, and
   * the initial dependent variable vector y. */
#ifdef PMC_USE_GPU
  flag = CVodeInit(sd->cvode_mem, phlex_gpu_solver_f, (realtype) 0.0, sd->y);
#else
  flag = CVodeInit(sd->cvode_mem, f, (realtype) 0.0, sd->y);
#endif
  check_flag_fail(&flag, "CVodeInit", 1);

  // Set the relative and absolute tolerances
  i_dep_var = 0;
  sd->abs_tol_nv = N_VNew_Serial(n_dep_var);
  for (int i_state=0; i_state<n_states; i_state++)
    for (int i_var=0; i_var<n_state_var; i_var++)
      if (var_type[i_var]==CHEM_SPEC_VARIABLE)
              NV_Ith_S(sd->abs_tol_nv, i_dep_var++) = (realtype) abs_tol[i_var];
  flag = CVodeSVtolerances(sd->cvode_mem, (realtype) rel_tol, sd->abs_tol_nv);
  check_flag_fail(&flag, "CVodeSVtolerances", 1);

  // Set the maximum number of iterations
  flag = CVodeSetMaxNumSteps(sd->cvode_mem, max_steps);
  check_flag_fail(&flag, "CVodeSetMaxNumSteps", 1);

  // Set the maximum number of convergence failures
  flag = CVodeSetMaxConvFails(sd->cvode_mem, max_conv_fails);
  check_flag_fail(&flag, "CVodeSetMaxConvFails", 1);

  // Get the structure of the Jacobian matrix
  sd->J = get_jac_init(sd);
  sd->model_data.J_init = SUNMatClone(sd->J);
  SUNMatCopy(sd->J, sd->model_data.J_init);

  // Create a KLU SUNLinearSolver
  sd->ls = SUNKLU(sd->y, sd->J);
  check_flag_fail((void*) sd->ls, "SUNKLU", 0);

  // Attach the linear solver and Jacobian to the CVodeMem object
  flag = CVDlsSetLinearSolver(sd->cvode_mem, sd->ls, sd->J);
  check_flag_fail(&flag, "CVDlsSetLinearSolver", 1);

  // Set the Jacobian function to Jac
#ifdef PMC_USE_GPU
  flag = CVDlsSetJacFn(sd->cvode_mem, phlex_gpu_solver_Jac);
#else
  flag = CVDlsSetJacFn(sd->cvode_mem, Jac);
#endif
  check_flag_fail(&flag, "CVDlsSetJacFn", 1);

#ifndef PMC_USE_DOUBLE_PRECISION
  // Set up a working Jacobian data array
  sd->model_data.jac = (PMC_SOLVER_C_FLOAT*) malloc( SM_NNZ_S(sd->J) *
                                          sizeof(PMC_SOLVER_C_FLOAT) );
  if (sd->model_data.jac==NULL) {
    printf("\n\nERROR allocating space for working Jacobian\n\n");
    exit(1);
  }

  // Save the size of the Jacobian data array
  sd->model_data.jac_size = SM_NNZ_S(sd->J);
#endif

#ifdef PMC_USE_GPU
  // Set up the device data for solving with GPUs
  phlex_gpu_solver_new( &(sd->model_data) );
#endif

#endif
}

/** \brief Solve for a given timestep
 *
 * \param solver_data A pointer to the initialized solver data
 * \param state A pointer to the state array
 * \param env A pointer to the array of environmental conditions
 * \param t_initial Initial time (s)
 * \param t_final (s)
 * \return Flag indicating PHLEX_SOLVER_SUCCESS or PHLEX_SOLVER_FAIL
 */
int solver_run(void *solver_data, PMC_C_FLOAT *state, PMC_C_FLOAT *env,
                PMC_C_FLOAT t_initial, PMC_C_FLOAT t_final)
{
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;

  // Update the dependent variables
  for (int i_spec=0, i_dep_var=0, i_state=0; i_state < sd->model_data.n_states;
      i_state++)
    for (int j_spec=0; j_spec<sd->model_data.n_state_var; i_spec++, j_spec++)
      if (sd->model_data.var_type[j_spec]==CHEM_SPEC_VARIABLE) 
        NV_Ith_S(sd->y,i_dep_var++) = (realtype) state[i_spec];

  // Update model data pointers
  sd->model_data.state = state;
  sd->model_data.env = env;

  // Update data for new environmental state
  // (This is set up to assume the environmental variables do not change during
  //  solving. This can be changed in the future if necessary.)
  aero_rep_update_env_state(&(sd->model_data), env);
  sub_model_update_env_state(&(sd->model_data), env);
#ifdef PMC_USE_GPU
  phlex_gpu_solver_update_env_state(&(sd->model_data), env);
#else
  rxn_update_env_state(&(sd->model_data), env);
#endif

  // Reinitialize the solver
  int flag = CVodeReInit(sd->cvode_mem, t_initial, sd->y);
  check_flag_fail(&flag, "CVodeReInit", 1);

  // Run the solver
  realtype t_rt = (realtype) t_initial;
  if (!sd->no_solve) {
    flag = CVode(sd->cvode_mem, (realtype) t_final, sd->y, &t_rt, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1)==PHLEX_SOLVER_FAIL) {
      for( int i_spec=0; i_spec<sd->model_data.n_state_var *
                         sd->model_data.n_states; i_spec++ )
        printf("\nspec %d conc = %le", i_spec, state[i_spec]);
      solver_print_stats(sd->cvode_mem);
      return PHLEX_SOLVER_FAIL;
    }
  }

  // Update the species concentrations on the state array
  for (int i_spec=0, i_dep_var=0, i_state=0; i_state < sd->model_data.n_states;
      i_state++)
    for (int j_spec=0; j_spec<sd->model_data.n_state_var; i_spec++, j_spec++)
      if (sd->model_data.var_type[j_spec]==CHEM_SPEC_VARIABLE) 
        state[i_spec] = (PMC_C_FLOAT) NV_Ith_S(sd->y,i_dep_var++);

  // Re-run the pre-derivative calculations to update equilibrium species
  sub_model_calculate(&(sd->model_data));
  rxn_pre_calc(&(sd->model_data));

  //solver_print_stats(sd->cvode_mem);

  return PHLEX_SOLVER_SUCCESS;
#else
  return PHLEX_SOLVER_FAIL;
#endif
}

#ifdef PMC_USE_SUNDIALS
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
  PMC_SOLVER_C_FLOAT *deriv_data;

  // Update the state array with the current dependent variable values
  // Signal a recoverable error (positive return value) for negative 
  // concentrations.
  for (int i_spec=0, i_dep_var=0, i_state=0; i_state < sd->model_data.n_states;
      i_state++) {
    for (int j_spec=0; j_spec<sd->model_data.n_state_var; i_spec++, j_spec++) {
      if (md->var_type[j_spec]==CHEM_SPEC_VARIABLE) {
        if (NV_DATA_S(y)[i_dep_var] < 0.0) return 1;
        md->state[i_spec] = (PMC_C_FLOAT) (NV_DATA_S(y)[i_dep_var++]);
      }
    }
  }

#ifdef PMC_USE_DOUBLE_PRECISION
  deriv_data = NV_DATA_S(deriv);
#else
  deriv_data = md->deriv;  
#endif

  // Initialize the derivative
  for (int i_spec=0; i_spec<NV_LENGTH_S(deriv); i_spec++)
          deriv_data[i_spec] = (PMC_SOLVER_C_FLOAT) ZERO;

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Run pre-derivative calculations
  rxn_pre_calc(md);

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);
  
  // Calculate the time derivative f(t,y)
  rxn_calc_deriv(md, deriv_data, (PMC_C_FLOAT) time_step);

#ifndef PMC_USE_DOUBLE_PRECISION
  for (int i_spec=0; i_spec<NV_LENGTH_S(deriv); i_spec++)
          NV_DATA_S(deriv)[i_spec] = (realtype) (deriv_data[i_spec]);
#endif

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
  PMC_SOLVER_C_FLOAT *J_data;

  // Update the state array with the current dependent variable values
  for (int i_spec=0, i_dep_var=0, i_state=0; i_state < sd->model_data.n_states;
      i_state++) {
    for (int j_spec=0; j_spec<sd->model_data.n_state_var; i_spec++, j_spec++) {
      if (md->var_type[j_spec]==CHEM_SPEC_VARIABLE) {
        if (NV_DATA_S(y)[i_dep_var] < 0.0) return 1;
        md->state[i_spec] = NV_DATA_S(y)[i_dep_var++];
      }
    }
  }

  // TODO Figure out how to keep the Jacobian from being redimensioned
  // Reset the Jacobian dimensions
  if (SM_NNZ_S(J)<SM_NNZ_S(md->J_init)) {
    SM_INDEXVALS_S(J) = realloc(SM_INDEXVALS_S(J),
              SM_NNZ_S(md->J_init)*sizeof(sunindextype));
    if (SM_INDEXVALS_S(J)==NULL) {
      printf("\n\nERROR allocating space for sparse matrix index values\n\n");
      exit(1);
    }
    SM_DATA_S(J) = realloc(SM_DATA_S(J), SM_NNZ_S(md->J_init)*sizeof(realtype));
    if (SM_DATA_S(J)==NULL) {
      printf("\n\nERROR allocating space for sparse matrix data\n\n");
      exit(1);
    }
  }
  SM_NNZ_S(J) = SM_NNZ_S(md->J_init);
#ifdef PMC_USE_DOUBLE_PRECISION
  J_data = SM_DATA_S(J);
#else
  J_data = md->jac;
#endif
  for (int i=0; i<SM_NNZ_S(J); i++) {
    J_data[i] = (PMC_SOLVER_C_FLOAT) 0.0;
    (SM_INDEXVALS_S(J))[i] = (SM_INDEXVALS_S(md->J_init))[i];
  }
  for (int i=0; i<=SM_NP_S(J); i++) {
    (SM_INDEXPTRS_S(J))[i] = (SM_INDEXPTRS_S(md->J_init))[i];
  } 

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Run pre-Jacobian calculations
  rxn_pre_calc(md);

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);
  
  // Calculate the Jacobian
  rxn_calc_jac(md, J_data, time_step);

#ifndef PMC_USE_DOUBLE_PRECISION
  for (int i=0; i<SM_NNZ_S(J); i++)
    SM_DATA_S(J)[i] = (realtype) (J_data[i]);
#endif

  return (0);

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
  pmc_bool **jac_struct;	/* structure of Jacobian with flags to indicate
				 * elements that could be used. */
  sunindextype n_jac_elem; 	/* number of potentially non-zero Jacobian
                                   elements */

  // Number of variables on the state array (these are the ids the reactions
  // are initialized with)
  int n_state_var = solver_data->model_data.n_state_var;

  // Set up the 2D array of flags
  jac_struct = (pmc_bool**) malloc(sizeof(int*) * n_state_var);
  if (jac_struct==NULL) {
    printf("\n\nERROR allocating space for jacobian structure array\n\n");
    exit(1);
  }
  for (int i_spec=0; i_spec < n_state_var; i_spec++) {
    jac_struct[i_spec] = (pmc_bool*) malloc(sizeof(int) * n_state_var);
    if (jac_struct[i_spec]==NULL) {
      printf("\n\nERROR allocating space for jacobian structure array "
                "row %d\n\n", i_spec);
      exit(1);
    }
    for (int j_spec=0; j_spec < n_state_var; j_spec++)
            jac_struct[i_spec][j_spec] = pmc_false;
  }

  // Fill in the 2D array of flags with Jacobian elements used by the
  // mechanism reactions
  rxn_get_used_jac_elem(&(solver_data->model_data), jac_struct);

  // Determine the number of non-zero Jacobian elements
  n_jac_elem = 0;
  for (int i=0; i<n_state_var; i++)
    for (int j=0; j<n_state_var; j++)
      if (jac_struct[i][j]==pmc_true &&
	  solver_data->model_data.var_type[i]==CHEM_SPEC_VARIABLE &&
	  solver_data->model_data.var_type[j]==CHEM_SPEC_VARIABLE) n_jac_elem++;

  // Scale the number of Jacobian elements by the number of states to solve at once
  n_jac_elem *= solver_data->model_data.n_states;

  // Initialize the sparse matrix
  int n_dep_var = NV_LENGTH_S(solver_data->y);
  SUNMatrix M = SUNSparseMatrix(n_dep_var, n_dep_var, n_jac_elem, CSC_MAT);

  // Set the column and row indices
  int i_col=0, i_elem=0;
  for (int i_state=0; i_state<solver_data->model_data.n_states; i_state++) {
    for (int i=0; i<n_state_var; i++) {
      if (solver_data->model_data.var_type[i]!=CHEM_SPEC_VARIABLE) continue;
      (SM_INDEXPTRS_S(M))[i_col] = i_elem;
      for (int j=0, i_row=i_state*(n_dep_var/solver_data->model_data.n_states); 
          j<n_state_var; j++) {
        if (solver_data->model_data.var_type[j]!=CHEM_SPEC_VARIABLE) continue;
        if (jac_struct[j][i]==pmc_true) {
	  (SM_DATA_S(M))[i_elem] = (realtype) 1.0;
	  (SM_INDEXVALS_S(M))[i_elem++] = i_row;
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
  for (int i_ind=0, i_jac_elem=0; i_ind < n_state_var; i_ind++)
    for (int i_dep=0; i_dep < n_state_var; i_dep++)
      if (solver_data->model_data.var_type[i_dep]==CHEM_SPEC_VARIABLE &&
          solver_data->model_data.var_type[i_ind]==CHEM_SPEC_VARIABLE &&
	  jac_struct[i_dep][i_ind]==pmc_true) {
	jac_ids[i_dep][i_ind] = i_jac_elem++;
      } else {
	jac_ids[i_dep][i_ind] = -1;
      }

  // Update the ids in the data blocks
  n_dep_var /= solver_data->model_data.n_states;
  n_jac_elem /= solver_data->model_data.n_states;
  rxn_update_ids(&(solver_data->model_data), n_dep_var, n_jac_elem, deriv_ids, jac_ids);
  aero_rep_update_ids(&(solver_data->model_data), n_dep_var, n_jac_elem, deriv_ids, jac_ids);
  sub_model_update_ids(&(solver_data->model_data), n_dep_var, n_jac_elem, deriv_ids, jac_ids);

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

  printf("\nSUNDIALS Solver Statistics:\n");
  printf("number of steps = %-6ld RHS evals = %-6ld LS setups = %-6ld\n", nst,
            nfe, nsetups);
  printf("error test fails = %-6ld LS iters = %-6ld NLS iters = %-6ld\n", netf,
            nni, ncfn);
  printf("NL conv fails = %-6ld Jac evals = %-6ld RHS evals = %-6ld G evals ="
            " %-6ld\n", ncfn, nje, nfeLS, nge);
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

  // free the linear solver
  SUNLinSolFree(sd->ls);
#endif

  // Free the allocated ModelData
  model_free(sd->model_data);

  // free the SolverData object
  free(sd);

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
#ifndef PMC_USE_DOUBLE_PRECISION
  free(model_data.deriv);
  free(model_data.jac);
#endif
#ifdef PMC_USE_GPU
  // FIXME For some reason, freeing host pointers causes a segmentation
  // fault when no reactions are solved on the GPUs
  //phlex_gpu_solver_free( model_data.model_dev_data );
#endif
}

/** \brief Free update data
 *
 * \param update_data Object to free
 */
void solver_free_update_data(void *update_data)
{
  free(update_data);
}

