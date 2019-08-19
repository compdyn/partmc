
#ifndef PHLEX_GPU_SOLVER_H_
#define PHLEX_GPU_SOLVER_H_
//#ifndef RXN_SOLVER_H
//#define RXN_SOLVER_H
//#include "phlex_common.h"
#include <cuda.h>
//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "../rxns.h"
//#include "phlex_solver.h"
//#include "phlex_common.h"

/* SUNDIALS Header files with a description of contents used */
#ifdef PMC_USE_SUNDIALS
#include <cvode/cvode.h>                 /* Protoypes for CVODE fcts., consts.  */
#include <cvode/cvode_impl.h>            /* CVodeMem structure                  */
#include <cvode/cvode_direct.h>          /* CVDls interface                     */
#include <nvector/nvector_serial.h>      /* Serial N_Vector types, fcts, macros */
#include <sundials/sundials_math.h>      /* SUNDIALS math function macros       */
#include <sundials/sundials_types.h>     /* definition of types                 */
#include <sunlinsol/sunlinsol_klu.h>     /* KLU SUNLinearSolver                 */
#include <sunmatrix/sunmatrix_sparse.h>  /* sparse SUNMatrix                    */
//#include <stdbool.h> //This produces segmentation fault on comment typedef
#endif

/* Math constants */
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Number of environmental parameters */
#define PMC_NUM_ENV_PARAM_ 2 // !!! Must match the value in phlex_state.f90 !!!

//Value to consider data size too big -> Memory optimization will change below and under the limit
#define DATA_SIZE_LIMIT_OPT 2000

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define HANDLE_ERROR2( ) (HandleError2( __FILE__, __LINE__ ))

//Todo: use only one common structure fixing bug with enum bool

/* Model data structure */
typedef struct {
    int n_state_var;      // number of state variables per grid cell
    int n_dep_var;        // number of solver variables per grid cell
    int n_jac_elem;       // number of potentially non-zero Jacobian elements
    // per grid cell
    int n_cells;          // Number of cells to compute simultaneously
    double *abs_tol;      // pointer to array of state variable absolute
    // integration tolerances
    int *var_type;	// pointer to array of state variable types (solver,
    // constant, PSSA)
#ifdef PMC_USE_SUNDIALS
    SUNMatrix J_init;	// sparse Jacobian matrix with used elements
                        // initialized to 1.0
#endif
    double *state;	// Pointer to the state array
    double *env;		// Pointer to the environmental state array
    double *rate_constants; // Pointer to the rate constants state array
    void *rxn_data;	// Pointer to reaction parameters
    void *nxt_rxn;	// Pointer to element of rxn_data in which to store next
    // set of reaction data
    void *aero_phase_gpu_data;// Pointer to aerosol phase parameters
    void *nxt_aero_phase; // Pointer to element of aero_phase_gpu_data in which to store
    // the next set of aerosol phase data
    void *aero_rep_gpu_data;	// Pointer to aerosol representation parameters
    void *nxt_aero_rep;	// Pointer to element of aero_rep_gpu_data in which to store
    // the next set of aerosol representation data
    void *sub_model_gpu_data; // Pointer to the sub model parameters
    void *nxt_sub_model;  // Pointer to the element of sub_model_gpu_data in which to
    // store the next set of sub model data
} ModelDatagpu;

/* Solver data structure */
typedef struct {
#ifdef PMC_USE_SUNDIALS
    N_Vector abs_tol_nv;  // abosolute tolerance vector
  N_Vector y;		// vector of solver variables
  SUNLinearSolver ls;   // linear solver
  N_Vector deriv;       // used to calculate the derivative outside the solver
  SUNMatrix J;          // Jacobian matrix
  SUNMatrix J_guess;    // Jacobian matrix for improving guesses sent to linear
                        // solver
  bool curr_J_guess;    // Flag indicating the Jacobian used by the guess helper
                        // is current
  realtype J_guess_t;   // Last time (t) for which J_guess was calculated
  int Jac_eval_fails;   // Number of Jacobian evaluation failures
#ifdef PMC_DEBUG
  booleantype debug_out;// Output debugging information during solving
  booleantype eval_Jac; // Evalute Jacobian data during solving
#endif
#endif
    void *cvode_mem;	// CVodeMem object
    ModelDatagpu model_data; // Model data (used during initialization and solving)
    bool no_solve;        // Flag to indicate whether to run the solver needs to be
    // run. Set to true when no reactions are present.
    double init_time_step;// Initial time step (s)
} SolverDatagpu;


void solver_new_gpu_cu(int n_dep_var,
     int n_state_var, int n_rxn,
     int n_rxn_int_param, int n_rxn_float_param, int n_cells);
void rxn_update_env_state_gpu(ModelDatagpu *model_data, double *env);
void solveRxncpu(ModelDatagpu *model_data, double *deriv_data,
                 double time_step, int *int_data, double *float_data, int deriv_length, int n_rxn);
void rxn_calc_deriv_gpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step);
void rxn_calc_jac_gpu(ModelDatagpu *model_data, SUNMatrix jac, realtype time_step);
void free_gpu_cu();
void print_gpu_specs();
void solver_set_data_gpu(ModelDatagpu *model_data);

#endif
