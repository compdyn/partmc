
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

//TTODO: Try change max shared memory per block allowed on cuda config
// with cudaDeviceSetCacheConfig, maybe we can reach improvement using 16kb
//instead of default 48kb

#define MAX_N_GPU_THREAD 1024
#define MAX_SHARED_MEMORY_BLOCK_DOUBLE 1000

//#define MAX_N_GPU_BLOCK 10

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

/* Model data structure */
typedef struct {
    int n_state_var;	// number of state variables (>=NV_LENGTH_S(y))
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
    bool use_adj;         // Flag to indicate whether state adjustments exist
    double *state_adj;    // Adjustments to the state array applied prior to
    // calculating rates, derivatives, etc. Used for fast
    // reactions that essentially go to completion during
    // the solver timestep.
    bool scale_adj;       // Flag to indicate state adjustments in state_adj should
    // be scaled by relative contributions from multiple rxns.
    double *rel_flux;     // Used to calculate relative contributions of each rxn
    // to state adjustments. (For scaling when more than one
    // rxn rapidly depletes the same species.)
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
#endif
    void *cvode_mem;	// CVodeMem object
    ModelDatagpu model_data; // Model data (used during initialization and solving)
    bool no_solve;        // Flag to indicate whether to run the solver needs to be
    // run. Set to true when no reactions are present.
    double init_time_step;// Initial time step (s)
} SolverDatagpu;


//void printfCPP();
void printfCUDA(int aggg);
void solver_new_gpu_cu(SolverDatagpu *sd, int n_dep_var,
     int n_state_var, int *var_type, int n_rxn,
     int n_rxn_int_param, int n_rxn_float_param, int n_aero_phase,
     int n_aero_phase_int_param, int n_aero_phase_float_param,
     int n_aero_rep, int n_aero_rep_int_param, int n_aero_rep_float_param,
     int n_sub_model, int n_sub_model_int_param, int n_sub_model_float_param);
void solver_update_gpu(ModelDatagpu *md);
void rxn_calc_deriv_gpu_cu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step);
void free_gpu_cu();

#endif
