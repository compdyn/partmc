/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for common constants and structures
 *
*/
/** \file
 * \brief Header file for common constants and structures
*/
#ifndef PHLEX_COMMON_H
#define PHLEX_COMMON_H

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
#endif

/* Math constants */
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define SMALL 1.0e-30
#define TINY 1.0e-60
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* Number of environmental parameters */
#define PMC_NUM_ENV_PARAM_ 2 // !!! Must match the value in phlex_state.f90 !!!

/* boolean definition */
//CUDA/C++ already has bool definition, so is not necessary for GPU files
#ifndef PHLEX_GPU_SOLVER_H_
typedef enum {false, true} bool;
#endif

/* Model data structure */
typedef struct {
  int n_state_var;      // number of state variablesper grid cell
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
  void *aero_phase_data;// Pointer to aerosol phase parameters
  void *nxt_aero_phase; // Pointer to element of aero_phase_data in which to store
                        // the next set of aerosol phase data
  void *aero_rep_data;	// Pointer to aerosol representation parameters
  void *nxt_aero_rep;	// Pointer to element of aero_rep_data in which to store
  			// the next set of aerosol representation data
  void *sub_model_data; // Pointer to the sub model parameters
  void *nxt_sub_model;  // Pointer to the element of sub_model_data in which to
                        // store the next set of sub model data
} ModelData;

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
  ModelData model_data; // Model data (used during initialization and solving)
  bool no_solve;        // Flag to indicate whether to run the solver needs to be
                        // run. Set to true when no reactions are present.
  double init_time_step;// Initial time step (s)
} SolverData;

#endif
