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
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* boolean definition */
typedef enum {false, true} bool;

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
  void *aero_phase_data;// Pointer to aerosol phase parameters
  void *nxt_aero_phase; // Pointer to element of aero_phase_data in which to store
                        // the next set of aerosol phase data
  void *aero_rep_data;	// Pointer to aerosol representation parameters
  void *nxt_aero_rep;	// Pointer to element of aero_rep_data in which to store
  			// the next set of aerosol representation data
  void *sub_model_data; // Pointer to the sub model parameters
  void *nxt_sub_model;  // Pointer to the element of sub_model_data in which to
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
} ModelData;

/* Solver data structure */
typedef struct {
#ifdef PMC_USE_SUNDIALS
  N_Vector abs_tol_nv;  // abosolute tolerance vector
  N_Vector y;		// vector of solver variables
  SUNLinearSolver ls;   // linear solver
  N_Vector deriv;       // used to calculate the derivative outside the solver
  SUNMatrix J;          // Jacobian matrix
#endif
  void *cvode_mem;	// CVodeMem object
  ModelData model_data; // Model data (used during initialization and solving)
  bool no_solve;        // Flag to indicate whether to run the solver needs to be
                        // run. Set to true when no reactions are present.
  double init_time_step;// Initial time step (s)
} SolverData;

#endif
