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

// State variable types (Must match parameters defined in pmc_chem_spec_data
// module)
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

/* Math constants */
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define SMALL 1.0e-30
#define TINY 1.0e-60
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* boolean definition */
typedef enum {false, true} bool;

/* Jacobian map */
typedef struct {
  int solver_id; // solver Jacobian id
  int rxn_id;    // reaction Jacobian id
  int param_id;  // sub model Jacobian id
} JacMap;

/* Model data structure */
typedef struct {
  int n_state_var;	        // number of state variables (>=NV_LENGTH_S(y))
  double *abs_tol;              // pointer to array of state variable absolute
                                // integration tolerances
  int *var_type;	        // pointer to array of state variable types (solver,
                                // constant, PSSA)
#ifdef PMC_USE_SUNDIALS
  SUNMatrix J_init; // sparse solver Jacobian matrix with used elements
                    // initialized to 1.0
  SUNMatrix J_rxn;  // Matrix for Jacobian contributions from reactions
  SUNMatrix J_params;           // Matrix for Jacobian contributions from sub model
                                // parameter calculations
#endif
  JacMap *jac_map;         // Array of Jacobian mapping elements
  JacMap *jac_map_params;  // Array of Jacobian mapping elements to account for
                           // sub-model interdependence. If sub-model parameter
                           // i_dep depends on sub-model parameter i_ind, and
                           // j_ind is a dependency (variable or parameter) of
                           // i_ind, then:
                           // solver_id = jac_id[i_dep][j_ind]
                           // rxn_id    = jac_id[i_dep][i_ind]
                           // param_id  = jac_id[i_ind][j_ind]
  int n_mapped_values;     // Number of Jacobian map elements
  int n_mapped_params;     // Number of Jacobian map elements for sub models
  double *state;           // State array
  double *env;             // Environmental state array

  int n_added_rxns;              // The number of reactions whose data has been
                                 // added to the reaction data arrays
  int *rxn_int_data;             // Pointer to the reaction integer parameters
  double *rxn_float_data;        // Pointer to the reaction floating-point
                                 // parameters
  int *nxt_rxn_int;              // Pointer to the next available integer in
                                 // rxn_int_data
  double *nxt_rxn_float;         // Pointer to the next available floating-point
                                 // number in rxn_float_data
  int **rxn_int_ptrs;            // Array of pointers to integer data for each
                                 // reaction
  double **rxn_float_ptrs;       // Array of pointers to floating-point data for
                                 // each reaction
  int n_added_aero_phases;       // The number of aerosol phases whose data has
                                 // been added to the aerosol phase data arrays
  int *aero_phase_int_data;      // Pointer to the aerosol phase integer parameters
  double *aero_phase_float_data; // Pointer to the aerosol phase floating-point
                                 // parameters
  int *nxt_aero_phase_int;       // Pointer the next available integer in
                                 // aero_phase_int_data
  double *nxt_aero_phase_float;  // Pointer to the next available floating-point
                                 // number in aero_phase_float_data
  int **aero_phase_int_ptrs;     // Array of pointers to the integer data for
                                 // each aerosol phase
  double **aero_phase_float_ptrs;// Array of pointers to the floating-point data
                                 // for each aerosol phase
  int n_added_aero_reps;         // The number of aerosol representations whose
                                 // data has been added to the aerosol
                                 // representation data arrays
  int *aero_rep_int_data;        // Pointer to the aerosol representation integer
                                 // parameters
  double *aero_rep_float_data;   // Pointer to the aerosol representation
                                 // floating-point parameters
  int *nxt_aero_rep_int;         // Pointer to the next available integer in
                                 // aero_rep_int_data
  double *nxt_aero_rep_float;    // Pointer to the next available floating-point
                                 // number in aero_rep_float_data
  int **aero_rep_int_ptrs;       // Array of pointers to integer data for each
                                 // aerosol representation
  double **aero_rep_float_ptrs;  // Array of pointers to floating-point data for
                                 // each aerosol representation
  int n_added_sub_models;        // The number of sub models whose data has been
                                 // added to the sub model data arrays
  int *sub_model_int_data;       // Pointer to sub model integer parameters
  double *sub_model_float_data;  // Pointer to sub model floating-point parameters
  int *nxt_sub_model_int;        // Pointer to the next available integer in
                                 // sub_model_int_data
  double *nxt_sub_model_float;   // Pointer to the next available floating-point
                                 // number in sub_model_float_data
  int **sub_model_int_ptrs;      // Array of pointers to integer data for each
                                 // sub model
  double **sub_model_float_ptrs; // Array of pointers to floating-point data for
                                 // each sub model
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
