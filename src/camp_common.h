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
#ifndef CAMP_COMMON_H
#define CAMP_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Jacobian.h"
#include "time_derivative.h"

// todo cant found cvode_impl in docker, put under #gpu flags or something

/* SUNDIALS Header files with a description of contents used */
#ifdef PMC_USE_SUNDIALS
#include <cvode/cvode.h>        /* Protoypes for CVODE fcts., consts.  */
#include <cvode/cvode_direct.h> /* CVDls interface                     */
#include <cvode/cvode_impl.h>   /* CVodeMem structure                  */

// todo ifndef GPU
//#ifdef PMC_USE_GPU
//  #include <nvector/nvector_cuda.h>      /* access to cuda N_Vector */
//  #include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver */
//#else
#include <nvector/nvector_serial.h>  /* Serial N_Vector types, fcts, macros */
#include <sunlinsol/sunlinsol_klu.h> /* KLU SUNLinearSolver                 */
//#endif

#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h> /* sparse SUNMatrix                    */

#include <sundials/sundials_math.h>  /* SUNDIALS math function macros       */
#include <sundials/sundials_types.h> /* definition of types                 */

#ifdef PMC_USE_GPU
#include <cvode/cvode_direct_impl.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#endif

#endif

#ifdef PMC_USE_MPI
#include <mpi.h>
#endif

#ifdef PMC_USE_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuda/cuda_structs.h"
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

/* Number of environmental parameters */
#define PMC_NUM_ENV_PARAM_ 2  // !!! Must match the value in camp_state.f90 !!!

/* boolean definition */
// CUDA/C++ already has bool definition: Avoid issues disabling it for GPU
#ifndef CAMP_GPU_SOLVER_H_
#ifndef CVODE_gpu2_SOLVER_H_
#ifndef CAMP_GPU_CUSOLVER_H_
#ifndef ITSOLVERGPU_H
typedef enum { false, true } bool;
#endif
#endif
#endif
#endif

/* Jacobian map */
typedef struct {
  int solver_id;  // solver Jacobian id
  int rxn_id;     // reaction Jacobian id
  int param_id;   // sub model Jacobian id
} JacMap;

/* Model data structure */
typedef struct {
  int n_per_cell_state_var;        // number of state variables per grid cell
  int n_per_cell_dep_var;          // number of solver variables per grid cell
  int n_per_cell_rxn_jac_elem;     // number of potentially non-zero
                                   // reaction Jacobian elements
  int n_per_cell_param_jac_elem;   // number of potentially non-zero
                                   // parameter Jacobian elements
  int n_per_cell_solver_jac_elem;  // number of potentially non-zero
                                   // solver Jacobian elements
  int n_cells;                     // number of cells to compute simultaneously
  double *abs_tol;  // pointer to array of state variable absolute
                    // integration tolerances
  int *var_type;    // pointer to array of state variable types (solver,
                    // constant, PSSA)
#ifdef PMC_USE_SUNDIALS
  SUNMatrix J_init;    // sparse solver Jacobian matrix with used elements
                       // initialized to 1.0
  SUNMatrix J_init2;   // sparse solver Jacobian matrix with used elements
  SUNMatrix J_rxn;     // Matrix for Jacobian contributions from reactions
  SUNMatrix J_params;  // Matrix for Jacobian contributions from sub model
                       // parameter calculations
  SUNMatrix J_solver;  // Solver Jacobian
  N_Vector J_state;    // Last state used to calculate the Jacobian
  N_Vector J_deriv;    // Last derivative used to calculate the Jacobian
  N_Vector J_tmp;      // Working vector (size of J_state and J_deriv)
  N_Vector J_tmp2;     // Working vector (size of J_state and J_deriv)
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

  int grid_cell_id;         // Index of the current grid cell
  double *grid_cell_state;  // Pointer to the current grid cell being solved
                            // on the total_state array
  double *total_state;      // Total (multi-cell) state array
  double *grid_cell_env;    // Pointer to the current grid cell being solved
                            // on the total_env state array
  double *total_env;        // Total (multi-cell) environmental state array
  double *grid_cell_rxn_env_data;  // Environment-dependent parameters for the
                                   // current grid cell
  double *rxn_env_data;            // Total (multi-cell) reaction environment-
                                   // dependent parameters
  double *grid_cell_aero_rep_env_data;
  // Environment-dependent parameters for the
  // current grid cell
  double *aero_rep_env_data;  // Total (multi-cell) aerosol representation
                              // environment-dependent parameters
  double *grid_cell_sub_model_env_data;
  // Environment-dependent parameters for the
  // current grid cell
  double *sub_model_env_data;  // Total (multi-cell) sub-model environment-
                               // dependent parameters

  int n_rxn;                 // Number of reactions
  int n_added_rxns;          // The number of reactions whose data has been
                             // added to the reaction data arrays
  int *rxn_int_data;         // Pointer to the reaction integer parameters
  double *rxn_float_data;    // Pointer to the reaction floating-point
                             // parameters
  int *rxn_int_indices;      // Array of indices of integer data
  int *rxn_float_indices;    // Array of indices of float data
  int *rxn_env_idx;          // Array of offsets for the environment-
                             // dependent data for each reaction from the
                             // beginning of the environmental dependent data
                             // for the current grid cell
  int n_rxn_env_data;        // Number of reaction environmental parameters
                             // from all reactions
  int n_aero_phase;          // Number of aerosol phases
  int n_added_aero_phases;   // The number of aerosol phases whose data has
                             // been added to the aerosol phase data arrays
  int *aero_phase_int_data;  // Pointer to the aerosol phase integer parameters
  double *aero_phase_float_data;  // Pointer to the aerosol phase floating-point
                                  // parameters
  int *aero_phase_int_indices;    // Array of indices of integer data
  int *aero_phase_float_indices;  // Array of indices of float data
  int n_aero_rep;                 // Number of aerosol representations
  int n_added_aero_reps;          // The number of aerosol representations whose
                                  // data has been added to the aerosol
                                  // representation data arrays
  int *aero_rep_int_data;       // Pointer to the aerosol representation integer
                                // parameters
  double *aero_rep_float_data;  // Pointer to the aerosol representation
                                // floating-point parameters
  int *aero_rep_int_indices;    // Array of indices of integer data
  int *aero_rep_float_indices;  // Array of indices of float data
  int *aero_rep_env_idx;        // Array of offsets for the environment-
                          // dependent data for each aerosol representation
                          // from the beginning of the environment-
                          // dependent data for the current grid cell
  int n_aero_rep_env_data;  // Number of aerosol representation environmental
                            // parameters for all aerosol representations
  int n_sub_model;          // Number of sub models
  int n_added_sub_models;   // The number of sub models whose data has been
                            // added to the sub model data arrays
  int *sub_model_int_data;  // Pointer to sub model integer parameters
  double
      *sub_model_float_data;   // Pointer to sub model floating-point parameters
  int *sub_model_int_indices;  // Array of indices of integer data
  int *sub_model_float_indices;  // Array of indices of float data
  int *sub_model_env_idx;        // Array of offsets for the environment-
                                 // dependent data for each sub model from the
                                 // beginning of the environment-dependent data
                                 // for the current grid cell
  int n_sub_model_env_data;      // Number of sub model environmental parameters
                                 // from all sub models
  int counterMD;

  //#ifdef CUDA_FOUND
  // GPU data
// Gpu definitions
#ifdef PMC_USE_GPU
  int max_n_gpu_thread;
  int max_n_gpu_blocks;
  int *index_deriv_state;
  int *index_deriv_state_gpu;
  double *deriv_gpu_data;
  double *deriv_aux;
  double *jac_gpu_data;  // todo set this pointer to bicg.dA and everyone will
                         // be happy
  double *jac_aux;
  int *indexvals_gpu;
  int *indexptrs_gpu;
  size_t deriv_size;
  size_t jac_size;
  size_t state_size;
  size_t env_size;
  size_t rxn_env_data_size;
  size_t rxn_env_data_idx_size;
  size_t index_deriv_state_size;
  int small_data;
  bool implemented_all;
  int *int_pointer_gpu;
  double *double_pointer_gpu;
  double *state_gpu;
  double *env_gpu;
  double *rxn_env_data_gpu;
  int *rxn_env_data_idx_gpu;
  int model_data_id;  // Id of the modelData object
  int counterDeriv2;
  // cudaStream_t *stream_gpu;
#endif

} ModelData;

/* Solver data structure */
typedef struct {
#ifdef PMC_USE_SUNDIALS
  N_Vector abs_tol_nv;  // abosolute tolerance vector
  N_Vector y;           // vector of solver variables
  N_Vector y2;
  SUNLinearSolver ls;         // linear solver
  TimeDerivative time_deriv;  // CAMP derivative structure for use in
                              // calculating deriv
  Jacobian jac;               // CAMP Jacobian structure for use in
                              // calculating the Jacobian
  N_Vector deriv;  // used to calculate the derivative outside the solver
  N_Vector deriv2;
  SUNMatrix J;  // Jacobian matrix
  SUNMatrix J2;
  SUNMatrix J_guess;  // Jacobian matrix for improving guesses sent to linear
                      // solver
  SUNMatrix J_guess2;
  bool curr_J_guess;   // Flag indicating the Jacobian used by the guess helper
                       // is current
  realtype J_guess_t;  // Last time (t) for which J_guess was calculated
  int Jac_eval_fails;  // Number of Jacobian evaluation failures
  // realtype gamma;           // Scale factor jacobian (M=I-gJ)
  int solver_flag;       // Last flag returned by a call to CVode()
  int output_precision;  // Flag indicating whether to output precision loss
  int use_deriv_est;     // Flag indicating whether to use an estimated
                         // derivative in the f() calculations
#ifdef PMC_DEBUG
  booleantype debug_out;  // Output debugging information during solving
  booleantype eval_Jac;   // Evalute Jacobian data during solving
  int counterDeriv;       // Total calls to f()
  int counterJac;         // Total calls to Jac()
  clock_t timeDeriv;      // Compute time for calls to f()
  clock_t timeJac;        // Compute time for calls to Jac()
  double
      max_loss_precision;  // Maximum loss of precision during last call to f()
#endif

  int counterDerivTotal;  // Total calls to f()
  int counterDerivGPU;    // todo set as counterDeriv and fix old counterDeriv
  int counterJacGPU;
  int counterSolve;
  int counterFail;
  double timeCVode;
  double timeCVodeTotal;
  double timeF;

#endif
#ifdef PMC_USE_GPU
  itsolver bicg;
#ifdef CHECK_GPU_LINSOLVE  // todo fix name
  double max_error_linsolver;
  int max_error_linsolver_i;
  int n_linsolver_i;
#endif
#endif
  void *cvode_mem;       // CVodeMem object
  ModelData model_data;  // Model data (used during initialization and solving)
  bool no_solve;  // Flag to indicate whether to run the solver needs to be
                  // run. Set to true when no reactions are present.
  double init_time_step;  // Initial time step (s)
  char **spec_names;      // Species names
} SolverData;

#endif
