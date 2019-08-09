/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 */
/** \file
 * \brief c function tests for the ZSR aerosol water sub model
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../test_common.h"
#include "../../src/sub_model_solver.h"
#include "../../src/sub_models.h"
#include "../../src/phlex_common.h"

// index for the test aerosol representation
#define AERO_REP_IDX 1

// index for the test phase
#define AERO_PHASE_IDX 0

// Number of state variables in the test
#define N_STATE_VAR 13

// Number of Jacobian elements used
// (gas-phase water, Cl_m, Ca_pp for each phase)
#define N_JAC_ELEM 9

// Test concentrations
#define CONC_H2O_G 8000.0
#define CONC_NA    2.5
#define CONC_CL    5.3
#define CONC_CA    1.3


/** \brief Run c function tests
 *
 * \param solver_data Pointer to the solver data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \return 0 if tests pass; otherwise the number of test failures
 */
int run_sub_model_zsr_c_tests(void *solver_data, double *state, double *env)
{
  int ret_val = 0;

#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  int n_solver_var = NV_LENGTH_S(sd->y);
  N_Vector solver_state = N_VNew_Serial(n_solver_var);
  SUNMatrix J = md->J_params;

  // solver state is the same as the model state except for the
  // three aerosol water species
  ret_val += ASSERT_MSG(n_solver_var==N_STATE_VAR-3,
                        "Bad number of state variables");

  md->state = state;

  ret_val += ASSERT_MSG(SM_NNZ_S(J)==N_JAC_ELEM,
                        "Wrong number of flagged Jac elements");

  // Set species concentrations
  NV_DATA_S(solver_state)[0] = state[0]  = CONC_H2O_G;
  NV_DATA_S(solver_state)[1] = state[1]  = CONC_NA;
  NV_DATA_S(solver_state)[2] = state[2]  = CONC_CL;
  NV_DATA_S(solver_state)[3] = state[3]  = CONC_CA;
                               state[4]  = 0.0;
  NV_DATA_S(solver_state)[4] = state[5]  = CONC_NA;
  NV_DATA_S(solver_state)[5] = state[6]  = CONC_CL;
  NV_DATA_S(solver_state)[6] = state[7]  = CONC_CA;
                               state[8]  = 0.0;
  NV_DATA_S(solver_state)[7] = state[9]  = CONC_NA;
  NV_DATA_S(solver_state)[8] = state[10] = CONC_CL;
  NV_DATA_S(solver_state)[9] = state[11] = CONC_CA;
                               state[12] = 0.0;

  // Update the environmental and concentration states
  sub_model_update_env_state(md, env);

  // Run the calculation tests


  // free allocated memory
  N_VDestroy(solver_state);
#endif

  return ret_val;
}

