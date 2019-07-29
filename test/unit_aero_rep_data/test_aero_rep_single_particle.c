/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 */
/** \file
 * \brief c function tests for the single particle aerosol representation
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../test_common.h"
#include "../../src/aero_rep_solver.h"
#include "../../src/aero_reps.h"
#include "../../src/phlex_common.h"

// index for the test aerosol representation
#define AERO_REP_IDX 0

// index for the test phase
#define AERO_PHASE_IDX 1 // (phase 2)

// number of Jacobian elements used for the test phase
#define N_JAC_ELEM 8

// Test concentrations (ug/m3)
#define CONC_1A 1.0
#define CONC_1B 2.0
#define CONC_1C 3.0
#define CONC_2C 4.0
#define CONC_2D 5.0
#define CONC_2E 6.0
#define CONC_3B 7.0
#define CONC_3E 8.0

// Molecular weight of test species (must match json file)
#define MW_A 1.0
#define MW_B 11.0
#define MW_C 36.2
#define MW_D 42.1
#define MW_E 52.3
#define MW_F 623.2
#define MW_G 72.3

// Density of test species (must match json file)
#define DENSITY_A 1.0
#define DENSITY_B 2.0
#define DENSITY_C 3.0
#define DENSITY_D 4.0
#define DENSITY_E 5.0
#define DENSITY_F 6.0
#define DENSITY_G 7.0

/** \brief Run c function tests
 *
 * \param solver_data Pointer to the solver data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \return 0 if tests pass; otherwise number of test failures
 */
int run_aero_rep_single_particle_c_tests(void *solver_data, double *state, double *env) {

  int ret_val = 0;

  SolverData *sd = (SolverData*) solver_data;
  ModelData * model_data = &(sd->model_data);
  void * aero_rep_data = model_data->aero_rep_data;
  int n_solver_var = NV_LENGTH_S(sd->y);
  N_Vector solver_state = N_VNew_Serial(n_solver_var);

  model_data->state = state;

  bool *jac_struct = malloc(sizeof(bool) * n_solver_var);
  ret_val += ASSERT_MSG(jac_struct!=NULL, "jac_struct not allocated");
  if (ret_val>0) return ret_val;

  int aero_phase_idx = AERO_PHASE_IDX;   // phase 2
  int aero_rep_idx   = AERO_REP_IDX;     // only one aero rep in the test

  int n_jac_elem = aero_rep_get_used_jac_elem(model_data, aero_rep_idx,
                       aero_phase_idx, jac_struct);

  ret_val += ASSERT_MSG(n_jac_elem==N_JAC_ELEM, "Bad number of Jac elements");

  // tests are for bin 2
  NV_DATA_S(solver_state)[0] = state[0] = CONC_1A; // phase one, species a
  NV_DATA_S(solver_state)[1] = state[1] = CONC_1B; // phase one, species a
  NV_DATA_S(solver_state)[2] = state[2] = CONC_1C; // phase one, species a
  NV_DATA_S(solver_state)[3] = state[3] = CONC_2C; // phase one, species a
  NV_DATA_S(solver_state)[4] = state[4] = CONC_2D; // phase one, species a
  NV_DATA_S(solver_state)[5] = state[5] = CONC_2E; // phase one, species a
  NV_DATA_S(solver_state)[6] = state[6] = CONC_3B; // phase one, species a
  NV_DATA_S(solver_state)[7] = state[7] = CONC_3E; // phase one, species a

  // Update the environmental and concentration states
  aero_rep_update_env_state(model_data, env);
  aero_rep_update_state(model_data);

  // Run the property tests

  return ret_val;
}
