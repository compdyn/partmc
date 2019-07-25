/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 */
/** \file
 * \brief c function tests for the modal/binned aerosol representation
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../test_common.h"
#include "../../src/aero_rep_solver.h"
#include "../../src/aero_reps.h"
#include "../../src/phlex_common.h"

/** \brief Test the effective radius function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_effective_radius(ModelData * model_data, N_Vector state) {

  int ret_val = 0;

  return ret_val;
}

/** \brief Run c function tests
 *
 * \param solver_data Pointer to solver data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \return 0 if tests pass; otherwise number of test failures
 */
int run_aero_rep_modal_c_tests(void *solver_data, double *state, double *env) {

  int ret_val = 0;

  SolverData *sd = (SolverData*) solver_data;
  ModelData * model_data = &(sd->model_data);
  void * aero_rep_data = model_data->aero_rep_data;
  int n_solver_var = NV_LENGTH_S(sd->y);
  N_Vector solver_state = N_VNew_Serial(n_solver_var);

  bool *jac_struct = malloc(sizeof(bool) * n_solver_var);
  ret_val += PMC_ASSERT_MSG(jac_struct!=NULL, "jac_struct not allocated");
  if (ret_val>0) return ret_val;

  int aero_phase_idx = 5; // bin 4 phase one
  int aero_rep_idx   = 0; // only one aero rep in the test

  int n_jac_elem = aero_rep_get_used_jac_elem(model_data, aero_rep_idx,
                        aero_phase_idx, jac_struct);
  free(jac_struct);

  ret_val += PMC_ASSERT_MSG(n_jac_elem==5, "Bad number of Jac elements");

  // tests are for bin 4
  NV_DATA_S(solver_state)[17] = 1.0; // phase one, species a
  NV_DATA_S(solver_state)[18] = 2.0; // phase one, species b
  NV_DATA_S(solver_state)[19] = 3.0; // phase one, species c
  NV_DATA_S(solver_state)[38] = 4.0; // last phase, species b
  NV_DATA_S(solver_state)[39] = 5.0; // last phase, species e

  ret_val += test_effective_radius(model_data, solver_state);

  return ret_val;
}
