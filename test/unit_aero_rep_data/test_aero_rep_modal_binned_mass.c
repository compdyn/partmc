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

// index for the test aerosol representation
#define AERO_REP_IDX 0

// index for test phase (bin 4 phase 1)
#define AERO_PHASE_IDX 9

// number of Jacobian elements used for test phase
#define N_JAC_ELEM 5

// Test concentrations (ug/m3)
#define CONC_1A 1.0
#define CONC_1B 2.0
#define CONC_1C 3.0
#define CONC_3B 4.0
#define CONC_3E 5.0

// Molecular weight of test species (must match json file)
#define MW_1A 11.2
#define MW_1B 21.2
#define MW_1C 31.2

/** \brief Test the effective radius function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_effective_radius(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double eff_rad = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;

  aero_rep_get_effective_radius(model_data, AERO_REP_IDX,
                                AERO_PHASE_IDX, &eff_rad, &(partial_deriv[1]));

  ret_val += ASSERT_MSG(fabs(eff_rad-6.3353e-8)<1.0e-12,
                        "Bad effective radius");

  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  for( int i = 1; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] == 999.9,
                        "Bad Jacobian index (end+1)");

  return ret_val;
}

/** \brief Test the total aerosol phase mass function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_aero_phase_mass(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double phase_mass = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;

  aero_rep_get_aero_phase_mass(model_data, AERO_REP_IDX, AERO_PHASE_IDX,
                               &phase_mass, &(partial_deriv[1]));

  double mass = CONC_1A + CONC_1B + CONC_1C;
  ret_val += ASSERT_MSG(fabs(phase_mass-mass) < 1.0e-10*mass, "Bad phase mass");

  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  for( int i = 1; i < 4; ++i ) {
    ret_val += ASSERT_MSG(partial_deriv[i] == ONE,
                          "Bad Jacobian element");
  }
  for( int i = 4; i < N_JAC_ELEM+2; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == 999.9,
                        "Bad Jacobian index (end+1)");

  return ret_val;
}

/** \brief Test the aerosol phase average molecular weight function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_aero_phase_avg_MW(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double avg_mw = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;

  aero_rep_get_aero_phase_avg_MW(model_data, AERO_REP_IDX, AERO_PHASE_IDX,
                                 &avg_mw, &(partial_deriv[1]));

  ret_val += ASSERT_MSG(fabs(avg_mw-21.44548402) < 1.0e-8, "Bad phase avg MW");

  // MW = mass_total / moles_total
  // d_MW / d_y = 1 / moles_total - mass_total / ( moles_total^2 * MW_y )
  double mass = CONC_1A + CONC_1B + CONC_1C;
  double moles = CONC_1A / MW_1A + CONC_1B / MW_1B + CONC_1C / MW_1C;
  double dMW_dA = ONE / moles - mass / (moles * moles * MW_1A);
  double dMW_dB = ONE / moles - mass / (moles * moles * MW_1B);
  double dMW_dC = ONE / moles - mass / (moles * moles * MW_1C);

  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  ret_val += ASSERT_MSG(fabs(partial_deriv[1]-dMW_dA) < fabs(1.0e-10*dMW_dA),
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(fabs(partial_deriv[2]-dMW_dB) < fabs(1.0e-10*dMW_dB),
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(fabs(partial_deriv[3]-dMW_dC) < fabs(1.0e-10*dMW_dC),
                          "Bad Jacobian element");
  for( int i = 4; i < N_JAC_ELEM+2; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == 999.9,
                        "Bad Jacobian index (end+1)");

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

  model_data->state = state;

  bool *jac_struct = malloc(sizeof(bool) * n_solver_var);
  ret_val += ASSERT_MSG(jac_struct!=NULL, "jac_struct not allocated");
  if (ret_val>0) return ret_val;

  int aero_phase_idx = AERO_PHASE_IDX; // bin 4 phase one
  int aero_rep_idx   = AERO_REP_IDX; // only one aero rep in the test

  int n_jac_elem = aero_rep_get_used_jac_elem(model_data, aero_rep_idx,
                        aero_phase_idx, jac_struct);
  free(jac_struct);

  ret_val += ASSERT_MSG(n_jac_elem==N_JAC_ELEM, "Bad number of Jac elements");

  // tests are for bin 4
  NV_DATA_S(solver_state)[17] = state[17] = CONC_1A; // phase one, species a
  NV_DATA_S(solver_state)[18] = state[18] = CONC_1B; // phase one, species b
  NV_DATA_S(solver_state)[19] = state[19] = CONC_1C; // phase one, species c
  NV_DATA_S(solver_state)[38] = state[38] = CONC_3B; // last phase, species b
  NV_DATA_S(solver_state)[39] = state[39] = CONC_3E; // last phase, species e

  // Update the environmental and concentration states
  aero_rep_update_env_state(model_data, env);
  aero_rep_update_state(model_data);

  // Run the property tests
  ret_val += test_effective_radius(model_data, solver_state);
  ret_val += test_aero_phase_mass(model_data, solver_state);
  ret_val += test_aero_phase_avg_MW(model_data, solver_state);

  return ret_val;
}
