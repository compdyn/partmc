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

// Externally set properties
#define PART_NUM_CONC 1.23e3
#define PART_RADIUS 2.43e-7

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

  ret_val += ASSERT_MSG(fabs(eff_rad-PART_RADIUS) < 1.0e-6*PART_RADIUS,
                        "Bad effective radius");

  ret_val += ASSERT_MSG(partial_deriv[0] = 999.9,
                        "Bad Jacobian (-1)");
  for( int i = 1; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] = 999.9,
                        "Bad Jacobian (end+1)");

  return ret_val;
}

/** \brief Test the number concentration function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_number_concentration(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double num_conc = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;

  aero_rep_get_number_conc(model_data, AERO_REP_IDX,
                           AERO_PHASE_IDX, &num_conc, &(partial_deriv[1]));

  ret_val += ASSERT_MSG(fabs(num_conc-PART_NUM_CONC) < 1.0e-10*PART_NUM_CONC,
                        "Bad number concentration");

  ret_val += ASSERT_MSG(partial_deriv[0] = 999.9,
                        "Bad Jacobian (-1)");
  for( int i = 1; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] = 999.9,
                        "Bad Jacobian (end+1)");

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

  double mass = CONC_2C + CONC_2D + CONC_2E;

  ret_val += ASSERT_MSG(fabs(phase_mass-mass) < 1.0e-10*mass,
                        "Bad aerosol phase mass");

  ret_val += ASSERT_MSG(partial_deriv[0] = 999.9,
                        "Bad Jacobian (-1)");
  for( int i = 1; i < 4; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  for( int i = 4; i < 7; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ONE,
                          "Bad Jacobian element");
  for( int i = 7; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] = 999.9,
                        "Bad Jacobian (end+1)");

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

  double mass = CONC_2C + CONC_2D + CONC_2E;
  double moles = CONC_2C / MW_C + CONC_2D / MW_D + CONC_2E / MW_E;
  double avg_mw_real = mass / moles;
  double dMW_dC = ONE / moles - mass / (moles * moles * MW_C);
  double dMW_dD = ONE / moles - mass / (moles * moles * MW_D);
  double dMW_dE = ONE / moles - mass / (moles * moles * MW_E);

  ret_val += ASSERT_MSG(fabs(avg_mw-avg_mw_real) < 1.0e-10*avg_mw_real,
                        "Bad average MW");

  ret_val += ASSERT_MSG(partial_deriv[0] = 999.9,
                        "Bad Jacobian (-1)");
  for( int i = 1; i < 4; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(fabs(partial_deriv[4]-dMW_dC) < 1.0e-10*fabs(dMW_dC),
                        "Bad Jacobian (-1)");
  ret_val += ASSERT_MSG(fabs(partial_deriv[5]-dMW_dD) < 1.0e-10*fabs(dMW_dD),
                        "Bad Jacobian (-1)");
  ret_val += ASSERT_MSG(fabs(partial_deriv[6]-dMW_dE) < 1.0e-10*fabs(dMW_dE),
                        "Bad Jacobian (-1)");
  for( int i = 7; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] = 999.9,
                        "Bad Jacobian (end+1)");

  return ret_val;
}

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
  ret_val += test_effective_radius(model_data, solver_state);
  ret_val += test_aero_phase_mass(model_data, solver_state);
  ret_val += test_aero_phase_avg_MW(model_data, solver_state);
  ret_val += test_number_concentration(model_data, solver_state);

  return ret_val;
}
