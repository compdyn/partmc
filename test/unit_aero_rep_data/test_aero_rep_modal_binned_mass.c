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
#include "../../src/camp_common.h"

// index for the test aerosol representation
#define AERO_REP_IDX 0

// index for test phase
#define AERO_PHASE_IDX 9   // (bin 4 phase 1)
#define AERO_PHASE_IDX_2 1 // (mode 1 phase 2)

// number of Jacobian elements used for test phase
#define N_JAC_ELEM 5   // (bin 4 phase 1)
#define N_JAC_ELEM_2 6 // (mode 1 phase 2)

// Test concentrations (kg/m3)
// (bin 4)
#define CONC_1A 1.0
#define CONC_1B 2.0
#define CONC_1C 3.0
#define CONC_3B 4.0
#define CONC_3E 5.0
// (mode 1)
#define CONC_2_1A 6.0
#define CONC_2_1B 7.0
#define CONC_2_1C 8.0
#define CONC_2_2C 9.0
#define CONC_2_2D 10.0
#define CONC_2_2E 11.0

// Molecular weight of test species (must match json file)
#define MW_A 11.2
#define MW_B 21.2
#define MW_C 31.2
#define MW_D 41.2
#define MW_E 51.2
#define MW_F 61.2

// Density of test species (must match json file)
#define DENSITY_A 1.0
#define DENSITY_B 2.0
#define DENSITY_C 3.0
#define DENSITY_D 4.0
#define DENSITY_E 5.0
#define DENSITY_F 6.0

/** \brief Test the effective radius function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
#ifdef PMC_USE_SUNDIALS
int test_effective_radius(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double partial_deriv_2[N_JAC_ELEM_2+2];
  double eff_rad = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;
  for( int i = 0; i < N_JAC_ELEM_2+2; ++i ) partial_deriv_2[i] = 999.9;

  aero_rep_get_effective_radius__m(model_data, AERO_REP_IDX,
                                AERO_PHASE_IDX, &eff_rad, &(partial_deriv[1]));

  double dp_bin4 = pow(10.0,(log10(1.0e-6) - log10(8.0e-9)) / 7.0 * 3.0 + log10(8.0e-9));
  double real_rad = dp_bin4 / 2.0;
  ret_val += ASSERT_MSG(fabs(eff_rad-real_rad)<1.0e-10*real_rad,
                        "Bad effective radius");

  double real_rad_2 = 1.2e-6 / 2.0 * exp(2.0 * log(1.2) * log(1.2));
  aero_rep_get_effective_radius__m(model_data, AERO_REP_IDX,
                                AERO_PHASE_IDX_2, &eff_rad, &(partial_deriv_2[1]));
  ret_val += ASSERT_MSG(fabs(eff_rad-real_rad_2)<1.0e-10*real_rad_2,
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

/** \brief Test the number concentration function
 *
 * \param model_data Pointer to the model data
 * \param state Solver state
 */
int test_number_conc(ModelData * model_data, N_Vector state) {

  int ret_val = 0;
  double partial_deriv[N_JAC_ELEM+2];
  double partial_deriv_2[N_JAC_ELEM_2+2];
  double number_conc = -999.9;
  double number_conc_2 = -999.9;

  for( int i = 0; i < N_JAC_ELEM+2; ++i ) partial_deriv[i] = 999.9;
  for( int i = 0; i < N_JAC_ELEM_2+2; ++i ) partial_deriv_2[i] = 999.9;

  aero_rep_get_number_conc__n_m3(model_data, AERO_REP_IDX, AERO_PHASE_IDX,
                           &number_conc, &(partial_deriv[1]));
  aero_rep_get_number_conc__n_m3(model_data, AERO_REP_IDX, AERO_PHASE_IDX_2,
                           &number_conc_2, &(partial_deriv_2[1]));

  double dp_bin4 = pow(10.0,(log10(1.0e-6) - log10(8.0e-9)) / 7.0 * 3.0 + log10(8.0e-9));
  double vp_bin4  = 4.0/3.0*M_PI * pow(dp_bin4/2.0, 3.0);

  double real_number_conc = (CONC_1A / DENSITY_A +
                             CONC_1B / DENSITY_B +
                             CONC_1C / DENSITY_C +
                             CONC_3B / DENSITY_B +
                             CONC_3E / DENSITY_E) / vp_bin4;

  double vp_mode1 = M_PI/6.0 * pow(1.2e-6, 3.0) * exp(9.0/2.0 * log(1.2) * log(1.2));

  double real_number_conc_2 = (CONC_2_1A / DENSITY_A +
                               CONC_2_1B / DENSITY_B +
                               CONC_2_1C / DENSITY_C +
                               CONC_2_2C / DENSITY_C +
                               CONC_2_2D / DENSITY_D +
                               CONC_2_2E / DENSITY_E) / vp_mode1;

  ret_val += ASSERT_MSG(fabs(number_conc-real_number_conc) <
                        1.0e-10 * real_number_conc,
                        "Bad number concentration");
  ret_val += ASSERT_MSG(fabs(number_conc_2-real_number_conc_2) <
                        1.0e-10 * real_number_conc_2,
                        "Bad number concentration");

  double real_partial;

  // (bin 4)
  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  real_partial = ONE / DENSITY_A / vp_bin4;
  ret_val += ASSERT_MSG(fabs(partial_deriv[1]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_B / vp_bin4;
  ret_val += ASSERT_MSG(fabs(partial_deriv[2]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_C / vp_bin4;
  ret_val += ASSERT_MSG(fabs(partial_deriv[3]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_B / vp_bin4;
  ret_val += ASSERT_MSG(fabs(partial_deriv[4]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_E / vp_bin4;
  ret_val += ASSERT_MSG(fabs(partial_deriv[5]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] == 999.9,
                        "Bad Jacobian index (-1)");

  // (mode 1)
  ret_val += ASSERT_MSG(partial_deriv_2[0] == 999.9,
                        "Bad Jacobian index (-1)");
  real_partial = ONE / DENSITY_A / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[1]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_B / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[2]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_C / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[3]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_C / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[4]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_D / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[5]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  real_partial = ONE / DENSITY_E / vp_mode1;
  ret_val += ASSERT_MSG(fabs(partial_deriv_2[6]-real_partial) <
                        1.0e-10 * fabs(real_partial),
                        "Bad Jacobian element");
  ret_val += ASSERT_MSG(partial_deriv_2[N_JAC_ELEM_2+1] == 999.9,
                        "Bad Jacobian index (-1)");

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

  aero_rep_get_aero_phase_mass__kg_m3(model_data, AERO_REP_IDX, AERO_PHASE_IDX,
                               &phase_mass, &(partial_deriv[1]));

  double mass = CONC_1A + CONC_1B + CONC_1C;
  ret_val += ASSERT_MSG(fabs(phase_mass-mass) < 1.0e-10*mass, "Bad phase mass");

  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  for( int i = 1; i < 4; ++i ) {
    ret_val += ASSERT_MSG(partial_deriv[i] == ONE,
                          "Bad Jacobian element");
  }
  for( int i = 4; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                        "Bad Jacobian index (end+1)");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] == 999.9,
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

  aero_rep_get_aero_phase_avg_MW__kg_mol(model_data, AERO_REP_IDX, AERO_PHASE_IDX,
                                 &avg_mw, &(partial_deriv[1]));

  ret_val += ASSERT_MSG(fabs(avg_mw-21.44548402) < 1.0e-8, "Bad phase avg MW");

  // MW = mass_total / moles_total
  // d_MW / d_y = 1 / moles_total - mass_total / ( moles_total^2 * MW_y )
  double mass = CONC_1A + CONC_1B + CONC_1C;
  double moles = CONC_1A / MW_A + CONC_1B / MW_B + CONC_1C / MW_C;
  double dMW_dA = ONE / moles - mass / (moles * moles * MW_A);
  double dMW_dB = ONE / moles - mass / (moles * moles * MW_B);
  double dMW_dC = ONE / moles - mass / (moles * moles * MW_C);

  ret_val += ASSERT_MSG(partial_deriv[0] == 999.9,
                        "Bad Jacobian index (-1)");
  ret_val += ASSERT_MSG(fabs(partial_deriv[1]-dMW_dA) < fabs(1.0e-10*dMW_dA),
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(fabs(partial_deriv[2]-dMW_dB) < fabs(1.0e-10*dMW_dB),
                          "Bad Jacobian element");
  ret_val += ASSERT_MSG(fabs(partial_deriv[3]-dMW_dC) < fabs(1.0e-10*dMW_dC),
                          "Bad Jacobian element");
  for( int i = 4; i < N_JAC_ELEM+1; ++i )
    ret_val += ASSERT_MSG(partial_deriv[i] == ZERO,
                        "Bad Jacobian index");
  ret_val += ASSERT_MSG(partial_deriv[N_JAC_ELEM+1] == 999.9,
                        "Bad Jacobian index (end+1)");

  return ret_val;
}
#endif

/** \brief Run c function tests
 *
 * \param solver_data Pointer to solver data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \return 0 if tests pass; otherwise number of test failures
 */
int run_aero_rep_modal_c_tests(void *solver_data, double *state, double *env) {

  int ret_val = 0;

#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData*) solver_data;
  ModelData * model_data = &(sd->model_data);
  int n_solver_var = NV_LENGTH_S(sd->y);
  N_Vector solver_state = N_VNew_Serial(n_solver_var);

  model_data->grid_cell_id    = 0;
  model_data->total_state     = state;
  model_data->grid_cell_state = model_data->total_state;
  model_data->total_env       = env;
  model_data->grid_cell_env   = model_data->total_env;

  bool *jac_struct = malloc(sizeof(bool) * n_solver_var);
  ret_val += ASSERT_MSG(jac_struct!=NULL, "jac_struct not allocated");
  if (ret_val>0) return ret_val;

  int aero_phase_idx   = AERO_PHASE_IDX;   // bin 4 phase one
  int aero_phase_idx_2 = AERO_PHASE_IDX_2; // mode 1 phase two
  int aero_rep_idx     = AERO_REP_IDX;     // only one aero rep in the test

  int n_jac_elem   = aero_rep_get_used_jac_elem(model_data, aero_rep_idx,
                        aero_phase_idx, jac_struct);
  int n_jac_elem_2 = aero_rep_get_used_jac_elem(model_data, aero_rep_idx,
                        aero_phase_idx_2, jac_struct);
  free(jac_struct);

  ret_val += ASSERT_MSG(n_jac_elem==N_JAC_ELEM, "Bad number of Jac elements");
  ret_val += ASSERT_MSG(n_jac_elem_2==N_JAC_ELEM_2, "Bad number of Jac elements");

  // tests are for bin 4
  NV_DATA_S(solver_state)[17] = state[17] = CONC_1A;   // phase one, species a
  NV_DATA_S(solver_state)[18] = state[18] = CONC_1B;   // phase one, species b
  NV_DATA_S(solver_state)[19] = state[19] = CONC_1C;   // phase one, species c
  NV_DATA_S(solver_state)[38] = state[38] = CONC_3B;   // last phase, species b
  NV_DATA_S(solver_state)[39] = state[39] = CONC_3E;   // last phase, species e
  NV_DATA_S(solver_state)[ 0] = state[ 0] = CONC_2_1A; // phase one, species a
  NV_DATA_S(solver_state)[ 1] = state[ 1] = CONC_2_1B; // phase one, species b
  NV_DATA_S(solver_state)[ 2] = state[ 2] = CONC_2_1C; // phase one, species c
  NV_DATA_S(solver_state)[ 3] = state[ 3] = CONC_2_2C; // phase two, species c
  NV_DATA_S(solver_state)[ 4] = state[ 4] = CONC_2_2D; // phase two, species d
  NV_DATA_S(solver_state)[ 5] = state[ 5] = CONC_2_2E; // phase two, species e

  // Set the environment-dependent parameter pointer to the first grid cell
  model_data->grid_cell_aero_rep_env_data = model_data->aero_rep_env_data;

  // Update the environmental and concentration states
  aero_rep_update_env_state(model_data);
  aero_rep_update_state(model_data);

  // Run the property tests
  ret_val += test_effective_radius(model_data, solver_state);
  ret_val += test_aero_phase_mass(model_data, solver_state);
  ret_val += test_aero_phase_avg_MW(model_data, solver_state);
  ret_val += test_number_conc(model_data, solver_state);

  N_VDestroy(solver_state);
#endif

  return ret_val;
}
