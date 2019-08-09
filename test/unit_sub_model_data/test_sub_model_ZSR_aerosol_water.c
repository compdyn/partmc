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
// (RH ~ 51%)
#define CONC_H2O_G 16000.0
#define CONC_NA    2.5
#define CONC_CL    5.3
#define CONC_CA    1.3

// Molecular weights
#define MW_H2O 0.01801
#define MW_NA 0.0229898
#define MW_CL 0.035453
#define MW_CA 0.040078

/** \brief Check the Jacobian calculations
 *
 * \param solver_data Pointer to the solver data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \param J Jacobian to use for testing
 * \return 0 if tests pass; otherwise the number of test failures
 */
#ifdef PMC_USE_SUNDIALS
int test_sub_model_zsr_jac_calc(void *solver_data, double *state, double *env,
    SUNMatrix J)
{
  int ret_val = 0;
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  double temperature_K = env[0];
  double pressure_Pa   = env[1];

  // Reset the Jacobian
  SUNMatZero(J);

  // variable to hold calculated Jac elements
  double J_H2O = 0.0;
  double J_NA  = 0.0;
  double J_CA  = 0.0;
  double J_CL  = 0.0;

  // ppm -> RH conversion used in sub model
  double t_steam = 373.15;
  double a = 1.0 - t_steam/temperature_K;
  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  double water_vp = 101325.0 * exp(a);          // (Pa)
  double ppm_to_RH = pressure_Pa / water_vp / 1.0e6; // (1/ppm)

  // Water activity and d_aw / d_[H2O_g]
  double a_w = ppm_to_RH * CONC_H2O_G;
  double d_aw_d_wg = ppm_to_RH;

  // Jacobson ion pair
  double Y0 = -1.918004e2;
  double Y1 =  2.001540e3;
  double Y2 = -8.557205e3;
  double Y3 =  1.987670e4;
  double Y4 = -2.717192e4;
  double Y5 =  2.187103e4;
  double Y6 = -9.591577e3;
  double Y7 =  1.763672e3;
  double molality = Y0 + Y1*a_w + Y2*pow(a_w,2) + Y3*pow(a_w,3) + Y4*pow(a_w,4) +
                    Y5*pow(a_w,5) + Y6*pow(a_w,6) + Y7*pow(a_w,7);
  double d_molal_d_wg = Y1 + 2.0*Y2*a_w + 3.0*Y3*pow(a_w,2) + 4.0*Y4*pow(a_w,3) +
                        5.0*Y5*pow(a_w,4) + 6.0*Y6*pow(a_w,5) + 7.0*Y7*pow(a_w,6);
  d_molal_d_wg *= d_aw_d_wg;

  double cation = CONC_CA / MW_CA / 1000.0;
  double d_cation_d_C = 1.0 / MW_CA / 1000.0;
  double anion = CONC_CL / 2.0 / MW_CL / 1000.0;
  double d_anion_d_A = 1.0 / 2.0 / MW_CL / 1000.0;

  if (cation>anion) {
    J_H2O += -2.0 * anion / pow(molality, 3) * 1000.0 * d_molal_d_wg;
    J_CL  +=  1.0 / pow(molality, 2) * 1000.0 * d_anion_d_A;
  } else {
    J_H2O += -2.0 * cation / pow(molality, 3) * 1000.0 * d_molal_d_wg;
    J_CA  +=  1.0 / pow(molality, 2) * 1000.0 * d_cation_d_C;
  }

  // EQSAM ion pair
  double NW = 2.0;
  double ZW = 0.67;
  double MW = 0.0585;

  molality = NW * 55.51 * 18.01 / MW / 1000.0 * (1.0/a_w-1.0);
  molality = pow(molality, ZW);
  d_molal_d_wg = NW * 55.01 * 18.01 / MW / 1000.0 / pow(a_w-1.0,2) * d_aw_d_wg;
  d_molal_d_wg = ZW * pow(molality, ZW-1.0) * d_molal_d_wg;

  J_H2O += -1.0 * CONC_CL / MW_CL / pow(molality,2) * d_molal_d_wg;
  J_CL  += 1.0 / MW_CL / molality;

  // Call the sub-model Jac calculation function
  sub_model_get_jac_contrib(md, J, 0.0);

  ret_val += ASSERT_MSG(SM_DATA_S(J)[0]==J_H2O, "gas-phase water Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[1]==J_H2O, "gas-phase water Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[2]==J_H2O, "gas-phase water Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[3]==J_CL,  "Cl- Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[5]==J_CL,  "Cl- Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[7]==J_CL,  "Cl- Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[4]==J_CA,  "Ca++ Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[6]==J_CA,  "Ca++ Jac element");
  ret_val += ASSERT_MSG(SM_DATA_S(J)[8]==J_CA,  "Ca++ Jac element");

  return ret_val;
}
#endif

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

  // Update the environmental state
  sub_model_update_env_state(md, env);

  // Run the calculation tests
  ret_val += test_sub_model_zsr_jac_calc(solver_data, state, env, J);

  // free allocated memory
  N_VDestroy(solver_state);
#endif

  return ret_val;
}

