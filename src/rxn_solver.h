/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for reaction functions
 *
*/
/** \file
 * \brief Header file for solver functions
*/
#ifndef RXN_SOLVER_H_
#define RXN_SOLVER_H_
#include "phlex_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents used */

// arrhenius
void * rxn_arrhenius_get_used_jac_elem(void *rxn_data, bool **jac_struct);
void * rxn_arrhenius_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data);
void * rxn_arrhenius_update_env_state(realtype *env_data, void *rxn_data);
void * rxn_arrhenius_calc_deriv_contrib(realtype *state, realtype *deriv, void *rxn_data);
void * rxn_arrhenius_calc_jac_contrib(realtype *state, realtype *J, void *rxn_data);
void * rxn_arrhenius_skip(void *rxn_data);
void * rxn_arrhenius_print(void *rxn_data);

// CMAQ_H2O2
void * rxn_CMAQ_H2O2_get_used_jac_elem(void *rxn_data, bool **jac_struct);
void * rxn_CMAQ_H2O2_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data);
void * rxn_CMAQ_H2O2_update_env_state(realtype *env_data, void *rxn_data);
void * rxn_CMAQ_H2O2_calc_deriv_contrib(realtype *state, realtype *deriv, void *rxn_data);
void * rxn_CMAQ_H2O2_calc_jac_contrib(realtype *state, realtype *J, void *rxn_data);
void * rxn_CMAQ_H2O2_skip(void *rxn_data);
void * rxn_CMAQ_H2O2_print(void *rxn_data);

// CMAQ_OH_HNO3
void * rxn_CMAQ_OH_HNO3_get_used_jac_elem(void *rxn_data, bool **jac_struct);
void * rxn_CMAQ_OH_HNO3_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data);
void * rxn_CMAQ_OH_HNO3_update_env_state(realtype *env_data, void *rxn_data);
void * rxn_CMAQ_OH_HNO3_calc_deriv_contrib(realtype *state, realtype *deriv, void *rxn_data);
void * rxn_CMAQ_OH_HNO3_calc_jac_contrib(realtype *state, realtype *J, void *rxn_data);
void * rxn_CMAQ_OH_HNO3_skip(void *rxn_data);
void * rxn_CMAQ_OH_HNO3_print(void *rxn_data);

// photolysis
void * rxn_photolysis_get_used_jac_elem(void *rxn_data, bool **jac_struct);
void * rxn_photolysis_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data);
void * rxn_photolysis_update_env_state(realtype *env_data, void *rxn_data);
void * rxn_photolysis_calc_deriv_contrib(realtype *state, realtype *deriv, void *rxn_data);
void * rxn_photolysis_calc_jac_contrib(realtype *state, realtype *J, void *rxn_data);
void * rxn_photolysis_set_photo_rate(int photo_id, realtype base_rate, void *rxn_data);
void * rxn_photolysis_skip(void *rxn_data);
void * rxn_photolysis_print(void *rxn_data);

// troe
void * rxn_troe_get_used_jac_elem(void *rxn_data, bool **jac_struct);
void * rxn_troe_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data);
void * rxn_troe_update_env_state(realtype *env_data, void *rxn_data);
void * rxn_troe_calc_deriv_contrib(realtype *state, realtype *deriv, void *rxn_data);
void * rxn_troe_calc_jac_contrib(realtype *state, realtype *J, void *rxn_data);
void * rxn_troe_skip(void *rxn_data);
void * rxn_troe_print(void *rxn_data);


#endif
