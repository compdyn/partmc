/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for reaction functions
 *
 *
*/
/** \file
 * \brief Header file for reaction solver functions
*/
#ifndef RXNS_H_
#define RXNS_H_
#include "camp_gpu_solver.h"

//#define PMC_USE_SUNDIALS

// aqueous_equilibrium
void * rxn_gpu_aqueous_equilibrium_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_aqueous_equilibrium_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_aqueous_equilibrium_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_aqueous_equilibrium_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_aqueous_equilibrium_get_float_pointer(void *rxn_data);
void * rxn_gpu_aqueous_equilibrium_skip(
          void *rxn_data);
void * rxn_gpu_aqueous_equilibrium_print(
          void *rxn_data);

#ifdef PMC_USE_SUNDIALS
//__device__ double rxn_aqueous_equilibrium_calc_overall_rate(int *rxn_data,
//     double *double_pointer_gpu, double *rate_constants, double *state,
//     double react_fact, double prod_fact, double water, int i_phase, int n_rxn2)
void rxn_cpu_aqueous_equilibrium_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_aqueous_equilibrium_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// arrhenius
void * rxn_gpu_arrhenius_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_arrhenius_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_arrhenius_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_arrhenius_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_arrhenius_get_float_pointer(void *rxn_data);
void * rxn_gpu_arrhenius_skip(
          void *rxn_data);
void * rxn_gpu_arrhenius_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_arrhenius_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_arrhenius_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_arrhenius_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// CMAQ_H2O2
void * rxn_gpu_CMAQ_H2O2_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_CMAQ_H2O2_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_CMAQ_H2O2_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_CMAQ_H2O2_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_CMAQ_H2O2_get_float_pointer(void *rxn_data);
void * rxn_gpu_CMAQ_H2O2_skip(
          void *rxn_data);
void * rxn_gpu_CMAQ_H2O2_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_CMAQ_H2O2_calc_deriv_contrib(double *rate_constants, double *state, double *deriv,
          void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_CMAQ_H2O2_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// CMAQ_OH_HNO3
void * rxn_gpu_CMAQ_OH_HNO3_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_CMAQ_OH_HNO3_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_CMAQ_OH_HNO3_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_CMAQ_OH_HNO3_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_CMAQ_OH_HNO3_get_float_pointer(void *rxn_data);
void * rxn_gpu_CMAQ_OH_HNO3_skip(
          void *rxn_data);
void * rxn_gpu_CMAQ_OH_HNO3_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_CMAQ_OH_HNO3_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_CMAQ_OH_HNO3_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// condensed_phase_arrhenius
void * rxn_gpu_condensed_phase_arrhenius_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_condensed_phase_arrhenius_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_condensed_phase_arrhenius_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_condensed_phase_arrhenius_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_condensed_phase_arrhenius_get_float_pointer(void *rxn_data);
void * rxn_gpu_condensed_phase_arrhenius_skip(
          void *rxn_data);
void * rxn_gpu_condensed_phase_arrhenius_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_condensed_phase_arrhenius_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_condensed_phase_arrhenius_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// emission
void * rxn_gpu_emission_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_emission_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_emission_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_emission_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_emission_update_data(
          void *update_data, void *rxn_data);
void * rxn_gpu_emission_get_float_pointer(void *rxn_data);
void * rxn_gpu_emission_skip(
          void *rxn_data);
void * rxn_gpu_emission_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_emission_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_emission_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_emission_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif
void * rxn_gpu_emission_create_rate_update_data();
void rxn_gpu_emission_set_rate_update_data(
          void *update_data, int rxn_id, double base_rate);

// first_order_loss
void * rxn_gpu_first_order_loss_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_first_order_loss_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_first_order_loss_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_first_order_loss_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_first_order_loss_update_data(
          void *update_data, void *rxn_data);
void * rxn_gpu_first_order_loss_get_float_pointer(void *rxn_data);
void * rxn_gpu_first_order_loss_skip(
          void *rxn_data);
void * rxn_gpu_first_order_loss_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_first_order_loss_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_first_order_loss_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_first_order_loss_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif
void * rxn_gpu_first_order_loss_create_rate_update_data();
void rxn_gpu_first_order_loss_set_rate_update_data(
          void *update_data, int rxn_id, double base_rate);

// HL_phase_transfer
__device__ void rxn_gpu_HL_phase_transfer_update_env_state(double *rate_constants,
           int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_HL_phase_transfer_get_float_pointer(void *rxn_data);
void * rxn_gpu_HL_phase_transfer_skip(
          void *rxn_data);
void * rxn_gpu_HL_phase_transfer_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_HL_phase_transfer_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#ifndef FORCE_CPU
__device__
#endif
void rxn_gpu_HL_phase_transfer_calc_deriv_contrib(
          ModelData *model_data, double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
#ifndef FORCE_CPU
__device__
#endif
void rxn_gpu_HL_phase_transfer_calc_jac_contrib(
          ModelData *model_data, double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// photolysis
void * rxn_gpu_photolysis_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_photolysis_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_photolysis_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_photolysis_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_photolysis_update_data(
          void *update_data, void *rxn_data);
void * rxn_gpu_photolysis_get_float_pointer(void *rxn_data);
void * rxn_gpu_photolysis_skip(
          void *rxn_data);
void * rxn_gpu_photolysis_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_photolysis_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_photolysis_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_photolysis_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif
void * rxn_gpu_photolysis_create_rate_update_data();
void rxn_gpu_photolysis_set_rate_update_data(
          void *update_data, int photo_id, double base_rate);

// SIMPOL_phase_transfer
__device__ void rxn_gpu_SIMPOL_phase_transfer_update_env_state(double *rate_constants,
           int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_SIMPOL_phase_transfer_get_float_pointer(void *rxn_data);
void * rxn_gpu_SIMPOL_phase_transfer_skip(
          void *rxn_data);
void * rxn_gpu_SIMPOL_phase_transfer_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_SIMPOL_phase_transfer_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#ifndef FORCE_CPU
__device__
#endif
void rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
          ModelData *model_data, double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
#ifndef FORCE_CPU
__device__
#endif
void rxn_gpu_SIMPOL_phase_transfer_calc_jac_contrib(
          ModelData *model_data, double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// troe
void * rxn_gpu_troe_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_troe_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_troe_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_troe_pre_calc(ModelData *model_data, void *rxn_data);
void * rxn_gpu_troe_get_float_pointer(void *rxn_data);
void * rxn_gpu_troe_skip(
          void *rxn_data);
void * rxn_gpu_troe_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_troe_calc_deriv_contrib(double *rate_constants, double *state, double *deriv,
          void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_troe_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_troe_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif

// wet_deposition
void * rxn_gpu_wet_deposition_get_used_jac_elem(
          void *rxn_data, bool **jac_struct);
void * rxn_gpu_wet_deposition_update_ids(
          ModelData *model_data, int *deriv_ids, int **jac_ids, void *rxn_data);
__device__ void rxn_gpu_wet_deposition_update_env_state(double *rate_constants,
          int n_rxn2, double *double_pointer_gpu, double *env_data, void *rxn_data);
void * rxn_gpu_wet_deposition_pre_calc(
          ModelData *model_data, void *rxn_data);
void * rxn_gpu_wet_deposition_update_data(
          void *update_data, void *rxn_data);
void * rxn_gpu_wet_deposition_get_float_pointer(void *rxn_data);
void * rxn_gpu_wet_deposition_skip(
          void *rxn_data);
void * rxn_gpu_wet_deposition_print(
          void *rxn_data);
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_wet_deposition_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_wet_deposition_calc_deriv_contrib(
          double *rate_constants, double *state, realtype *deriv, void *rxn_data,
          double * double_pointer_gpu, double time_step, int n_rxn);
__device__ void rxn_gpu_wet_deposition_calc_jac_contrib(
          double *rate_constants, double *state, realtype *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn);
#endif
void * rxn_gpu_wet_deposition_create_rate_update_data();
void rxn_gpu_wet_deposition_set_rate_update_data(
          void *update_data, int rxn_id, double base_rate);


#endif
