/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for sub model calculations
 */
/** \file
 * \brief Header file for sub model functions
 */
#ifndef SUB_MODELS_H
#define SUB_MODELS_H
#include "phlex_common.h"

// PD-FiTE activity
void sub_model_PDFiTE_get_used_jac_elem(
          int *sub_model_int_data, double *sub_model_float_data,
          bool **jac_struct);
void sub_model_PDFiTE_update_ids(
          int *sub_model_int_data, double *sub_model_float_data,
          int *deriv_ids, int **jac_ids);
void sub_model_PDFiTE_update_env_state(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
void sub_model_PDFiTE_calculate(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
#ifdef PMC_USE_SUNDIALS
void sub_model_PDFiTE_get_jac_contrib(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data, realtype *J,
          double time_step);
#endif
void sub_model_PDFiTE_print(
          int *sub_model_int_data, double *sub_model_float_data);

// UNIFAC
void sub_model_UNIFAC_get_used_jac_elem(
          int *sub_model_int_data, double *sub_model_float_data,
          bool **jac_struct);
void sub_model_UNIFAC_update_ids(
          int *sub_model_int_data, double *sub_model_float_data,
          int *deriv_ids, int **jac_ids);
void sub_model_UNIFAC_update_env_state(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
void sub_model_UNIFAC_calculate(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
#ifdef PMC_USE_SUNDIALS
void sub_model_UNIFAC_get_jac_contrib(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data, realtype *J,
          double time_step);
#endif
void sub_model_UNIFAC_print(
          int *sub_model_int_data, double *sub_model_float_data);

// ZSR_aerosol_water
void sub_model_ZSR_aerosol_water_get_used_jac_elem(
          int *sub_model_int_data, double *sub_model_float_data,
          bool **jac_struct);
void sub_model_ZSR_aerosol_water_update_ids(
          int *sub_model_int_data, double *sub_model_float_data,
          int *deriv_ids, int **jac_ids);
void sub_model_ZSR_aerosol_water_update_env_state(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
void sub_model_ZSR_aerosol_water_calculate(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data);
#ifdef PMC_USE_SUNDIALS
void sub_model_ZSR_aerosol_water_get_jac_contrib(
          int *sub_model_int_data, double *sub_model_float_data,
          double *sub_model_env_data, ModelData *model_data, realtype *J,
          double time_step);
#endif
void sub_model_ZSR_aerosol_water_print(
          int *sub_model_int_data, double *sub_model_float_data);

#endif
