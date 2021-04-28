/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for reaction functions
 *
 * TODO Automatically generate rxn_solver.c and rxn_solver.h code
 * maybe using cmake?
 *
 */
/** \file
 * \brief Header file for reaction solver functions
 */
#ifndef RXNS_H_
#define RXNS_H_
#include "Jacobian.h"
#include "camp_common.h"

// aqueous_equilibrium
void rxn_aqueous_equilibrium_get_used_jac_elem(int *rxn_int_data,
                                               double *rxn_float_data,
                                               Jacobian *jac);
void rxn_aqueous_equilibrium_update_ids(ModelData *model_data, int *deriv_ids,
                                        Jacobian jac, int *rxn_int_data,
                                        double *rxn_float_data);
void rxn_aqueous_equilibrium_update_env_state(ModelData *model_data,
                                              int *rxn_int_data,
                                              double *rxn_float_data,
                                              double *rxn_env_data);
void rxn_aqueous_equilibrium_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_aqueous_equilibrium_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_aqueous_equilibrium_calc_jac_contrib(ModelData *model_data,
                                              Jacobian jac, int *rxn_int_data,
                                              double *rxn_float_data,
                                              double *rxn_env_data,
                                              realtype time_step);
#endif

// arrhenius
void rxn_arrhenius_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
                                     Jacobian *jac);
void rxn_arrhenius_update_ids(ModelData *model_data, int *deriv_ids,
                              Jacobian jac, int *rxn_int_data,
                              double *rxn_float_data);
void rxn_arrhenius_update_env_state(ModelData *model_data, int *rxn_int_data,
                                    double *rxn_float_data,
                                    double *rxn_env_data);
void rxn_arrhenius_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_arrhenius_calc_deriv_contrib(ModelData *model_data,
                                      TimeDerivative time_deriv,
                                      int *rxn_int_data, double *rxn_float_data,
                                      double *rxn_env_data, realtype time_step);
void rxn_arrhenius_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                    int *rxn_int_data, double *rxn_float_data,
                                    double *rxn_env_data, realtype time_step);
#endif

// CMAQ_H2O2
void rxn_CMAQ_H2O2_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
                                     Jacobian *jac);
void rxn_CMAQ_H2O2_update_ids(ModelData *model_data, int *deriv_ids,
                              Jacobian jac, int *rxn_int_data,
                              double *rxn_float_data);
void rxn_CMAQ_H2O2_update_env_state(ModelData *model_data, int *rxn_int_data,
                                    double *rxn_float_data,
                                    double *rxn_env_data);
void rxn_CMAQ_H2O2_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_CMAQ_H2O2_calc_deriv_contrib(ModelData *model_data,
                                      TimeDerivative time_deriv,
                                      int *rxn_int_data, double *rxn_float_data,
                                      double *rxn_env_data, realtype time_step);
void rxn_CMAQ_H2O2_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                    int *rxn_int_data, double *rxn_float_data,
                                    double *rxn_env_data, realtype time_step);
#endif

// CMAQ_OH_HNO3
void rxn_CMAQ_OH_HNO3_get_used_jac_elem(int *rxn_int_data,
                                        double *rxn_float_data, Jacobian *jac);
void rxn_CMAQ_OH_HNO3_update_ids(ModelData *model_data, int *deriv_ids,
                                 Jacobian jac, int *rxn_int_data,
                                 double *rxn_float_data);
void rxn_CMAQ_OH_HNO3_update_env_state(ModelData *model_data, int *rxn_int_data,
                                       double *rxn_float_data,
                                       double *rxn_env_data);
void rxn_CMAQ_OH_HNO3_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_CMAQ_OH_HNO3_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_CMAQ_OH_HNO3_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                       int *rxn_int_data,
                                       double *rxn_float_data,
                                       double *rxn_env_data,
                                       realtype time_step);
#endif

// condensed_phase_arrhenius
void rxn_condensed_phase_arrhenius_get_used_jac_elem(int *rxn_int_data,
                                                     double *rxn_float_data,
                                                     Jacobian *jac);
void rxn_condensed_phase_arrhenius_update_ids(ModelData *model_data,
                                              int *deriv_ids, Jacobian jac,
                                              int *rxn_int_data,
                                              double *rxn_float_data);
void rxn_condensed_phase_arrhenius_update_env_state(ModelData *model_data,
                                                    int *rxn_int_data,
                                                    double *rxn_float_data,
                                                    double *rxn_env_data);
void rxn_condensed_phase_arrhenius_print(int *rxn_int_data,
                                         double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_condensed_phase_arrhenius_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_condensed_phase_arrhenius_calc_jac_contrib(
    ModelData *model_data, Jacobian jac, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
#endif

// emission
void rxn_emission_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
                                    Jacobian *jac);
void rxn_emission_update_ids(ModelData *model_data, int *deriv_ids,
                             Jacobian jac, int *rxn_int_data,
                             double *rxn_float_data);
void rxn_emission_update_env_state(ModelData *model_data, int *rxn_int_data,
                                   double *rxn_float_data,
                                   double *rxn_env_data);
bool rxn_emission_update_data(void *update_data, int *rxn_int_data,
                              double *rxn_float_data, double *rxn_env_data);
void rxn_emission_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_emission_calc_deriv_contrib(ModelData *model_data,
                                     TimeDerivative time_deriv,
                                     int *rxn_int_data, double *rxn_float_data,
                                     double *rxn_env_data, realtype time_step);
void rxn_emission_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                   int *rxn_int_data, double *rxn_float_data,
                                   double *rxn_env_data, realtype time_step);
#endif
void *rxn_emission_create_rate_update_data();
void rxn_emission_set_rate_update_data(void *update_data, int rxn_id,
                                       double base_rate);

// first_order_loss
void rxn_first_order_loss_get_used_jac_elem(int *rxn_int_data,
                                            double *rxn_float_data,
                                            Jacobian *jac);
void rxn_first_order_loss_update_ids(ModelData *model_data, int *deriv_ids,
                                     Jacobian jac, int *rxn_int_data,
                                     double *rxn_float_data);
void rxn_first_order_loss_update_env_state(ModelData *model_data,
                                           int *rxn_int_data,
                                           double *rxn_float_data,
                                           double *rxn_env_data);
bool rxn_first_order_loss_update_data(void *update_data, int *rxn_int_data,
                                      double *rxn_float_data,
                                      double *rxn_env_data);
void rxn_first_order_loss_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_first_order_loss_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_first_order_loss_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                           int *rxn_int_data,
                                           double *rxn_float_data,
                                           double *rxn_env_data,
                                           realtype time_step);
#endif
void *rxn_first_order_loss_create_rate_update_data();
void rxn_first_order_loss_set_rate_update_data(void *update_data, int rxn_id,
                                               double base_rate);

// HL_phase_transfer
void rxn_HL_phase_transfer_get_used_jac_elem(ModelData *model_data,
                                             int *rxn_int_data,
                                             double *rxn_float_data,
                                             Jacobian *jac);
void rxn_HL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
                                      Jacobian jac, int *rxn_int_data,
                                      double *rxn_float_data);
void rxn_HL_phase_transfer_update_env_state(ModelData *model_data,
                                            int *rxn_int_data,
                                            double *rxn_float_data,
                                            double *rxn_env_data);
void rxn_HL_phase_transfer_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_HL_phase_transfer_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_HL_phase_transfer_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                            int *rxn_int_data,
                                            double *rxn_float_data,
                                            double *rxn_env_data,
                                            realtype time_step);
#endif

// photolysis
void rxn_photolysis_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
                                      Jacobian *jac);
void rxn_photolysis_update_ids(ModelData *model_data, int *deriv_ids,
                               Jacobian jac, int *rxn_int_data,
                               double *rxn_float_data);
void rxn_photolysis_update_env_state(ModelData *model_data, int *rxn_int_data,
                                     double *rxn_float_data,
                                     double *rxn_env_data);
bool rxn_photolysis_update_data(void *update_data, int *rxn_int_data,
                                double *rxn_float_data, double *rxn_env_data);
void rxn_photolysis_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_photolysis_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_photolysis_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                     int *rxn_int_data, double *rxn_float_data,
                                     double *rxn_env_data, realtype time_step);
#endif
void *rxn_photolysis_create_rate_update_data();
void rxn_photolysis_set_rate_update_data(void *update_data, int photo_id,
                                         double base_rate);

// SIMPOL_phase_transfer
void rxn_SIMPOL_phase_transfer_get_used_jac_elem(ModelData *model_data,
                                                 int *rxn_int_data,
                                                 double *rxn_float_data,
                                                 Jacobian *jac);
void rxn_SIMPOL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
                                          Jacobian jac, int *rxn_int_data,
                                          double *rxn_float_data);
void rxn_SIMPOL_phase_transfer_update_env_state(ModelData *model_data,
                                                int *rxn_int_data,
                                                double *rxn_float_data,
                                                double *rxn_env_data);
void rxn_SIMPOL_phase_transfer_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_SIMPOL_phase_transfer_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_SIMPOL_phase_transfer_calc_jac_contrib(ModelData *model_data,
                                                Jacobian jac, int *rxn_int_data,
                                                double *rxn_float_data,
                                                double *rxn_env_data,
                                                realtype time_step);
#endif

// ternary_chemical_activation
void rxn_ternary_chemical_activation_get_used_jac_elem(int *rxn_int_data,
                                                       double *rxn_float_data,
                                                       Jacobian *jac);
void rxn_ternary_chemical_activation_update_ids(ModelData *model_data,
                                                int *deriv_ids, Jacobian jac,
                                                int *rxn_int_data,
                                                double *rxn_float_data);
void rxn_ternary_chemical_activation_update_env_state(ModelData *model_data,
                                                      int *rxn_int_data,
                                                      double *rxn_float_data,
                                                      double *rxn_env_data);
void rxn_ternary_chemical_activation_print(int *rxn_int_data,
                                           double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_ternary_chemical_activation_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_ternary_chemical_activation_calc_jac_contrib(
    ModelData *model_data, Jacobian jac, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
#endif

// troe
void rxn_troe_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
                                Jacobian *jac);
void rxn_troe_update_ids(ModelData *model_data, int *deriv_ids, Jacobian jac,
                         int *rxn_int_data, double *rxn_float_data);
void rxn_troe_update_env_state(ModelData *model_data, int *rxn_int_data,
                               double *rxn_float_data, double *rxn_env_data);
void rxn_troe_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_troe_calc_deriv_contrib(ModelData *model_data,
                                 TimeDerivative time_deriv, int *rxn_int_data,
                                 double *rxn_float_data, double *rxn_env_data,
                                 realtype time_step);
void rxn_troe_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                               int *rxn_int_data, double *rxn_float_data,
                               double *rxn_env_data, realtype time_step);
#endif

// wennberg_tunneling
void rxn_wennberg_tunneling_get_used_jac_elem(int *rxn_int_data,
                                              double *rxn_float_data,
                                              Jacobian *jac);
void rxn_wennberg_tunneling_update_ids(ModelData *model_data, int *deriv_ids,
                                       Jacobian jac, int *rxn_int_data,
                                       double *rxn_float_data);
void rxn_wennberg_tunneling_update_env_state(ModelData *model_data,
                                             int *rxn_int_data,
                                             double *rxn_float_data,
                                             double *rxn_env_data);
bool rxn_wennberg_tunneling_update_data(void *update_data, int *rxn_int_data,
                                        double *rxn_float_data,
                                        double *rxn_env_data);
void rxn_wennberg_tunneling_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_wennberg_tunneling_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_wennberg_tunneling_calc_jac_contrib(ModelData *model_data,
                                             Jacobian jac, int *rxn_int_data,
                                             double *rxn_float_data,
                                             double *rxn_env_data,
                                             realtype time_step);
#endif
void *rxn_wennberg_tunneling_create_rate_update_data();
void rxn_wennberg_tunneling_set_rate_update_data(void *update_data, int rxn_id,
                                                 double base_rate);

// wet_deposition
void rxn_wet_deposition_get_used_jac_elem(int *rxn_int_data,
                                          double *rxn_float_data,
                                          Jacobian *jac);
void rxn_wet_deposition_update_ids(ModelData *model_data, int *deriv_ids,
                                   Jacobian jac, int *rxn_int_data,
                                   double *rxn_float_data);
void rxn_wet_deposition_update_env_state(ModelData *model_data,
                                         int *rxn_int_data,
                                         double *rxn_float_data,
                                         double *rxn_env_data);
bool rxn_wet_deposition_update_data(void *update_data, int *rxn_int_data,
                                    double *rxn_float_data,
                                    double *rxn_env_data);
void rxn_wet_deposition_print(int *rxn_int_data, double *rxn_float_data);
#ifdef PMC_USE_SUNDIALS
void rxn_wet_deposition_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step);
void rxn_wet_deposition_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                         int *rxn_int_data,
                                         double *rxn_float_data,
                                         double *rxn_env_data,
                                         realtype time_step);
#endif
void *rxn_wet_deposition_create_rate_update_data();
void rxn_wet_deposition_set_rate_update_data(void *update_data, int rxn_id,
                                             double base_rate);

#endif
