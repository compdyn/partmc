/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for aero_rep_solver.c
 *
 */
/** \file
 * \brief Header file for abstract aerosol representation functions
 */
#ifndef AERO_REP_SOLVER_H
#define AERO_REP_SOLVER_H
#include "camp_common.h"

/** Public aerosol representation functions **/

/* Solver functions */
int aero_rep_get_used_jac_elem(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx,
    bool *jac_struct);
void aero_rep_get_dependencies(
    ModelData *model_data, bool *state_flags);
void aero_rep_update_env_state(
    ModelData *model_data);
void aero_rep_update_state(
    ModelData *model_data);
void aero_rep_get_effective_radius(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx, double *radius,
    double *partial_deriv);
void aero_rep_get_number_conc(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx,
    double *number_conc, double *partial_deriv);
int aero_rep_get_aero_conc_type(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx);
void aero_rep_get_aero_phase_mass(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx,
    double *aero_phase_mass, double *partial_deriv);
void aero_rep_get_aero_phase_avg_MW(
    ModelData *model_data, int aero_rep_idx, int aero_phase_idx,
    double *aero_phase_avg_MW, double *partial_deriv);
void aero_rep_print_data(
    void *solver_data);

/* Setup functions */
void aero_rep_add_condensed_data(
    int aero_rep_type, int n_int_param, int n_float_param, int n_env_param,
    int *int_param, double *float_param, void *solver_data);

/* Update data functions */
void aero_rep_update_data(
    int cell_id, int *aero_rep_id, int update_aero_rep_type, void *update_data,
    void *solver_data);
void aero_rep_free_update_data(
    void *update_data);

#endif
