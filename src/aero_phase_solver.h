/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for aero_phase_solver.c
 *
 */
/** \file
 * \brief Header file for aerosol phase functions
 */
#ifndef AERO_PHASE_SOLVER_H
#define AERO_PHASE_SOLVER_H
#include "camp_common.h"

/* Public aerosol phase functions*/

/* Solver functions */
int aero_phase_get_used_jac_elem(ModelData *model_data, int aero_phase_idx,
                                 int state_var_id, bool *jac_struct);
void aero_phase_get_mass(ModelData *model_data, int aero_phase_idx,
                         double *state_var, double *mass, double *MW,
                         double *jac_elem_mass, double *jac_elem_MW);
void aero_phase_get_volume(ModelData *model_data, int aero_phase_idx,
                           double *state_var, double *volume, double *jac_elem);
void *aero_phase_find(ModelData *model_data, int int_aero_phase_idx);
void aero_phase_print_data(void *solver_data);

/* Setup functions */
void aero_phase_add_condensed_data(int n_int_param, int n_float_param,
                                   int *int_param, double *float_param,
                                   void *solver_data);

#endif
