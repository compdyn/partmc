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
#include "phlex_common.h"

/* Public aerosol phase functions*/

/* Solver functions */
void * aero_phase_get_mass(ModelData *model_data, int aero_phase_idx,
          double *state_var, double *mass, double *MW);
void * aero_phase_get_volume(ModelData *model_data, int aero_phase_idx,
          double *state_var, double *volume);
void * aero_phase_find(ModelData *model_data, int int_aero_phase_idx);
void * aero_phase_skip(void *aero_phase_data);
void aero_phase_print_data(void *solver_data);

/* Setup functions */
void aero_phase_add_condensed_data(int n_int_param, int n_float_param,
          int *int_param, double *float_param, void *solver_data);

#endif
