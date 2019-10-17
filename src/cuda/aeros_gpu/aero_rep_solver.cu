/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aerosol representation-specific functions for use by the solver
 *
 */
/** \file
 * \brief Aerosol representation functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aero_rep_solver_gpu.h"
#include "aero_reps_gpu.h"

// Aerosol representations (Must match parameters defined in pmc_aero_rep_factory
#define AERO_REP_SINGLE_PARTICLE   1
#define AERO_REP_MODAL_BINNED_MASS 2

extern "C" {

  /*

#ifndef FORCE_CPU
__device__
#endif
void aero_rep_gpu_get_effective_radius(ModelData *model_data, int aero_rep_idx,
                                   int aero_phase_idx, double *radius, double *partial_deriv)
{

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
          model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle radius and set of partial derivatives
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_modal_binned_mass_get_effective_radius(
              model_data, aero_phase_idx, radius, partial_deriv,
              aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_single_particle_get_effective_radius(
              model_data, aero_phase_idx, radius, partial_deriv,
              aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return;
}

#ifndef FORCE_CPU
__device__
#endif
void aero_rep_gpu_get_number_conc(ModelData *model_data, int aero_rep_idx,
                              int aero_phase_idx, double *number_conc, double *partial_deriv)
{

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
          model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_modal_binned_mass_get_number_conc(
              model_data, aero_phase_idx, number_conc, partial_deriv,
              aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_single_particle_get_number_conc(
              model_data, aero_phase_idx, number_conc, partial_deriv,
              aero_rep_int_data, aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return;
}


#ifndef FORCE_CPU
__device__
#endif
int aero_rep_gpu_get_aero_conc_type(ModelData *model_data, int aero_rep_idx,
                                int aero_phase_idx)
{

  // Initialize the aerosol concentration type
  int aero_conc_type = 0;

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
          model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the type of aerosol concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_modal_binned_mass_get_aero_conc_type(
              aero_phase_idx, &aero_conc_type, aero_rep_int_data,
              aero_rep_float_data, aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_single_particle_get_aero_conc_type(
              aero_phase_idx, &aero_conc_type, aero_rep_int_data,
              aero_rep_float_data, aero_rep_env_data);
      break;
  }
  return aero_conc_type;
}


#ifndef FORCE_CPU
__device__
#endif
void aero_rep_gpu_get_aero_phase_mass(ModelData *model_data, int aero_rep_idx,
                                  int aero_phase_idx, double *aero_phase_mass, double *partial_deriv)
{

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
          model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_modal_binned_mass_get_aero_phase_mass(
              model_data, aero_phase_idx, aero_phase_mass,
              partial_deriv, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_single_particle_get_aero_phase_mass(
              model_data, aero_phase_idx, aero_phase_mass,
              partial_deriv, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
      break;
  }
}


#ifndef FORCE_CPU
__device__
#endif
void aero_rep_gpu_get_aero_phase_avg_MW(ModelData *model_data, int aero_rep_idx,
                                    int aero_phase_idx, double *aero_phase_avg_MW, double *partial_deriv)
{

  // Get pointers to the aerosol data
  int *aero_rep_int_data =
          &(model_data->aero_rep_int_data[model_data->aero_rep_int_indices[aero_rep_idx]]);
  double *aero_rep_float_data =
          &(model_data->aero_rep_float_data[model_data->aero_rep_float_indices[aero_rep_idx]]);
  double *aero_rep_env_data   = &(model_data->grid_cell_aero_rep_env_data[
          model_data->aero_rep_env_idx[aero_rep_idx]]);

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_int_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_modal_binned_mass_get_aero_phase_avg_MW(
              model_data, aero_phase_idx, aero_phase_avg_MW,
              partial_deriv, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_single_particle_get_aero_phase_avg_MW(
              model_data, aero_phase_idx, aero_phase_avg_MW,
              partial_deriv, aero_rep_int_data, aero_rep_float_data,
              aero_rep_env_data);
      break;
  }
}

   */

}

