/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Aerosol representation-specific functions for use by the solver
 *
 */
/** \file
 * \brief Aerosol representation functions
 */
extern "C" {
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aero_rep_solver_gpu.h"
#include "aero_reps_gpu.h"

// Aerosol representations (Must match parameters defined in pmc_aero_rep_gpu_factory
#define AERO_REP_SINGLE_PARTICLE   1
#define AERO_REP_MODAL_BINNED_MASS 2


/** \brief Get the effective particle radius, \f$r_{eff}\f$ (m)
 *
 * Calculates effective particle radius \f$r_{eff}\f$ (m), as well as the set of
 * \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_gpu_idx Index of aerosol representation to use for calculation
 * \param aero_phase_gpu_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param radius Pointer to hold effective particle radius (m)
 * \return A pointer to a set of partial derivatives
 *         \f$\frac{\partial r_{eff}}{\partial y}\f$, or a NULL pointer if no
 *         partial derivatives exist
 */
__device__ void *aero_rep_gpu_get_effective_radius(ModelDatagpu *model_data, int aero_rep_gpu_idx,
                                    int aero_phase_gpu_idx, double *radius) {

  // Set up a pointer for the partial derivatives
  double *partial_deriv = NULL;

  // Get the number of aerosol representations
  int *aero_rep_gpu_data = (int *) (model_data->aero_rep_gpu_data);
  int n_aero_rep = *(aero_rep_gpu_data++);

  // Loop through the aerosol representations to find the one requested
  for (int i_aero_rep = 0; i_aero_rep < aero_rep_gpu_idx; i_aero_rep++) {

    // Get the aerosol representation type
    int aero_rep_gpu_type = *(aero_rep_gpu_data++);

    // Advance the pointer to the next aerosol representation
    switch (aero_rep_gpu_type) {
        case AERO_REP_MODAL_BINNED_MASS : {
          aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_skip(
                  (void *) aero_rep_gpu_data);
        break;
      }
      case AERO_REP_SINGLE_PARTICLE : {
        aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_skip(
                (void *) aero_rep_gpu_data);
        break;
      }
    }
  }

  // Get the aerosol representation type
  int aero_rep_gpu_type = *(aero_rep_gpu_data++);

  // Get the particle radius and set of partial derivatives
  switch (aero_rep_gpu_type) {
    case AERO_REP_MODAL_BINNED_MASS : {
      aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_get_effective_radius(
              aero_phase_gpu_idx, radius, partial_deriv, (void *) aero_rep_gpu_data);
      break;
    }
    case AERO_REP_SINGLE_PARTICLE : {
      aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_get_effective_radius(
              aero_phase_gpu_idx, radius, partial_deriv, (void *) aero_rep_gpu_data);
      break;
    }
  }
  return partial_deriv;
}

/** \brief Get the particle number concentration \f$n\f$ (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * Calculates particle number concentration, \f$n\f$
 * (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$), as well as the set of
 * \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_gpu_idx Index of aerosol representation to use for calculation
 * \param aero_phase_gpu_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param number_conc Pointer to hold calculated number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \return A pointer to a set of partial derivatives
 *         \f$\frac{\partial n}{\partial y}\f$, or a NULL pointer if no partial
 *         derivatives exist
 */
__device__ void *aero_rep_gpu_get_number_conc(ModelDatagpu *model_data, int aero_rep_gpu_idx,
                               int aero_phase_gpu_idx, double *number_conc) {

  // Set up a pointer for the partial derivatives
  double *partial_deriv = NULL;

  // Get the number of aerosol representations
  int *aero_rep_gpu_data = (int *) (model_data->aero_rep_gpu_data);
  int n_aero_rep = *(aero_rep_gpu_data++);

  // Loop through the aerosol representations to find the one requested
  for (int i_aero_rep = 0; i_aero_rep < aero_rep_gpu_idx; i_aero_rep++) {

    // Get the aerosol representation type
    int aero_rep_gpu_type = *(aero_rep_gpu_data++);

    // Advance the pointer to the next aerosol representation
    switch (aero_rep_gpu_type) {
      case AERO_REP_MODAL_BINNED_MASS :
        aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_skip(
                (void *) aero_rep_gpu_data);
        break;
      case AERO_REP_SINGLE_PARTICLE :
        aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_skip(
                (void *) aero_rep_gpu_data);
        break;
    }
  }

  // Get the aerosol representation type
  int aero_rep_gpu_type = *(aero_rep_gpu_data++);

  // Get the particle number concentration
  switch (aero_rep_gpu_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_get_number_conc(
              aero_phase_gpu_idx, number_conc, partial_deriv,
              (void *) aero_rep_gpu_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_get_number_conc(
              aero_phase_gpu_idx, number_conc, partial_deriv,
              (void *) aero_rep_gpu_data);
      break;
  }
  return partial_deriv;
}

/** \brief Check whether aerosol concentrations are per-particle or total for each phase
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_gpu_idx Index of aerosol representation to use for calculation
 * \param aero_phase_gpu_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \return 0 for per-particle; 1 for total for each phase
 */
__device__ int aero_rep_gpu_get_aero_conc_type(ModelDatagpu *model_data, int aero_rep_gpu_idx,
                                int aero_phase_gpu_idx) {

  // Initialize the aerosol concentration type
  int aero_conc_type = 0;

  // Get the number of aerosol representations
  int *aero_rep_gpu_data = (int *) (model_data->aero_rep_gpu_data);
  int n_aero_rep = *(aero_rep_gpu_data++);

  // Loop through the aerosol representations to find the one requested
  for (int i_aero_rep = 0; i_aero_rep < aero_rep_gpu_idx; i_aero_rep++) {

    // Get the aerosol representation type
    int aero_rep_gpu_type = *(aero_rep_gpu_data++);

    // Advance the pointer to the next aerosol representation
    switch (aero_rep_gpu_type) {
      case AERO_REP_MODAL_BINNED_MASS :
        aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_skip(
                (void *) aero_rep_gpu_data);
        break;
      case AERO_REP_SINGLE_PARTICLE :
        aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_skip(
                (void *) aero_rep_gpu_data);
        break;
    }
  }

  // Get the aerosol representation type
  int aero_rep_gpu_type = *(aero_rep_gpu_data++);

  // Get the type of aerosol concentration
  switch (aero_rep_gpu_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_gpu_data = (int *) aero_rep_gpu_modal_binned_mass_get_aero_conc_type(
              aero_phase_gpu_idx, &aero_conc_type, (void *) aero_rep_gpu_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_gpu_data = (int *) aero_rep_gpu_single_particle_get_aero_conc_type(
              aero_phase_gpu_idx, &aero_conc_type, (void *) aero_rep_gpu_data);
      break;
  }
  return aero_conc_type;
}


/** \brief Get the total mass of an aerosol phase in this representation \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * Calculates total aerosol phase mass, \f$m\f$
 * (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$), as well as the set of
 * \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are variables on the
 * solver state array.
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_idx Index of aerosol representation to use for calculation
 * \param aero_phase_idx Index of the aerosol phase within the aerosol
 *                       representation
 * \param aero_phase_mass Pointer to hold calculated aerosol-phase mass,
 *                        \f$m\f$
 *                        (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 * \param aero_phase_avg_MW Pointer to hold calculated average MW in the
 *                          aerosol phase (\f$\mbox{\si{\kilogram\per\mole}}\f$)
 * \return A pointer to a set of partial derivatives
 *         \f$\frac{\partial m}{\partial y}\f$, or a NULL pointer if no partial
 *         derivatives exist
 */
__device__ void * aero_rep_gpu_get_aero_phase_mass(ModelDatagpu *model_data, int aero_rep_idx,
                                    int aero_phase_idx, double *aero_phase_mass,
                                    double *aero_phase_avg_MW)
{

  // Set up a pointer for the partial derivatives
  double *partial_deriv = NULL;

  // Get the number of aerosol representations
  int *aero_rep_data = (int*) (model_data->aero_rep_gpu_data);
  int n_aero_rep = *(aero_rep_data++);

  // Loop through the aerosol representations to find the one requested
  for (int i_aero_rep=0; i_aero_rep<aero_rep_idx; i_aero_rep++) {

    // Get the aerosol representation type
    int aero_rep_type = *(aero_rep_data++);

    // Advance the pointer to the next aerosol representation
    switch (aero_rep_type) {
      case AERO_REP_MODAL_BINNED_MASS :
        aero_rep_data = (int*) aero_rep_gpu_modal_binned_mass_skip(
                (void*) aero_rep_data);
        break;
      case AERO_REP_SINGLE_PARTICLE :
        aero_rep_data = (int*) aero_rep_gpu_single_particle_skip(
                (void*) aero_rep_data);
        break;
    }
  }

  // Get the aerosol representation type
  int aero_rep_type = *(aero_rep_data++);

  // Get the particle number concentration
  switch (aero_rep_type) {
    case AERO_REP_MODAL_BINNED_MASS :
      aero_rep_data = (int*) aero_rep_gpu_modal_binned_mass_get_aero_phase_mass(
              aero_phase_idx, aero_phase_mass, aero_phase_avg_MW,
              partial_deriv, (void*) aero_rep_data);
      break;
    case AERO_REP_SINGLE_PARTICLE :
      aero_rep_data = (int*) aero_rep_gpu_single_particle_get_aero_phase_mass(
              aero_phase_idx, aero_phase_mass, aero_phase_avg_MW,
              partial_deriv, (void*) aero_rep_data);
      break;
  }
  return partial_deriv;
}


}

