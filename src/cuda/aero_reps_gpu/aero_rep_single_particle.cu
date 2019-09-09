/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Single particle aerosol representation functions
 *
 */
/** \file
 * \brief Single particle aerosol representation functions
 */
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_phase_solver_gpu.h"
#include "../aeros_gpu/aero_reps_gpu.h"
#include "../camp_gpu_solver.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define UPDATE_RADIUS 0
#define UPDATE_NUMBER 1

#define NUM_PHASE_ int_data[0]
#define AERO_REP_ID_ int_data[1]
#define RADIUS_ float_data[0]
#define NUMBER_CONC_ float_data[1]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 2
#define PHASE_STATE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PHASE_MODEL_DATA_ID_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PHASE_MASS_(x) (float_data[NUM_FLOAT_PROP_+x])
#define PHASE_AVG_MW_(x) (float_data[NUM_FLOAT_PROP_+NUM_PHASE_+x])
#define INT_DATA_SIZE_ (NUM_INT_PROP_+2*NUM_PHASE_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+2*NUM_PHASE_)


/** \brief Get the effective particle radius \f$r_{eff}\f$ (m)
 *
 * The single particle radius is set by the aerosol micro-physics model prior to
 * solving the chemistry. Thus, all \f$\frac{\partial r_{eff}}{\partial y}\f$
 * are zero. Also, there is only one set of particles in the single particle
 * representation, so the phase index is not used.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 *                       (not used)
 * \param radius Effective particle radius (m)
 * \param partial_deriv \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$
 *                      are species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void *aero_rep_gpu_single_particle_get_effective_radius(int aero_phase_idx,
                                                    double *radius, double *partial_deriv, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  *radius = RADIUS_;

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the particle number concentration \f$n\f$ (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * This single particle number concentration is set by the aerosol model prior
 * to solving the chemistry. Thus, all \f$\frac{\partial n}{\partial y}\f$ are
 * zero. Also, there is only one set of particles in the single particle
 * representation, so the phase index is not used.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 *                       (not used)
 * \param number_conc Particle number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \param partial_deriv \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void *aero_rep_gpu_single_particle_get_number_conc(int aero_phase_idx,
                                               double *number_conc, double *partial_deriv, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  *number_conc = NUMBER_CONC_;

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the type of aerosol concentration used.
 *
 * Single particle concentrations are per-particle.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_conc_type Pointer to int that will hold the concentration type
 *                       code
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void *aero_rep_gpu_single_particle_get_aero_conc_type(int aero_phase_idx,
                                                  int *aero_conc_type, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  *aero_conc_type = 0;

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the total mass in an aerosol phase \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * The single particle mass is set for each new state as the sum of the masses
 * of the aerosol phases that compose the particle
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_phase_mass Total mass in the aerosol phase, \f$m\f$
 *                        (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 * \param aero_phase_avg_MW Average molecular weight in the aerosol phase
 *                          (\f$\mbox{\si{\kilogram\per\mole}}\f$)
 * \param partial_deriv \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void * aero_rep_gpu_single_particle_get_aero_phase_mass(int aero_phase_idx,
                                                    double *aero_phase_mass, double *aero_phase_avg_MW,
                                                    double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *aero_phase_mass = PHASE_MASS_(aero_phase_idx);
  *aero_phase_avg_MW = PHASE_AVG_MW_(aero_phase_idx);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Advance the aerosol representation data pointer to the next aerosol representation
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
__device__ void *aero_rep_gpu_single_particle_skip(void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef UPDATE_RADIUS
#undef UPDATE_NUMBER

#undef NUM_PHASE_
#undef AERO_REP_ID_
#undef RADIUS_
#undef NUMBER_CONC_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef PHASE_STATE_ID_
#undef PHASE_MODEL_DATA_ID_
#undef PHASE_MASS_
#undef PHASE_AVG_MW_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_

}
