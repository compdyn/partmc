/* Copyright (C) 2017-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Single particle aerosol representation functions
 *
 */
/** \file
 * \brief Single particle aerosol representation functions
 */
#include <stdio.h>
#include <stdlib.h>
#include "../aero_phase_solver.h"
#include "../aero_reps.h"
#include "../phlex_solver.h"

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
#define NUM_ENV_PARAM_ 0
#define PHASE_STATE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PHASE_MODEL_DATA_ID_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PHASE_NUM_JAC_ELEM_(x) int_data[NUM_INT_PROP_+2*NUM_PHASE_+x]
#define PHASE_MASS_(x) (float_data[NUM_FLOAT_PROP_+x])
#define PHASE_AVG_MW_(x) (float_data[NUM_FLOAT_PROP_+NUM_PHASE_+x])

/** \brief Flag Jacobian elements used in calcualtions of mass and volume
 *
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param model_data Pointer to the model data
 * \param aero_phase_idx Index of the aerosol phase to find elements for
 * \param jac_struct 1D array of flags indicating potentially non-zero
 *                   Jacobian elements. (The dependent variable should have
 *                   been chosen by the calling function.)
 * \return Number of Jacobian elements flagged
 */
int aero_rep_single_particle_get_used_jac_elem(ModelData *model_data,
          int aero_phase_idx, int *aero_rep_int_data,
          double *aero_rep_float_data, bool *jac_struct)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  int n_jac_elem = 0;

  // Each phase in a single particle has the same jac elements
  // (one for each species in each phase in the particle)
  for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
    PHASE_NUM_JAC_ELEM_(i_phase) =
      aero_phase_get_used_jac_elem( model_data,
              PHASE_MODEL_DATA_ID_(i_phase),
              PHASE_STATE_ID_(i_phase), jac_struct );
    n_jac_elem += PHASE_NUM_JAC_ELEM_(i_phase);
  }

  return n_jac_elem;
}

/** \brief Flag elements on the state array used by this aerosol representation
 *
 * The single particle aerosol representation functions do not use state array
 * values
 *
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param state_flags Array of flags indicating state array elements used
 */
void aero_rep_single_particle_get_dependencies(int *aero_rep_int_data,
          double *aero_rep_float_data, bool *state_flags)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  return;
}

/** \brief Update aerosol representation data for new environmental conditions
 *
 * The single particle aerosol representation does not use environmental
 * conditions
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_update_env_state(ModelData *model_data,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;
  double *env_data = model_data->grid_cell_env;

  return;
}

/** \brief Update aerosol representation data for a new state
 *
 * Updates the mass (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$) and
 * average MW (\f$\mbox{\si{\kilogram\per\mole}}\f$) for each aerosol phase in
 * the particle
 *
 * \param model_data Pointer to the model data, include the state array
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_update_state(ModelData *model_data,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  // Calculate the total aerosol phase masses
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {

    // Get a pointer to the phase on the state array
    double *state_var = (double*) (model_data->grid_cell_state);
    state_var += PHASE_STATE_ID_(i_phase);

    // Get the mass and average MW
    aero_phase_get_mass(model_data, PHASE_MODEL_DATA_ID_(i_phase), state_var,
               &(PHASE_MASS_(i_phase)), &(PHASE_AVG_MW_(i_phase)), NULL, NULL);
  }

  return;
}

/** \brief Get the effective particle radius \f$r_{eff}\f$ (m)
 *
 * The single particle radius is set by the aerosol micro-physics model prior to
 * solving the chemistry. Thus, all \f$\frac{\partial r_{eff}}{\partial y}\f$
 * are zero. Also, there is only one set of particles in the single particle
 * representation, so the phase index is not used.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 *                       (not used)
 * \param radius Effective particle radius (m)
 * \param partial_deriv \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$
 *                      are species on the state array
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_get_effective_radius(ModelData *model_data,
          int aero_phase_idx, double *radius, double *partial_deriv,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  *radius = RADIUS_;

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
  }

  return;
}

/** \brief Get the particle number concentration \f$n\f$ (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * This single particle number concentration is set by the aerosol model prior
 * to solving the chemistry. Thus, all \f$\frac{\partial n}{\partial y}\f$ are
 * zero. Also, there is only one set of particles in the single particle
 * representation, so the phase index is not used.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 *                       (not used)
 * \param number_conc Particle number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \param partial_deriv \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_get_number_conc(ModelData *model_data,
          int aero_phase_idx, double *number_conc, double *partial_deriv,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  *number_conc = NUMBER_CONC_;

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
  }

  return;
}

/** \brief Get the type of aerosol concentration used.
 *
 * Single particle concentrations are per-particle.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_conc_type Pointer to int that will hold the concentration type
 *                       code
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_get_aero_conc_type(int aero_phase_idx,
          int *aero_conc_type, int *aero_rep_int_data,
          double *aero_rep_float_data, double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  *aero_conc_type = 0;

  return;
}

/** \brief Get the total mass in an aerosol phase \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * The single particle mass is set for each new state as the sum of the masses
 * of the aerosol phases that compose the particle
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_phase_mass Total mass in the aerosol phase, \f$m\f$
 *                        (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 * \param partial_deriv \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_get_aero_phase_mass(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_mass, double *partial_deriv,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  *aero_phase_mass = PHASE_MASS_(aero_phase_idx);

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
      if (i_phase == aero_phase_idx) {
        double *state = (double*) (model_data->grid_cell_state);
        state += PHASE_STATE_ID_(i_phase);
        double mass, mw;
        aero_phase_get_mass(model_data, aero_phase_idx, state, &mass, &mw,
                            partial_deriv, NULL);
        partial_deriv += PHASE_NUM_JAC_ELEM_(i_phase);
      } else {
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
      }
    }
  }

  return;
}

/** \brief Get the average molecular weight in an aerosol phase
 **        \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
 *
 * The single particle mass is set for each new state as the sum of the masses
 * of the aerosol phases that compose the particle
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_phase_avg_MW Average molecular weight in the aerosol phase
 *                          (\f$\mbox{\si{\kilogram\per\mole}}\f$)
 * \param partial_deriv \f$\frac{\partial m}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_get_aero_phase_avg_MW(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_avg_MW, double *partial_deriv,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  *aero_phase_avg_MW = PHASE_AVG_MW_(aero_phase_idx);

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
      if (i_phase == aero_phase_idx) {
        double *state = (double*) (model_data->grid_cell_state);
        state += PHASE_STATE_ID_(i_phase);
        double mass, mw;
        aero_phase_get_mass(model_data, aero_phase_idx, state, &mass, &mw,
                            NULL, partial_deriv);
        partial_deriv += PHASE_NUM_JAC_ELEM_(i_phase);
      } else {
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
      }
    }
  }

  return;
}

/** \brief Update aerosol representation data
 *
 * Single particle aerosol representation update data is structured as follows:
 *
 *  - \b int aero_rep_id (Id of one or more aerosol representations set by the
 *       host model using the
 *       pmc_aero_rep_single_particle::aero_rep_single_particle_t::set_id
 *       function prior to initializing the solver.)
 *  - \b int update_type (Type of update to perform. Can be UPDATE_RADIUS or
 *       UPDATE_NUMBER.)
 *  - \b double new_value (Either the new radius (m) or the new number
 *       concentration (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$).)
 *
 * \param update_data Pointer to the updated aerosol representation data
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 * \param aero_rep_env_data Pointer to the aerosol representation
 *                          environment-dependent parameters
 */
void aero_rep_single_particle_update_data(void *update_data,
          int *aero_rep_int_data, double *aero_rep_float_data,
          double *aero_rep_env_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  int *aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(aero_rep_id[1]);
  double *new_value = (double*) &(update_type[1]);

  // Set the new radius or number concentration for matching aerosol
  // representations
  if (*aero_rep_id==AERO_REP_ID_ && AERO_REP_ID_!=0) {
    if (*update_type==UPDATE_RADIUS) {
      RADIUS_ = (double) *new_value;
    } else if (*update_type==UPDATE_NUMBER) {
      NUMBER_CONC_ = (double) *new_value;
    }
  }

  return;
}

/** \brief Print the Single Particle reaction parameters
 *
 * \param aero_rep_int_data Pointer to the aerosol representation integer data
 * \param aero_rep_float_data Pointer to the aerosol representation
 *                            floating-point data
 */
void aero_rep_single_particle_print(int *aero_rep_int_data,
          double *aero_rep_float_data)
{
  int *int_data = aero_rep_int_data;
  double *float_data = aero_rep_float_data;

  printf("\n\nSingle particle aerosol representation\n");

  return;
}

/** \brief Create update data for new particle radius
 *
 * \return Pointer to a new radius update data object
 */
void * aero_rep_single_particle_create_radius_update_data()
{
  int *update_data = (int*) malloc(2*sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for radius update data\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set radius update data
 *
 * \param update_data Pointer to an allocated radius update data object
 * \param aero_rep_id Id of the aerosol representation(s) to update
 * \param radius New particle radius
 */
void aero_rep_single_particle_set_radius_update_data(void *update_data,
          int aero_rep_id, double radius)
{
  int *new_aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(new_aero_rep_id[1]);
  double *new_radius = (double*) &(update_type[1]);
  *new_aero_rep_id = aero_rep_id;
  *update_type = UPDATE_RADIUS;
  *new_radius = radius;
}

/** \brief Create update data for new particle number
 *
 * \return Pointer to a new number update data object
 */
void * aero_rep_single_particle_create_number_update_data()
{
  int *update_data = (int*) malloc(2*sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for number update data\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set number update data
 *
 * \param update_data Pointer to an allocated number update data object
 * \param aero_rep_id Id of the aerosol representation(s) to update
 * \param number_conc New particle number
 */
void aero_rep_single_particle_set_number_update_data(void *update_data,
          int aero_rep_id, double number_conc)
{
  int *new_aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(new_aero_rep_id[1]);
  double *new_number_conc = (double*) &(update_type[1]);
  *new_aero_rep_id = aero_rep_id;
  *update_type = UPDATE_NUMBER;
  *new_number_conc = number_conc;
}
