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
#define PHASE_STATE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PHASE_MODEL_DATA_ID_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PHASE_NUM_JAC_ELEM_(x) int_data[NUM_INT_PROP_+2*NUM_PHASE_+x]
#define PHASE_MASS_(x) (float_data[NUM_FLOAT_PROP_+x])
#define PHASE_AVG_MW_(x) (float_data[NUM_FLOAT_PROP_+NUM_PHASE_+x])
#define INT_DATA_SIZE_ (NUM_INT_PROP_+3*NUM_PHASE_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+2*NUM_PHASE_)

/** \brief Flag Jacobian elements used in calcualtions of mass and volume
 *
 * \param model_data Pointer to the model data
 * \param aero_rep_data A pointer to the aerosol representation data
 * \param aero_phase_idx Index of the aerosol phase to find elements for
 * \param jac_struct 1D array of flags indicating potentially non-zero
 *                   Jacobian elements. (The dependent variable should have
 *                   been chosen by the calling function.)
 * \return Number of Jacobian elements flagged
 */
int aero_rep_single_particle_get_used_jac_elem(ModelData *model_data,
          int aero_phase_idx, void *aero_rep_data, bool *jac_struct)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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
 * \param aero_rep_data A pointer to the aerosol representation data
 * \param state_flags Array of flags indicating state array elements used
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
void * aero_rep_single_particle_get_dependencies(void *aero_rep_data,
          bool *state_flags)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update aerosol representation data for new environmental conditions
 *
 * The single particle aerosol representation does not use environmental
 * conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_update_env_state(double *env_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update aerosol representation data for a new state
 *
 * Updates the mass (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$) and
 * average MW (\f$\mbox{\si{\kilogram\per\mole}}\f$) for each aerosol phase in
 * the particle
 *
 * \param model_data Pointer to the model data, include the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_update_state(ModelData *model_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the total aerosol phase masses
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {

    // Get a pointer to the phase on the state array
    double *state_var = (double*) (model_data->state);
    state_var += PHASE_STATE_ID_(i_phase);

    // Get the mass and average MW
    aero_phase_get_mass(model_data, PHASE_MODEL_DATA_ID_(i_phase), state_var,
               &(PHASE_MASS_(i_phase)), &(PHASE_AVG_MW_(i_phase)), NULL, NULL);
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_get_effective_radius(ModelData *model_data,
          int aero_phase_idx, double *radius, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *radius = RADIUS_;

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_get_number_conc(ModelData *model_data,
          int aero_phase_idx, double *number_conc, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *number_conc = NUMBER_CONC_;

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase)
      for (int i_spec = 0; i_spec < PHASE_NUM_JAC_ELEM_(i_phase); ++i_spec)
        *(partial_deriv++) = ZERO;
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
void * aero_rep_single_particle_get_aero_conc_type(int aero_phase_idx,
          int *aero_conc_type, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *aero_conc_type = 0;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_get_aero_phase_mass(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_mass, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *aero_phase_mass = PHASE_MASS_(aero_phase_idx);

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
      if (i_phase == aero_phase_idx) {
        double *state = (double*) (model_data->state);
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

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_single_particle_get_aero_phase_avg_MW(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_avg_MW, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *aero_phase_avg_MW = PHASE_AVG_MW_(aero_phase_idx);

  if (partial_deriv) {
    for (int i_phase = 0; i_phase < NUM_PHASE_; ++i_phase) {
      if (i_phase == aero_phase_idx) {
        double *state = (double*) (model_data->state);
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

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *        representaiton data
 */
void * aero_rep_single_particle_update_data(void *update_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Single Particle reaction parameters
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
void * aero_rep_single_particle_print(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nSingle particle aerosol representation\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Advance the aerosol representation data pointer to the next aerosol representation
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
void * aero_rep_single_particle_skip(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
