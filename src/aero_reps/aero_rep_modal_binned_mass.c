/* Copyright (C) 2017-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Modal mass aerosol representation functions
 *
 */
/** \file
 * \brief Modal mass aerosol representation functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aero_phase_solver.h"
#include "../aero_reps.h"
#include "../phlex_solver.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define UPDATE_GMD 0
#define UPDATE_GSD 1

#define BINNED 1
#define MODAL 2

#define NUM_SECTION_ (int_data[0])
#define INT_DATA_SIZE_ (int_data[1])
#define FLOAT_DATA_SIZE_ (int_data[2])
#define AERO_REP_ID_ (int_data[3])
#define NUM_INT_PROP_ 4
#define NUM_FLOAT_PROP_ 0
#define MODE_INT_PROP_LOC_(x) (int_data[NUM_INT_PROP_+x]-1)
#define MODE_FLOAT_PROP_LOC_(x) (int_data[NUM_INT_PROP_+NUM_SECTION_+x]-1)
#define SECTION_TYPE_(x) (int_data[MODE_INT_PROP_LOC_(x)])

// For modes, NUM_BINS_ = 1
#define NUM_BINS_(x) (int_data[MODE_INT_PROP_LOC_(x)+1])

// Number of aerosol phases in this mode/bin set
#define NUM_PHASE_(x) (int_data[MODE_INT_PROP_LOC_(x)+2])

// Phase state and model data ids
#define PHASE_STATE_ID_(x,y,b) (int_data[MODE_INT_PROP_LOC_(x)+3+b*NUM_PHASE_(x)+y]-1)
#define PHASE_MODEL_DATA_ID_(x,y,b) (int_data[MODE_INT_PROP_LOC_(x)+3+NUM_BINS_(x)*NUM_PHASE_(x)+b*NUM_PHASE_(x)+y]-1)

// Number of Jacobian elements in a phase
#define PHASE_NUM_JAC_ELEM_(x,y,b) int_data[MODE_INT_PROP_LOC_(x)+3+2*NUM_BINS_(x)*NUM_PHASE_(x)+b*NUM_PHASE_(x)+y]

// GMD and bin diameter are stored in the same position - for modes, b=0
#define GMD_(x,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4])
#define BIN_DP_(x,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4])

// GSD - only used for modes, b=0
#define GSD_(x,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+1])

// Real-time number concentration - used for modes and bins - for modes, b=0
#define NUMBER_CONC_(x,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+2])

// Real-time effective radius - only used for modes, b=0
#define EFFECTIVE_RADIUS_(x,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+3])

// Real-time phase mass (ug/m^3) - used for modes and bins - for modes, b=0
#define PHASE_MASS_(x,y,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+4*NUM_BINS_(x)+b*NUM_PHASE_(x)+y])

// Real-time phase average MW (kg/mol) - used for modes and bins - for modes, b=0
#define PHASE_AVG_MW_(x,y,b) (float_data[MODE_FLOAT_PROP_LOC_(x)+(4+NUM_PHASE_(x))*NUM_BINS_(x)+b*NUM_PHASE_(x)+y])

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
int aero_rep_modal_binned_mass_get_used_jac_elem(ModelData *model_data,
          int aero_phase_idx, void *aero_rep_data, bool *jac_struct)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  int num_flagged_elem = 0;

  // Loop through the modes/bins flagging Jacobian elements used by each
  // aerosol phase
  for (int i_section=0; i_section<NUM_SECTION_ && aero_phase_idx>=0;
            i_section++) {
    for (int i_bin=0; i_bin<NUM_BINS_(i_section) && aero_phase_idx>=0;
              i_bin++) {
      for (int i_phase=0; i_phase<NUM_PHASE_(i_section) && aero_phase_idx>=0;
                i_phase++) {
        if (aero_phase_idx==0) {
          for (int j_phase=0; j_phase<NUM_PHASE_(i_section); j_phase++) {
            PHASE_NUM_JAC_ELEM_(i_section, j_phase, i_bin) =
              aero_phase_get_used_jac_elem( model_data,
                  PHASE_MODEL_DATA_ID_(i_section, j_phase, i_bin),
                  PHASE_STATE_ID_(i_section, j_phase, i_bin), jac_struct );
            num_flagged_elem += PHASE_NUM_JAC_ELEM_(i_section, j_phase, i_bin);
          }
        }
        aero_phase_idx-=1;
      }
    }
  }

  return num_flagged_elem;
}

/** \brief Flag elements on the state array used by this aerosol representation
 *
 * The modal mass aerosol representation functions do not use state array values
 *
 * \param aero_rep_data A pointer to the aerosol representation data
 * \param state_flags Array of flags indicating state array elements used
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
void * aero_rep_modal_binned_mass_get_dependencies(void *aero_rep_data,
          bool *state_flags)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update aerosol representation data for new environmental conditions
 *
 * The modal mass aerosol representation is not updated for new environmental
 * conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_modal_binned_mass_update_env_state(double *env_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update aerosol representation data for a new state
 *
 * The modal mass aerosol representation recalculates effective radius and
 * number concentration for each new state.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_modal_binned_mass_update_state(ModelData *model_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Loop through the modes and calculate effective radius and number
  // concentration
  for (int i_section=0; i_section<NUM_SECTION_; i_section++) {

    double volume, mass;
    switch (SECTION_TYPE_(i_section)) {

      // Mode
      case (MODAL) :

        // Sum the volumes of each species in the mode
        volume = 0.0;
        for (int i_phase=0; i_phase<NUM_PHASE_(i_section); i_phase++) {

          // Get a pointer to the phase on the state array
          double *state = (double*) (model_data->state);
          state += PHASE_STATE_ID_(i_section, i_phase, 0);

          // Set the aerosol-phase mass and average MW
          aero_phase_get_mass(model_data,
                    PHASE_MODEL_DATA_ID_(i_section, i_phase, 0),
                    state, &(PHASE_MASS_(i_section, i_phase, 0)),
                    &(PHASE_AVG_MW_(i_section, i_phase, 0)), NULL, NULL);

          // Get the phase volume
          double phase_volume = 0.0;
          aero_phase_get_volume(model_data,
                    PHASE_MODEL_DATA_ID_(i_section, i_phase, 0),
                    state, &phase_volume, NULL);
          volume += phase_volume;

        }

        // Calculate the number concentration based on the total mode volume
        // (see aero_rep_modal_binned_mass_get_number_conc for details)
        NUMBER_CONC_(i_section,0) = volume * 6.0 /
                (M_PI * pow(GMD_(i_section,0),3) *
                exp(9.0/2.0*pow(GSD_(i_section,0),2)));

        break;

      // Bins
      case (BINNED) :

        // Loop through the bins
        for (int i_bin=0; i_bin<NUM_BINS_(i_section); i_bin++) {

          // Sum the volumes of each species in the bin
          volume = 0.0;
          for (int i_phase=0; i_phase<NUM_PHASE_(i_section); i_phase++) {

            // Get a pointer to the phase on the state array
            double *state = (double*) (model_data->state);
            state += PHASE_STATE_ID_(i_section, i_phase, i_bin);

            // Set the aerosol-phase mass and average MW
            aero_phase_get_mass(model_data,
                      PHASE_MODEL_DATA_ID_(i_section, i_phase, i_bin),
                      state, &(PHASE_MASS_(i_section, i_phase, i_bin)),
                      &(PHASE_AVG_MW_(i_section, i_phase, i_bin)), NULL, NULL);

            // Get the phase volume
            double phase_volume = 0.0;
            aero_phase_get_volume(model_data,
                      PHASE_MODEL_DATA_ID_(i_section, i_phase, i_bin),
                      state, &phase_volume, NULL);
            volume += phase_volume;

          }

          // Calculate the number concentration based on the total bin volume
          // (see aero_rep_modal_binned_mass_get_number_conc for details)
          NUMBER_CONC_(i_section, i_bin) = volume * 3.0 / (4.0*M_PI) /
                  pow(BIN_DP_(i_section, i_bin)/2.0, 3);
        }

        break;
    }

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the effective particle radius \f$r_{eff}\f$ (m)
 *
 * The modal mass effective radius is calculated for a log-normal distribution
 * where the geometric mean diameter (\f$\tilde{D}_n\f$) and geometric standard
 * deviation (\f$\tilde{\sigma}_g\f$) are set by the aerosol model prior to
 * solving the chemistry. Thus, all \f$\frac{\partial r_{eff}}{\partial y}\f$
 * are zero. The effective radius is calculated according to the equation given
 * in Table 1 of Zender \cite Zender2002 :
 *
 * \f[
 *      r_{eff} = \frac{\tilde{D}_n}{2}*exp(9\tilde{\sigma}_g^2/2)
 * \f]
 * \f[
 *      r_{eff} = \frac{D_{eff}}{2}
 * \f]
 *
 * For bins, \f$r_{eff}\f$ is assumed to be the bin radius.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param radius Effective particle radius (m)
 * \param partial_deriv \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$
 *                       are species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_modal_binned_mass_get_effective_radius(ModelData *model_data,
          int aero_phase_idx, double *radius, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_section=0; i_section<NUM_SECTION_; i_section++) {
    for (int i_bin=0; i_bin<NUM_BINS_(i_section); i_bin++) {
      aero_phase_idx-=NUM_PHASE_(i_section);
      if (aero_phase_idx<0) {
        *radius = EFFECTIVE_RADIUS_(i_section, i_bin);
        // Effective radii are constant for bins and modes
        if (partial_deriv) {
          for (int i_phase = 0; i_phase < NUM_PHASE_(i_section); ++i_phase) {
            for (int i_elem = 0;
                 i_elem < PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
                 ++i_elem) {
              *(partial_deriv++) = ZERO;
            }
          }
        }
        i_section = NUM_SECTION_;
        break;
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the particle number concentration \f$n\f$ (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * The modal mass number concentration is calculated for a log-normal
 * distribution where the geometric mean diameter (\f$\tilde{D}_n\f$) and
 * geometric standard deviation (\f$\tilde{\sigma}_g\f$) are set by the aerosol
 * model prior to solving the chemistry. The number concentration is
 * calculated according to the equation given in Table 1 of Zender
 * \cite Zender2002 :
 * \f[
 *      n = N_0 = \frac{6V_0}{\pi}\tilde{D}_n^{-3}e^{-9\tilde{\sigma}_g^2/2}
 * \f]
 * \f[
 *      V_0 = \sum_i{\frac{m_i}{\rho_i}}
 * \f]
 * where \f$\rho_i\f$ and \f$m_i\f$ are the density and total mass of species
 * \f$i\f$ in the specified mode.
 *
 * The binned number concentration is calculated according to:
 * \f[
 *     n = V_0 / V_p
 * \f]
 * \f[
 *     V_p = \frac{4}{3}\pi r^{3}
 * \f]
 * where \f$r\f$ is the radius of the size bin and \f$V_0\f$ is defined as
 * above.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param number_conc Particle number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \param partial_deriv \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_modal_binned_mass_get_number_conc(ModelData *model_data,
          int aero_phase_idx, double *number_conc, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_section=0; i_section<NUM_SECTION_ && aero_phase_idx>=0;
       i_section++) {
    for (int i_bin=0; i_bin<NUM_BINS_(i_section) && aero_phase_idx>=0;
         i_bin++) {
      aero_phase_idx-=NUM_PHASE_(i_section);
      if (aero_phase_idx<0) {
        *number_conc = NUMBER_CONC_(i_section, i_bin);
        if (partial_deriv) {
          for (int i_phase = 0; i_phase < NUM_PHASE_(i_section); ++i_phase) {

            // Get a pointer to the phase on the state array
            double *state = (double*) (model_data->state);
            state += PHASE_STATE_ID_(i_section, i_phase, i_bin);

            // Get the aerosol phase volume
            double phase_volume = 0.0;
            aero_phase_get_volume(model_data,
                      PHASE_MODEL_DATA_ID_(i_section, i_phase, i_bin),
                      state, &phase_volume, partial_deriv);

            // Convert d_vol/d_conc to d_number/d_conc
            for (int i_elem = 0;
                 i_elem < PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
                 ++i_elem) {
              switch (SECTION_TYPE_(i_section)) {
                case (MODAL) :
                  *(partial_deriv++) *= 6.0 / (M_PI * pow(GMD_(i_section,0),3) *
                                        exp(9.0/2.0*pow(GSD_(i_section,0),2)));
                  break;
                case (BINNED) :
                  *(partial_deriv++) *= 3.0 / (4.0*M_PI) /
                                        pow(BIN_DP_(i_section, i_bin)/2.0, 3);
                  break;
              }
            }
          }
        }
        i_section = NUM_SECTION_;
        break;
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the type of aerosol concentration used.
 *
 * Modal mass concentrations are per-mode or per-bin.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_conc_type Pointer to int that will hold the concentration type
 *                       code
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
void * aero_rep_modal_binned_mass_get_aero_conc_type(int aero_phase_idx,
          int *aero_conc_type, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  *aero_conc_type = 1;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the total mass in an aerosol phase \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
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
void * aero_rep_modal_binned_mass_get_aero_phase_mass(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_mass, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_section=0; i_section<NUM_SECTION_ && aero_phase_idx>=0;
       ++i_section) {
    for (int i_bin=0; i_bin<NUM_BINS_(i_section) && aero_phase_idx>=0;
         ++i_bin) {
      if (aero_phase_idx<0 || aero_phase_idx>=NUM_PHASE_(i_section)) {
        aero_phase_idx -= NUM_PHASE_(i_section);
        continue;
      }
      for (int i_phase=0; i_phase<NUM_PHASE_(i_section); ++i_phase) {
        if (aero_phase_idx==0) {
          *aero_phase_mass = PHASE_MASS_(i_section, i_phase, i_bin);
          if (partial_deriv) {

            // Get a pointer to the phase on the state array
            double *state = (double*) (model_data->state);
            state += PHASE_STATE_ID_(i_section, i_phase, i_bin);

            // Get d_mass / d_conc
            double mass, mw;
            aero_phase_get_mass(model_data,
                      PHASE_MODEL_DATA_ID_(i_section, i_phase, i_bin),
                      state, &mass, &mw, partial_deriv, NULL);
            partial_deriv += PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
          }

        // Other phases present in the bin or mode do not contribute to
        // the aerosol phase mass
        } else if (partial_deriv) {
          for( int i_elem=0;
               i_elem < PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
               ++i_elem) {
            *(partial_deriv++) = ZERO;
          }
        }
        aero_phase_idx-=1;
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the average molecular weight in an aerosol phase
 **        \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
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
void * aero_rep_modal_binned_mass_get_aero_phase_avg_MW(ModelData *model_data,
          int aero_phase_idx, double *aero_phase_avg_MW, double *partial_deriv,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_section=0; i_section<NUM_SECTION_ && aero_phase_idx>=0;
       ++i_section) {
    for (int i_bin=0; i_bin<NUM_BINS_(i_section) && aero_phase_idx>=0;
         ++i_bin) {
      if (aero_phase_idx<0 || aero_phase_idx>=NUM_PHASE_(i_section)) {
        aero_phase_idx -= NUM_PHASE_(i_section);
        continue;
      }
      for (int i_phase=0; i_phase<NUM_PHASE_(i_section); ++i_phase) {
        if (aero_phase_idx==0) {
          *aero_phase_avg_MW = PHASE_AVG_MW_(i_section, i_phase, i_bin);
          if (partial_deriv) {

            // Get a pointer to the phase on the state array
            double *state = (double*) (model_data->state);
            state += PHASE_STATE_ID_(i_section, i_phase, i_bin);

            // Get d_MW / d_conc
            double mass, mw;
            aero_phase_get_mass(model_data,
                      PHASE_MODEL_DATA_ID_(i_section, i_phase, i_bin),
                      state, &mass, &mw, NULL, partial_deriv);
            partial_deriv += PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
          }

        // Other phases present in the bin/mode do not contribute to the
        // average MW of the aerosol phase
        } else if (partial_deriv) {
          for( int i_elem=0;
               i_elem < PHASE_NUM_JAC_ELEM_(i_section, i_phase, i_bin);
               ++i_elem) {
            *(partial_deriv++) = ZERO;
          }
        }
        aero_phase_idx-=1;
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update the aerosol representation data
 *
 * The model mass aerosol representation update data is structured as follows:
 *
 *  - \b int aero_rep_id (Id of one or more aerosol representations set by the
 *       host model using the
 *       pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::set_id
 *       function prior to initializing the solver.)
 *  - \b int update_type (Type of update to perform. Can be UPDATE_GMD or
 *       UPDATE_GSD.)
 *  - \b int section_id (Index of the mode to update.)
 *  - \b double new_value (Either the new GMD (m) or the new GSD (unitless).)
 *
 * \param update_data Pointer to the updated aerosol representation data
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *        representation data
 */
void * aero_rep_modal_binned_mass_update_data(void *update_data,
          void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  int *aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(aero_rep_id[1]);
  int *section_id = (int*) &(update_type[1]);
  double *new_value = (double*) &(section_id[1]);

  // Set the new GMD or GSD for matching aerosol representations
  if (*aero_rep_id==AERO_REP_ID_ && AERO_REP_ID_!=0) {
    if (*update_type==UPDATE_GMD) {
      GMD_(*section_id,0) = (double) *new_value;
    } else if (*update_type==UPDATE_GSD) {
      GSD_(*section_id,0) = (double) *new_value;
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the mass-only modal/binned reaction parameters
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
void * aero_rep_modal_binned_mass_print(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nModal/binned mass-only aerosol representation\n");
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
void * aero_rep_modal_binned_mass_skip(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Create update data for new GMD
 *
 * \return Pointer to a new GMD update data object
 */
void * aero_rep_modal_binned_mass_create_gmd_update_data()
{
  int *update_data = (int*) malloc(3*sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for GMD update data.\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set GMD update data
 *
 * \param update_data Pointer to an allocated GMD update data object
 * \param aero_rep_id Id of the aerosol representation(s) to update
 * \param section_id Id of the mode to update
 * \param gmd New mode GMD (m)
 */
void aero_rep_modal_binned_mass_set_gmd_update_data(void *update_data,
          int aero_rep_id, int section_id, double gmd)
{
  int *new_aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(new_aero_rep_id[1]);
  int *new_section_id = (int*) &(update_type[1]);
  double *new_GMD = (double*) &(new_section_id[1]);
  *new_aero_rep_id = aero_rep_id;
  *update_type = UPDATE_GMD;
  *new_section_id = section_id;
  *new_GMD = gmd;
}

/** \brief Create update data for new GSD
 *
 * \return Pointer to a new GSD update data object
 */
void * aero_rep_modal_binned_mass_create_gsd_update_data()
{
  int *update_data = (int*) malloc(3*sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for GSD update data.\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set GSD update data
 *
 * \param update_data Pointer to an allocated GSD update data object
 * \param aero_rep_id Id of the aerosol representation(s) to update
 * \param section_id Id of the mode to update
 * \param gsd New mode GSD (unitless)
 */
void aero_rep_modal_binned_mass_set_gsd_update_data(void *update_data,
          int aero_rep_id, int section_id, double gsd)
{
  int *new_aero_rep_id = (int*) update_data;
  int *update_type = (int*) &(new_aero_rep_id[1]);
  int *new_section_id = (int*) &(update_type[1]);
  double *new_GSD = (double*) &(new_section_id[1]);
  *new_aero_rep_id = aero_rep_id;
  *update_type = UPDATE_GSD;
  *new_section_id = section_id;
  *new_GSD = gsd;
}

#undef BINNED
#undef MODAL

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef UPDATE_GSD
#undef UPDATE_GMD

#undef NUM_SECTION_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
#undef AERO_REP_ID_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef MODE_INT_PROP_LOC_
#undef MODE_FLOAT_PROP_LOC_
#undef SECTION_TYPE_
#undef NUM_BINS_
#undef NUM_PHASE_
#undef PHASE_STATE_ID_
#undef PHASE_MODEL_DATA_ID_
#undef GMD_
#undef BIN_DP_
#undef GSD_
#undef NUMBER_CONC_
#undef EFFECTIVE_RADIUS_
#undef PHASE_MASS_
#undef PHASE_AVG_MW_
