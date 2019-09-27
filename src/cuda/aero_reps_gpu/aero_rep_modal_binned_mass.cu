/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Modal mass aerosol representation functions
 *
 */
/** \file
 * \brief Modal mass aerosol representation functions
 */
extern "C" {
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_phase_solver_gpu.h"
#include "../aeros_gpu/aero_reps_gpu.h"
#include "../camp_gpu_solver.h"

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
#define PHASE_STATE_ID_(x, y, b) (int_data[MODE_INT_PROP_LOC_(x)+3+b*NUM_PHASE_(x)+y]-1)
#define PHASE_MODEL_DATA_ID_(x, y, b) (int_data[MODE_INT_PROP_LOC_(x)+3+NUM_BINS_(x)*NUM_PHASE_(x)+b*NUM_PHASE_(x)+y]-1)

// GMD and bin diameter are stored in the same position - for modes, b=0
#define GMD_(x, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4])
#define BIN_DP_(x, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4])

// GSD - only used for modes, b=0
#define GSD_(x, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+1])

// Real-time number concentration - used for modes and bins - for modes, b=0
#define NUMBER_CONC_(x, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+2])

// Real-time effective radius - only used for modes, b=0
#define EFFECTIVE_RADIUS_(x, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+b*4+3])

// Real-time phase mass (ug/m^3) - used for modes and bins - for modes, b=0
#define PHASE_MASS_(x, y, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+4*NUM_BINS_(x)+b*NUM_PHASE_(x)+y])

// Real-time phase average MW (kg/mol) - used for modes and bins - for modes, b=0
#define PHASE_AVG_MW_(x, y, b) (float_data[MODE_FLOAT_PROP_LOC_(x)+(4+NUM_PHASE_(x))*NUM_BINS_(x)+b*NUM_PHASE_(x)+y])


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
 *      r_{eff} = \frac{\tilde{D}_n}{2}*exp(5\tilde{\sigma}_g^2/2)
 * \f]
 * \f[
 *      r_{eff} = \frac{D_{eff}}{2}
 * \f]
 *
 * For bins, \f$r_{eff}\f$ is assumed to be the bin radius.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param radius Effective particle radius (m)
 * \param partial_deriv \f$\frac{\partial r_{eff}}{\partial y}\f$ where \f$y\f$
 *                       are species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void *aero_rep_gpu_modal_binned_mass_get_effective_radius(int aero_phase_idx,
                                                      double *radius, double *partial_deriv, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  for (int i_section = 0; i_section < NUM_SECTION_; i_section++) {
    for (int i_bin = 0; i_bin < NUM_BINS_(i_section); i_bin++) {
      aero_phase_idx -= NUM_PHASE_(i_section);
      if (aero_phase_idx < 0) {
        *radius = EFFECTIVE_RADIUS_(i_section, i_bin);
        i_section = NUM_SECTION_;
        break;
      }
    }
  }

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the particle number concentration \f$n\f$ (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 *
 * The modal mass number concentration is calculated for a log-normal
 * distribution where the geometric mean diameter (\f$\tilde{D}_n\f$) and
 * geometric standard deviation (\f$\tilde{\sigma}_g\f$) are set by the aerosol
 * model prior to solving the chemistry. Thus, all
 * \f$\frac{\partial n}{\partial y}\f$ are zero. The number concentration is
 * calculated according to the equation given in Table 1 of Zender
 * \cite Zender2002 :
 * \f[
 *      n = N_0 = \frac{6V_0}{\pi}\tilde{D}_n^{-3}e^{-9\tilde{\sigma}_g^2/2}
 * \f]
 * \f[
 *      V_0 = \sum_i{\rho_im_i}
 * \f]
 * where \f$\rho_i\f$ and \f$m_i\f$ are the density and total mass of species
 * \f$i\f$ in the specified mode.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param number_conc Particle number concentration, \f$n\f$
 *                    (\f$\mbox{\si{\#\per\cubic\centi\metre}}\f$)
 * \param partial_deriv \f$\frac{\partial n}{\partial y}\f$ where \f$y\f$ are
 *                      the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation
 */
__device__ void *aero_rep_gpu_modal_binned_mass_get_number_conc(int aero_phase_idx,
                                                 double *number_conc, double *partial_deriv, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  for (int i_section = 0; i_section < NUM_SECTION_ && aero_phase_idx >= 0;
       i_section++) {
    for (int i_bin = 0; i_bin < NUM_BINS_(i_section) && aero_phase_idx >= 0;
         i_bin++) {
      aero_phase_idx -= NUM_PHASE_(i_section);
      if (aero_phase_idx < 0) {
        *number_conc = NUMBER_CONC_(i_section, i_bin);
        i_section = NUM_SECTION_;
        break;
      }
    }
  }

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
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
__device__ void *aero_rep_gpu_modal_binned_mass_get_aero_conc_type(int aero_phase_idx,
                                                    int *aero_conc_type, void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  *aero_conc_type = 1;

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the total mass in an aerosol phase \f$m\f$ (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$)
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
__device__ void * aero_rep_gpu_modal_binned_mass_get_aero_phase_mass(int aero_phase_idx,
                                                      double *aero_phase_mass, double *aero_phase_avg_MW,
                                                      double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_section=0; i_section<NUM_SECTION_ && aero_phase_idx>=0;
       i_section++) {
    for (int i_phase=0; i_phase<NUM_PHASE_(i_section) && aero_phase_idx>=0;
         i_phase++) {
      for (int i_bin=0; i_bin<NUM_BINS_(i_section) && aero_phase_idx>=0;
           i_bin++) {
        if (aero_phase_idx==0) {
          *aero_phase_mass = PHASE_MASS_(i_section, i_phase, i_bin);
          *aero_phase_avg_MW = PHASE_AVG_MW_(i_section, i_phase, i_bin);
        }
        aero_phase_idx-=1;
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Advance the aerosol representation data pointer to the next aerosol representation
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol
 *         representation data
 */
__device__ void *aero_rep_gpu_modal_binned_mass_skip(void *aero_rep_data) {
  int *int_data = (int *) aero_rep_data;
  double *float_data = (double *) &(int_data[INT_DATA_SIZE_]);

  return (void *) &(float_data[FLOAT_DATA_SIZE_]);
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

}
