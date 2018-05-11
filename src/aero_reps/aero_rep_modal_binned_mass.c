/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Modal mass aerosol representation functions
 *
 */
/** \file
 * \brief Modal mass aerosol representation functions
 */
#ifdef PMC_USE_SUNDIALS

#include "../aero_rep_solver.h"

// TODO Lookup environmental indicies during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define BINNED 1
#define MODAL 2

#define _NUM_SECTION_ (int_data[0])
#define _INT_DATA_SIZE_ (int_data[1])
#define _FLOAT_DATA_SIZE_ (int_data[2])
#define _NUM_INT_PARAM_ 3
#define _NUM_FLOAT_PARAM_ 0
#define _MODE_INT_PARAM_LOC_(x) (int_data[_NUM_INT_PARAM_+x]-1)
#define _MODE_FLOAT_PARAM_LOC_(x) (int_data[_NUM_INT_PARAM_+_NUM_SECTION_+x]-1)
#define _SECTION_TYPE_(x) (int_data[_MODE_INT_PARAM_LOC_(x)])

// For sections, _NUM_BINS_ = 1
#define _NUM_BINS_(x) (int_data[_MODE_INT_PARAM_LOC_(x)+1])

#define _NUM_PHASE_(x) (int_data[_MODE_INT_PARAM_LOC_(x)+2])
#define _PHASE_INT_PARAM_LOC_(x,y) (int_data[_MODE_INT_PARAM_LOC_(x)+3+y]-1)
#define _PHASE_FLOAT_PARAM_LOC_(x,y) (int_data[_MODE_INT_PARAM_LOC_(x)+3+_NUM_PHASE_(x)+y]-1)
#define _NUM_SPEC_(x,y) (int_data[_PHASE_INT_PARAM_LOC_(x,y)])

// Species state ids - for modes, b=0
#define _SPEC_STATE_ID_(x,y,b,z) (int_data[_PHASE_INT_PARAM_LOC_(x,y)+b*(_NUM_SPEC_(x,y))+1+z]-1)

// GMD and bin diameter are stored in the same position - for modes, b=0
#define _GMD_(x,b) (float_data[_MODE_FLOAT_PARAM_LOC_(x)+b*4])
#define _BIN_Dp_(x,b) (float_data[_MODE_FLOAT_PARAM_LOC_(x)+b*4])

// GSD - only used for modes, b=0
#define _GSD_(x,b) (float_data[_MODE_FLOAT_PARAM_LOC_(x)+b*4+1])

// Real-time number concentration - used for modes and bins - for modes, b=0
#define _NUMBER_CONC_(x,b) (float_data[_MODE_FLOAT_PARAM_LOC_(x)+b*4+2])

// Real-time effective radius - only used for modes, b=0
#define _EFFECTIVE_RADIUS_(x,b) (float_data[_MODE_FLOAT_PARAM_LOC_(x)+b*4+3])

// Real-time phase mass (ug/m^3) - used for modes and bins - for modes, b=0
#define _AERO_PHASE_MASS_(x,y,b) (float_data[_PHASE_FLOAT_PARAM_LOC_(x,y)+b])

// Species density
#define _DENSITY_(x,y,z) (float_data[_PHASE_FLOAT_PARAM_LOC_(x,y)+_NUM_BINS_(x)+z])

// Update types (These must match values in aero_rep_modal_binned_mass.F90)
#define UPDATE_GMD 0
#define UPDATE_GSD 1

/** \brief Flag elements on the state array used by this aerosol representation
 *
 * The modal mass aerosol representation functions do not use state array values
 *
 * \param aero_rep_data A pointer to the aerosol representation data
 * \param state_flags Array of flags indicating state array elements used
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_modal_binned_mass_get_dependencies(void *aero_rep_data, bool *state_flags)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update aerosol representation data for new environmental conditions
 *
 * The modal mass aerosol representation does not use environmental conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_update_env_state(double *env_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update aerosol representation data for a new state
 *
 * The modal mass aerosol representation recalculates effective radius and number
 * concentration for each new state. 
 *
 * \param model_data Pointer to the model data, including the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_update_state(ModelData *model_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Loop through the modes and calculate effective radius and number concentration
  for (int i_section=0; i_section<_NUM_SECTION_; i_section++) {

    realtype volume, mass;
    switch (_SECTION_TYPE_(i_section)) {
      
      // Mode
      case (MODAL) : 
    
        // Sum the volumes of each species in the mode
        volume = 0.0;
        for (int i_phase=0; i_phase<_NUM_PHASE_(i_section); i_phase++) {
          mass = 0.0;
          for (int i_spec=0; i_spec<_NUM_SPEC_(i_section, i_phase); i_spec++){
            volume += _DENSITY_(i_section, i_phase, i_spec) * 
                    model_data->state[_SPEC_STATE_ID_(i_section, i_phase, 0, i_spec)];
            mass += model_data->state[_SPEC_STATE_ID_(i_section, i_phase, 0, i_spec)];
          }

          // Set the aerosol-phase mass
          _AERO_PHASE_MASS_(i_section, i_phase, 0) = mass;
        }
 
        // Calculate effective radius
        _EFFECTIVE_RADIUS_(i_section,0) = exp(_GMD_(i_section,0) + 5.0/2.0*_GSD_(i_section,0));

        // Calculate the number concentration based on the total mode volume  
        _NUMBER_CONC_(i_section,0) = volume * 3.0 / (4.0*M_PI) * 
                exp(2.0*_GMD_(i_section,0) + 9.0/2.0*_GSD_(i_section,0));

        break;

      // Bins
      case (BINNED) :

        // Loop through the bins
        for (int i_bin=0; i_bin<_NUM_BINS_(i_section); i_bin++) {
          
          // Sum the volumes of each species in the bin
          volume = 0.0;
          for (int i_phase=0; i_phase<_NUM_PHASE_(i_section); i_phase++) {
            mass = 0.0;
            for (int i_spec=0; i_spec<_NUM_SPEC_(i_section, i_phase); i_spec++) {
              volume += _DENSITY_(i_section, i_phase, i_spec) * 
                      model_data->state[_SPEC_STATE_ID_(i_section, i_phase, i_bin, i_spec)];
              mass += model_data->state[_SPEC_STATE_ID_(i_section, i_phase, 0, i_spec)];
            }

            // Set the aerosol-phase mass
            _AERO_PHASE_MASS_(i_section, i_phase, i_bin) = mass;
          }

          // Calculate the number concentration based on the total bin volume
          _NUMBER_CONC_(i_section, i_bin) = volume * 3.0 / (4.0*M_PI) *
                  pow(_BIN_Dp_(i_section, i_bin)/1.0, 3.0);
        }

        break;
    }

  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the effective particle radius
 *
 * The modal mass effective radius is calculated for a log-normal distribution
 * where the geometric mean diameter (\f$\mu\f$) and geometric standard
 * deviation (\f$S\f$) are set by the aerosol model prior to solving the
 * chemistry. Thus, all dr/dy are zero. The effective radius is calculated as:
 *
 * \f$[
 *      r_e = e^{\mu+\frac{5}{2}S}
 * ]\f$
 *
 * FIXME Check and add reference
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param radius Effective particle radius (m)
 * \param partial_deriv dr/dy where y are species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_get_effective_radius(int aero_phase_idx, double *radius, 
		double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  for (int i_section=0; i_section<_NUM_SECTION_; i_section++) {
    for (int i_bin=0; i_bin<_NUM_BINS_(i_section); i_bin++) {
      aero_phase_idx-=_NUM_PHASE_(i_section);
      if (aero_phase_idx<0) {
        *radius = _EFFECTIVE_RADIUS_(i_section, i_bin);
        i_section = _NUM_SECTION_;
        break;
      }
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the particle number concentration
 *
 * The modal mass number concentration is calculated using the total mass in
 * the specified mode assuming a log-normal distribution of fixed geometric
 * mean diameter (\f$\mu\f$) and geometric standard deviation (\f$S\f$):
 * \f$[
 *      N = \frac{3}{4\pi}Ve^{-3\mu-\frac{9}{2}S}
 *      V = \sum_i{\rho_im_i}
 * \f$]
 * where \f$\rho_i\f$ and \f$m_i\f$ are the density and total mass of species
 * \f$i\f$ in the specified mode.
 *
 * FIXME Check and add reference
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param number_conc Particle number concentration (#/cm^3)
 * \param partial_deriv dn/dy where y are the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_get_number_conc(int aero_phase_idx, double *number_conc, 
		double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  for (int i_section=0; i_section<_NUM_SECTION_; i_section++) {
    for (int i_bin=0; i_bin<_NUM_BINS_(i_section); i_bin++) {
      aero_phase_idx-=_NUM_PHASE_(i_section);
      if (aero_phase_idx<0) {
        *number_conc = _NUMBER_CONC_(i_section, i_bin);
        i_section = _NUM_SECTION_;
        break;
      }
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Get the type of aerosol concentration type used.
 *
 * Modal mass concentrations are per-mode.
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_conc_type Pointer to int that will hold the concentration type code
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_get_aero_conc_type(int aero_phase_idx, int *aero_conc_type, 
		void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *aero_conc_type = 1;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}
  
/** Get the total mass in an aerosol phase
 *
 * \param aero_phase_idx Index of the aerosol phase within the representation
 * \param aero_phase_mass Total mass in the aerosol phase (ug/m^3)
 * \param partial_deriv dn/dy where y are the species on the state array
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_get_aero_phase_mass(int aero_phase_idx, double *aero_phase_mass,
		double *partial_deriv, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  for (int i_section=0; i_section<_NUM_SECTION_; i_section++) {
    for (int i_phase=0; i_phase<_NUM_PHASE_(i_section); i_section++) {
      for (int i_bin=0; i_bin<_NUM_BINS_(i_section); i_bin++) {
        if (aero_phase_idx==0) {
          *aero_phase_mass = _AERO_PHASE_MASS_(i_section, i_phase, i_bin);
          i_section = _NUM_SECTION_;
          break;
        }
        aero_phase_idx-=1;
      }
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** Update the aerosol representation data
 *
 *  The modal mass aerosol representation has two update types:
 *
 *  UPDATE_GMD : where the update data should point to an int indicating the
 *  mode id to update followed by a floating-point variable holding the new
 *  geometric mean diameter
 *
 *  UPDATE_GSD : where the update data should point to an int indicating the
 *  mode id to update followed by a floating-point variable holding the new
 *  geometric standard deviation
 *
 * \param update_type The type of update to perform
 * \param update_data Pointer to the data required for the update
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation
 */
void * aero_rep_modal_binned_mass_update_data(int update_type, void *update_data, void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  int *i_section = (int*)update_data;

  if (_SECTION_TYPE_(*i_section)==MODAL) {  
    switch (update_type) {
      case UPDATE_GMD :
        _GMD_(*i_section,0) = *((realtype*)(i_section+1));
        break;
      case UPDATE_GSD :
        _GSD_(*i_section,0) = *((realtype*)(i_section+1));
        break;
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the Modal Mass reaction parameters
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_modal_binned_mass_print(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nModal mass aerosol representation\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Advance the aerosol representation data pointer to the next aerosol representation
 *
 * \param aero_rep_data Pointer to the aerosol representation data
 * \return The aero_rep_data pointer advanced by the size of the aerosol representation data
 */
void * aero_rep_modal_binned_mass_skip(void *aero_rep_data)
{
  int *int_data = (int*) aero_rep_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef BINNED
#undef MODAL
#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_
#undef _NUM_SECTION_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_
#undef _NUM_INT_PARAM_
#undef _NUM_FLOAT_PARAM_
#undef _MODE_INT_PARAM_LOC_
#undef _MODE_FLOAT_PARAM_LOC_
#undef _SECTION_TYPE_
#undef _NUM_BINS_
#undef _NUM_PHASE_
#undef _PHASE_INT_PARAM_LOC_
#undef _PHASE_FLOAT_PARAM_LOC_
#undef _NUM_SPEC_
#undef _SPEC_STATE_ID_
#undef _GMD_
#undef _BIN_Dp_
#undef _GSD_
#undef _NUMBER_CONC_
#undef _EFFECTIVE_RADIUS_
#undef _AERO_PHASE_MASS_
#undef _DENSITY_

#endif
