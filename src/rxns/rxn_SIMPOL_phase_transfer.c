/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Phase Transfer reaction solver functions
 *
*/
/** \file
 * \brief Phase Transfer reaction solver functions
*/
#ifdef PMC_USE_SUNDIALS

#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

// Universal gas constant (J/mol/K)
#define _UNIV_GAS_CONST_ 8.314472
// Small number
#define _SMALL_NUMBER_ 1.0e-30

#define _del_H_ float_data[0]
#define _del_S_ float_data[1]
#define _Dg_ float_data[2]
#define _pre_c_rms_ float_data[3]
#define _B1_ float_data[4]
#define _B2_ float_data[5]
#define _B3_ float_data[6]
#define _B4_ float_data[7]
#define _c_rms_alpha_ float_data[8]
#define _equil_const_ float_data[9]
#define _CONV_ float_data[10]
#define _MW_ float_data[11]
#define _ug_m3_TO_ppm_ float_data[12]
#define _NUM_AERO_PHASE_ int_data[0]
#define _GAS_SPEC_ (int_data[1]-1)
#define _NUM_INT_PROP_ 2
#define _NUM_FLOAT_PROP_ 13
#define _AERO_SPEC_(x) (int_data[_NUM_INT_PROP_ + x]-1)
#define _AERO_ACT_ID_(x) (int_data[_NUM_INT_PROP_ + _NUM_AERO_PHASE_ + x])
#define _AERO_PHASE_ID_(x) (int_data[_NUM_INT_PROP_ + 2*(_NUM_AERO_PHASE_) + x]-1)
#define _AERO_REP_ID_(x) (int_data[_NUM_INT_PROP_ + 3*(_NUM_AERO_PHASE_) + x]-1)
#define _DERIV_ID_(x) (int_data[_NUM_INT_PROP_ + 4*(_NUM_AERO_PHASE_) + x])
#define _JAC_ID_(x) (int_data[_NUM_INT_PROP_ + 1 + 5*(_NUM_AERO_PHASE_) + x])
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+2+(8*_NUM_AERO_PHASE_))
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  jac_struct[_GAS_SPEC_][_GAS_SPEC_] = true;
  for (int i_aero_phase = 0; i_aero_phase < _NUM_AERO_PHASE_; i_aero_phase++) {
    jac_struct[_AERO_SPEC_(i_aero_phase)][_GAS_SPEC_] = true;
    jac_struct[_GAS_SPEC_][_AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[_AERO_SPEC_(i_aero_phase)][_AERO_SPEC_(i_aero_phase)] = true;
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}  

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data for finding sub model ids
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids, 
    int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Update the time derivative ids
  _DERIV_ID_(0) = deriv_ids[_GAS_SPEC_];
  for (int i=0; i < _NUM_AERO_PHASE_; i++)
	  _DERIV_ID_(i + 1) = deriv_ids[_AERO_SPEC_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  _JAC_ID_(i_jac++) = jac_ids[_GAS_SPEC_][_GAS_SPEC_];  
  for (int i_aero_phase = 0; i_aero_phase < _NUM_AERO_PHASE_; i_aero_phase++) {
      _JAC_ID_(i_jac++) = jac_ids[_AERO_SPEC_(i_aero_phase)][_GAS_SPEC_];
      _JAC_ID_(i_jac++) = jac_ids[_GAS_SPEC_][_AERO_SPEC_(i_aero_phase)];
      _JAC_ID_(i_jac++) = jac_ids[_AERO_SPEC_(i_aero_phase)][_AERO_SPEC_(i_aero_phase)];
    }

  // Find activity coefficient ids, if they exist
  // TODO Don't hard-code sub model ids
  for (int i_aero_phase = 0; i_aero_phase < _NUM_AERO_PHASE_; i_aero_phase++) {
    int aero_state_id = _AERO_SPEC_(i_aero_phase);
    _AERO_ACT_ID_(i_aero_phase) = 
      sub_model_get_parameter_id(model_data, 1, (void*) (&aero_state_id));
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Phase Transfer reaction this only involves recalculating the rate 
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the mass accomodation coefficient if the N* parameter
  // was provided, otherwise set it to 1.0
  realtype mass_acc = 1.0;
  if (_del_H_!=0.0 || _del_S_!=0.0) {
    realtype del_G = _del_H_ - _TEMPERATURE_K_ * _del_S_; 
    mass_acc = exp(-del_G/(_UNIV_GAS_CONST_ * _TEMPERATURE_K_));
    mass_acc = mass_acc / (1.0 + mass_acc);
  }

  // Save c_rms * mass_acc for use in mass transfer rate calc
  _c_rms_alpha_ = _pre_c_rms_ * sqrt(_TEMPERATURE_K_) * mass_acc;

  // SIMPOL.1 vapor pressure (Pa)
  realtype vp = _B1_ / _TEMPERATURE_K_
                + _B2_ + _B3_ * _TEMPERATURE_K_
                + _B4_ * log(_TEMPERATURE_K_);
  vp = 101325.0 * pow(10, vp);

  // Calculate the conversion from ug_x/m^3 -> ppm_x
  _ug_m3_TO_ppm_ = _CONV_ * _TEMPERATURE_K_ / _PRESSURE_PA_;

  // Calculate the partitioning coefficient K_eq (ppm_x/ug_x*ug_tot/kg_tot) such that for
  // partitioning species X at equilibrium:
  //   [X]_gas = [X]_aero * activity_coeff_X * K_eq * MW_tot_aero / [tot]_aero
  // where 'tot' indicates all species within an aerosol phase combined
  // with []_gas in (ppm) and []_aero in (ug/m^3)
  _equil_const_ = vp                    // (Pa_x*mol_tot/mol_x)
                  / _PRESSURE_PA_       // (1/Pa_air)
                  / _MW_                // (mol_x/kg_x)  
                  * 1.0e6;             // 1.0e6ppm_x*Pa_air/Pa_x * 1.0e-9kg_x/ug_x * 1.0e9ug_tot/kg_tot

  rxn_SIMPOL_phase_transfer_print( rxn_data );

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for phase_transfer reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_pre_calc(ModelData *model_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
		void *rxn_data)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<_NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    // FIXME check whether effective radius is the correct radius to use
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)
   
    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number concentration (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    realtype aero_phase_mass;
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_mass(
                  model_data,                   // model data
                  _AERO_REP_ID_(i_phase),       // aerosol representation index
                  _AERO_PHASE_ID_(i_phase),     // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  &aero_phase_avg_MW);          // average MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius < _SMALL_NUMBER_ || number_conc < _SMALL_NUMBER_
        || aero_phase_mass < _SMALL_NUMBER_) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*_Dg_) + 4.0*radius/(3.0*_c_rms_alpha_));
  
    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    // cond_rate / (  K_eq                                     * MW_tot           / mass_tot     )
    //     1/s      * ppm_x * mol_tot / ug_x * ug_tot / kg_tot * kg_tot / mol_tot * m^3 / ug_tot )
    realtype evap_rate = cond_rate / (_equil_const_ * aero_phase_avg_MW / aero_phase_mass);

    // Calculate gas-phase condensation rate (ppm/s)
    cond_rate *= state[_GAS_SPEC_];

    // Get the activity coefficient (if one exists)
    realtype act_coeff = 1.0;
    if (_AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_get_parameter_value(model_data, _AERO_ACT_ID_(i_phase));
    }

    // Calculate aerosol-phase evaporation rate (ug/m^3/s)
    evap_rate *= state[_AERO_SPEC_(i_phase)] * act_coeff;

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (_DERIV_ID_(0)>=0) {
      if (aero_conc_type==0) {
        // Scale the changes to the gas-phase by the number of particles for per-particle 
        // aerosol concentrations
        deriv[_DERIV_ID_(0)] += number_conc * (evap_rate * _ug_m3_TO_ppm_ - cond_rate);
      } else {
        // No scaling for aerosol concentrations with total mass per aerosol phase
        deriv[_DERIV_ID_(0)] += evap_rate * _ug_m3_TO_ppm_ - cond_rate;
      }
    }

    // Change in the aerosol-phase species is condensation - evaporation (ug/m^3/s)
    if (_DERIV_ID_(1+i_phase)>=0) deriv[_DERIV_ID_(1+i_phase)] += cond_rate / _ug_m3_TO_ppm_ - evap_rate;

  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param state Pointer to the state array
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_calc_jac_contrib(ModelData *model_data, realtype *J,
		void *rxn_data)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<_NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)
   
    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number concentration (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  _AERO_REP_ID_(i_phase),	// aerosol representation index
		  _AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // Get the total mass of the aerosol phase
    realtype aero_phase_mass;
    realtype aero_phase_avg_MW;
    aero_rep_get_aero_phase_mass(
                  model_data,                   // model data
                  _AERO_REP_ID_(i_phase),       // aerosol representation index
                  _AERO_PHASE_ID_(i_phase),     // aerosol phase index
                  &aero_phase_mass,             // total aerosol-phase mass
                  &aero_phase_avg_MW);          // average MW in the aerosol phase

    // If the radius, number concentration, or aerosol-phase mass are zero,
    // no transfer occurs
    if (radius < _SMALL_NUMBER_ || number_conc < _SMALL_NUMBER_
        || aero_phase_mass < _SMALL_NUMBER_) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*_Dg_) + 4.0*radius/(3.0*_c_rms_alpha_));
  
    // Calculate the evaporation rate constant (ppm_x*m^3/ug_x/s)
    // cond_rate / (  K_eq                                     * MW_tot           / mass_tot     )
    //     1/s      * ppm_x * mol_tot / ug_x * ug_tot / kg_tot * kg_tot / mol_tot * m^3 / ug_tot )
    realtype evap_rate = cond_rate / (_equil_const_ * aero_phase_avg_MW / aero_phase_mass);

    // Get the activity coefficient (if one exists)
    realtype act_coeff = 1.0;
    if (_AERO_ACT_ID_(i_phase)>-1) {
      act_coeff = sub_model_get_parameter_value(model_data, _AERO_ACT_ID_(i_phase));
    }

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (aero_conc_type==0) {
      // Scale the changes to the gas-phase by the number of particles for per-particle 
      // aerosol concentrations
      if (_JAC_ID_(1+i_phase*3+1)>=0) 
	      J[_JAC_ID_(1+i_phase*3+1)] += number_conc * evap_rate * _ug_m3_TO_ppm_ * act_coeff;
      if (_JAC_ID_(0)>=0) J[_JAC_ID_(0)] -= number_conc * cond_rate;
    } else {
      // No scaling for aerosol concentrations with total mass per aerosol phase
      if (_JAC_ID_(1+i_phase*3+1)>=0) J[_JAC_ID_(1+i_phase*3+1)] += evap_rate * _ug_m3_TO_ppm_ 
              * act_coeff;
      if (_JAC_ID_(0)>=0) J[_JAC_ID_(0)] -= cond_rate;
    }

    // Change in the aerosol-phase species is condensation - evaporation (ug/m^3/s)
    if (_JAC_ID_(1+i_phase*3)>=0) J[_JAC_ID_(1+i_phase*3)] += cond_rate / _ug_m3_TO_ppm_;
    if (_JAC_ID_(1+i_phase*3+2)>=0) J[_JAC_ID_(1+i_phase*3+2)] -= evap_rate * act_coeff; 
  
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nSIMPOL.1 Phase Transfer reaction\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Return the reaction rate for the current conditions
 *
 * Phase-transfer reactions have a rate for each aerosol phase they affect
 * TODO figure out how to include these reactions in the rate functions
 *
 * \param rxn_data Pointer to the reaction data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \param rate Pointer to a double value to store the calculated rate
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_SIMPOL_phase_transfer_get_rate(void *rxn_data, realtype *state, realtype *env, realtype *rate)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *rate = 0.0;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _UNIV_GAS_CONST_
#undef _SMALL_NUMBER_

#undef _del_H_
#undef _del_S_
#undef _Dg_
#undef _pre_c_rms_
#undef _B1_
#undef _B2_
#undef _B3_
#undef _B4_
#undef _c_rms_alpha_
#undef _equil_const_
#undef _CONV_
#undef _MW_
#undef _ug_m3_TO_ppm_
#undef _NUM_AERO_PHASE_
#undef _GAS_SPEC_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _AERO_SPEC_
#undef _AERO_ACT_ID_
#undef _AERO_PHASE_ID_
#undef _AERO_REP_ID_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
