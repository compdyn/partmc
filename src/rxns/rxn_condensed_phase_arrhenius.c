/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Condensed Phase Arrhenius reaction solver functions
 *
*/
/** \file
 * \brief Condensed Phase Arrhenius reaction solver functions
*/
#ifdef PMC_USE_SUNDIALS

#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

// Small number
#define _SMALL_NUMBER_ 1.0e-30

#define _NUM_REACT_ (int_data[0])
#define _NUM_PROD_ (int_data[1])
#define _NUM_AERO_PHASE_ (int_data[2])
#define _A_ (float_data[0])
#define _B_ (float_data[1])
#define _C_ (float_data[2])
#define _D_ (float_data[3])
#define _E_ (float_data[4])
#define _RATE_CONSTANT_ (float_data[5])
#define _NUM_INT_PROP_ 3
#define _NUM_FLOAT_PROP_ 6
#define _REACT_(x) (int_data[_NUM_INT_PROP_+x]-1)
#define _PROD_(x) (int_data[_NUM_INT_PROP_+_NUM_REACT_*_NUM_AERO_PHASE_+x]-1)
#define _WATER_(x) (int_data[_NUM_INT_PROP_+(_NUM_REACT_+_NUM_PROD_)*_NUM_AERO_PHASE_+x]-1)
#define _DERIV_ID_(x) (int_data[_NUM_INT_PROP_+(_NUM_REACT_+_NUM_PROD_+1)*_NUM_AERO_PHASE_+x])
#define _JAC_ID_(x) (int_data[_NUM_INT_PROP_+(2*(_NUM_REACT_+_NUM_PROD_)+1)*_NUM_AERO_PHASE_+x])
#define _yield_(x) (float_data[_NUM_FLOAT_PROP_+x])
#define _ugm3_TO_molm3_(x) (float_data[_NUM_FLOAT_PROP_+_NUM_REACT_+_NUM_PROD_+x])
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+((_NUM_REACT_+_NUM_PROD_)*(_NUM_REACT_+3)+1)*_NUM_AERO_PHASE_)
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_+2*(_NUM_PROD_+_NUM_REACT_))

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Loop over all the instances of the specified phase
  for (int i_phase = 0; i_phase < _NUM_AERO_PHASE_; i_phase++) {
    
    // Add dependence on reactants for reactants and products
    for (int i_react_ind = i_phase*_NUM_REACT_; i_react_ind < (i_phase+1)*_NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        jac_struct[_REACT_(i_react_dep)][_REACT_(i_react_ind)] = true;
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        jac_struct[_PROD_(i_prod_dep)][_REACT_(i_react_ind)] = true;
    }

    // Add dependence on aerosol-phase water for reactants and products in aqueous reactions
    if (_WATER_(i_phase)>=0) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        jac_struct[_REACT_(i_react_dep)][_WATER_(i_phase)] = true;
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        jac_struct[_PROD_(i_prod_dep)][_WATER_(i_phase)] = true;
    }

  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}  

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Update the time derivative ids
  for (int i_phase = 0, i_deriv = 0; i_phase < _NUM_AERO_PHASE_; i_phase++) {
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++)
      _DERIV_ID_(i_deriv++) = deriv_ids[_REACT_(i_phase*_NUM_REACT_+i_react)];
    for (int i_prod = 0; i_prod < _NUM_PROD_; i_prod++)
      _DERIV_ID_(i_deriv++) = deriv_ids[_PROD_(i_phase*_NUM_PROD_+i_prod)];
  }

  // Update the Jacobian ids
  for (int i_phase = 0, i_jac = 0; i_phase < _NUM_AERO_PHASE_; i_phase++) {
    
    // Add dependence on reactants for reactants and products 
    for (int i_react_ind = i_phase*_NUM_REACT_; i_react_ind < (i_phase+1)*_NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_react_dep)][_REACT_(i_react_ind)];
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_prod_dep)][_REACT_(i_react_ind)];
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
      if (_WATER_(i_phase)>=0) {
        _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_react_dep)][_WATER_(i_phase)];
      } else {
        _JAC_ID_(i_jac++) = -1;
      }
    for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
      if (_WATER_(i_phase)>=0) {
        _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_prod_dep)][_WATER_(i_phase)];
      } else {
        _JAC_ID_(i_jac++) = -1;
      }

  }
  
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Update reaction data for new environmental conditions
 *
 * For Condensed Phase Arrhenius reaction this only involves recalculating the 
 * forward rate constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the rate constant in (M or mol/m3)
  // k = A*exp(C/T) * (T/D)^B * (1+E*P)
  _RATE_CONSTANT_ = _A_ * SUNRexp(_C_/_TEMPERATURE_K_)
          * (_B_==ZERO ? ONE : SUNRpowerR(_TEMPERATURE_K_/_D_, _B_))
          * (_E_==ZERO ? ONE : (ONE + _E_*_PRESSURE_PA_));

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for condensed_phase_arrhenius reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_pre_calc(ModelData *model_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param state Pointer to the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
		void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<_NUM_AERO_PHASE_; i_phase++) {

    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    realtype unit_conv = 1.0;
    if (_WATER_(i_phase)>=0) {
      unit_conv = state[_WATER_(i_phase)] * 1.0e-9; // convert from ug/m3 -> L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction occurs
      if (unit_conv < _SMALL_NUMBER_) {
        i_deriv += _NUM_REACT_ + _NUM_PROD_;
        continue;
      }
      unit_conv = 1.0/unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    realtype rate = _RATE_CONSTANT_;
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      rate *= state[_REACT_(i_phase*_NUM_REACT_+i_react)] * 
              _ugm3_TO_molm3_(i_react) * unit_conv;
    }

    // Reactant change
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      if (_DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[_DERIV_ID_(i_deriv++)] -= rate /
	      (_ugm3_TO_molm3_(i_react) * unit_conv);
    }

    // Products change
    for (int i_prod = 0; i_prod < _NUM_PROD_; i_prod++) {
      if (_DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[_DERIV_ID_(i_deriv++)] += rate * _yield_(i_prod) /
	      (_ugm3_TO_molm3_(_NUM_REACT_+i_prod) * unit_conv);
    }

  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param state Pointer to the state array
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_calc_jac_contrib(ModelData *model_data, realtype *J,
		void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase=0, i_jac = 0; i_phase<_NUM_AERO_PHASE_; i_phase++) {

    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    realtype unit_conv = 1.0;
    if (_WATER_(i_phase)>=0) {
      unit_conv = state[_WATER_(i_phase)] * 1.0e-9; // convert from ug/m3 -> L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction occurs
      if (unit_conv < _SMALL_NUMBER_) {
        i_jac += (_NUM_REACT_ + _NUM_PROD_) * (_NUM_REACT_ + 1);
        continue;
      }
      unit_conv = 1.0/unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    realtype rate = _RATE_CONSTANT_;
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      rate *= state[_REACT_(i_phase*_NUM_REACT_+i_react)] * 
              _ugm3_TO_molm3_(i_react) * unit_conv;
    }

    // No Jac contributions to add if the rate is zero
    if (rate==0.0) {
      i_jac += (_NUM_REACT_ + _NUM_PROD_) * (_NUM_REACT_ + 1);
      continue;
    }

    // Add dependence on reactants for reactants and products 
    for (int i_react_ind = 0; i_react_ind < _NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < _NUM_REACT_; i_react_dep++) {
	if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] -= rate / state[_REACT_(i_phase*_NUM_REACT_+i_react_ind)] / 
	        (_ugm3_TO_molm3_(i_react_dep) * unit_conv);
      }
      for (int i_prod_dep = 0; i_prod_dep < _NUM_PROD_; i_prod_dep++) {
	if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] += rate / state[_REACT_(i_phase*_NUM_REACT_+i_react_ind)] / 
	        (_ugm3_TO_molm3_(_NUM_REACT_+i_prod_dep) * unit_conv);
      }
    }

    // Add dependence on aerosol-phase water for reactants and products in aqueous reactions
    for (int i_react_dep = 0; i_react_dep < _NUM_REACT_; i_react_dep++) {
      if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[_JAC_ID_(i_jac++)] += (_NUM_REACT_-1) *rate / state[_WATER_(i_phase)] /
	        (_ugm3_TO_molm3_(i_react_dep) * unit_conv);
    }
    for (int i_prod_dep = 0; i_prod_dep < _NUM_PROD_; i_prod_dep++) {
      if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[_JAC_ID_(i_jac++)] -= (_NUM_REACT_-1) * rate / state[_WATER_(i_phase)] /
	        (_ugm3_TO_molm3_(_NUM_REACT_+i_prod_dep) * unit_conv);
    }

  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the Condensed Phase Arrhenius reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nCondensed Phase Arrhenius reaction\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Return the reaction rate for the current conditions
 *
 * Condensed-phase Arrhenius reactions have a rate for each aerosol phase they affect
 * TODO figure out how to include these reactions in the rate functions
 *
 * \param rxn_data Pointer to the reaction data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \param rate Pointer to a double value to store the calculated rate
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_condensed_phase_arrhenius_get_rate(void *rxn_data, realtype *state, realtype *env, realtype *rate)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  *rate = 0.0;

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _SMALL_NUMBER_

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _NUM_AERO_PHASE_
#undef _A_
#undef _B_
#undef _C_
#undef _D_
#undef _E_
#undef _RATE_CONSTANT_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _REACT_
#undef _PROD_
#undef _WATER_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _yield_
#undef _ugm3_TO_molm3_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
