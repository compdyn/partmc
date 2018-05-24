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
#define _C_ (float_data[1])
#define _k_reverse_ (float_data[2])
#define _k_forward_ (float_data[3])
#define _NUM_INT_PROP_ 3
#define _NUM_FLOAT_PROP_ 4
#define _REACT_(x) (int_data[_NUM_INT_PROP_+x]-1)
#define _PROD_(x) (int_data[_NUM_INT_PROP_+_NUM_REACT_*_NUM_AERO_PHASE_+x]-1)
#define _WATER_(x) (int_data[_NUM_INT_PROP_+(_NUM_REACT_+_NUM_PROD_)*_NUM_AERO_PHASE_+x]-1)
#define _ACTIVITY_COEFF_(x) (int_data[_NUM_INT_PROP_+(_NUM_REACT_+_NUM_PROD_+1)*_NUM_AERO_PHASE_+x]-1)
#define _DERIV_ID_(x) (int_data[_NUM_INT_PROP_+(_NUM_REACT_+_NUM_PROD_+2)*_NUM_AERO_PHASE_+x])
#define _JAC_ID_(x) (int_data[_NUM_INT_PROP_+(2*(_NUM_REACT_+_NUM_PROD_)+2)*_NUM_AERO_PHASE_+x])
#define _mass_frac_TO_M_(x) (float_data[_NUM_FLOAT_PROP_+x])
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+((_NUM_REACT_+_NUM_PROD_)*(_NUM_REACT_+_NUM_PROD_+3)+2)*_NUM_AERO_PHASE_)
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_+_NUM_PROD_+_NUM_REACT_)

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
    
    // Add dependence on reactants for reactants and products (forward reaction) 
    for (int i_react_ind = i_phase*_NUM_REACT_; i_react_ind < (i_phase+1)*_NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        jac_struct[_REACT_(i_react_dep)][_REACT_(i_react_ind)] = true;
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        jac_struct[_PROD_(i_prod_dep)][_REACT_(i_react_ind)] = true;
    }

    // Add dependence on products for reactants and products (reverse reaction) 
    for (int i_prod_ind = i_phase*_NUM_PROD_; i_prod_ind < (i_phase+1)*_NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        jac_struct[_REACT_(i_react_dep)][_PROD_(i_prod_ind)] = true;
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        jac_struct[_PROD_(i_prod_dep)][_PROD_(i_prod_ind)] = true;
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
      jac_struct[_REACT_(i_react_dep)][_WATER_(i_phase)] = true;
    for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
      jac_struct[_PROD_(i_prod_dep)][_WATER_(i_phase)] = true;

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
    
    // Add dependence on reactants for reactants and products (forward reaction) 
    for (int i_react_ind = i_phase*_NUM_REACT_; i_react_ind < (i_phase+1)*_NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_react_dep)][_REACT_(i_react_ind)];
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_prod_dep)][_REACT_(i_react_ind)];
    }

    // Add dependence on products for reactants and products (reverse reaction) 
    for (int i_prod_ind = i_phase*_NUM_PROD_; i_prod_ind < (i_phase+1)*_NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_react_dep)][_PROD_(i_prod_ind)];
      for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
        _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_prod_dep)][_PROD_(i_prod_ind)];
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*_NUM_REACT_; i_react_dep < (i_phase+1)*_NUM_REACT_; i_react_dep++)
      _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_react_dep)][_WATER_(i_phase)];
    for (int i_prod_dep = i_phase*_NUM_PROD_; i_prod_dep < (i_phase+1)*_NUM_PROD_; i_prod_dep++)
      _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_prod_dep)][_WATER_(i_phase)];

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

  // Calculate the equilibrium constant (assumes reactant and product concentrations in M)
  realtype equil_const;
  if (_C_==0.0) {
    equil_const = _A_;
  } else {
    equil_const = _A_ * exp(_C_ * (1.0/_TEMPERATURE_K_ - 1.0/298.0));
  }

  // Set the forward rate constant
  _k_forward_ = equil_const * _k_reverse_;

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

    // If no aerosol water is present, no reaction occurs
    if (state[_WATER_(i_phase)] < _SMALL_NUMBER_) {
      i_deriv += _NUM_REACT_ + _NUM_PROD_;
      continue;
    }

    // Calculate the forward rate (M/s)
    realtype forward_rate = _k_forward_;
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      forward_rate *= state[_REACT_(i_phase*_NUM_REACT_+i_react)] * 
              _mass_frac_TO_M_(i_react) / state[_WATER_(i_phase)];
    }

    // Calculate the reverse rate (M/s)
    realtype reverse_rate = _k_reverse_;
    for (int i_prod = 0; i_prod < _NUM_PROD_; i_prod++) {
      reverse_rate *= state[_PROD_(i_phase*_NUM_PROD_+i_prod)] * 
              _mass_frac_TO_M_(_NUM_REACT_+i_prod) / state[_WATER_(i_phase)];
    }
    if (_ACTIVITY_COEFF_(i_phase)>=0) reverse_rate *= 
            state[_ACTIVITY_COEFF_(i_phase)];

    // Reactants change as (reverse - forward) (ug/m3/s)
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      if (_DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[_DERIV_ID_(i_deriv++)] += (reverse_rate - forward_rate) /
	      _mass_frac_TO_M_(i_react) * state[_WATER_(i_phase)];
    }

    // Products change as (forward - reverse) (ug/m3/s)
    for (int i_prod = 0; i_prod < _NUM_PROD_; i_prod++) {
      if (_DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[_DERIV_ID_(i_deriv++)] += (forward_rate - reverse_rate) /
	      _mass_frac_TO_M_(_NUM_REACT_+i_prod) * state[_WATER_(i_phase)];
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

    // If not aerosol water is present, no reaction occurs
    if (state[_WATER_(i_phase)] < _SMALL_NUMBER_) {
      i_jac += (_NUM_REACT_ + _NUM_PROD_) * (_NUM_REACT_ + _NUM_PROD_ + 1);
      continue;
    }

    // Calculate the forward rate (M/s)
    realtype forward_rate = _k_forward_;
    for (int i_react = 0; i_react < _NUM_REACT_; i_react++) {
      forward_rate *= state[_REACT_(i_phase*_NUM_REACT_+i_react)] * 
              _mass_frac_TO_M_(i_react) / state[_WATER_(i_phase)];
    }

    // Calculate the reverse rate (M/s)
    realtype reverse_rate = _k_reverse_;
    for (int i_prod = 0; i_prod < _NUM_PROD_; i_prod++) {
      reverse_rate *= state[_PROD_(i_phase*_NUM_PROD_+i_prod)] * 
              _mass_frac_TO_M_(_NUM_REACT_+i_prod) / state[_WATER_(i_phase)];
    }
    if (_ACTIVITY_COEFF_(i_phase)>=0) reverse_rate *= 
            state[_ACTIVITY_COEFF_(i_phase)];

    // No Jac contributions to add if the rates are zero
    if ((forward_rate - reverse_rate)==0.0) {
      i_jac += (_NUM_REACT_ + _NUM_PROD_) * (_NUM_REACT_ + _NUM_PROD_ + 1);
      continue;
    }

    // Add dependence on reactants for reactants and products (forward reaction) 
    for (int i_react_ind = 0; i_react_ind < _NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < _NUM_REACT_; i_react_dep++) {
	if (_JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] += (-forward_rate) / state[_REACT_(i_phase*_NUM_REACT_+i_react_ind)] / 
		_mass_frac_TO_M_(i_react_dep) * state[_WATER_(i_phase)];
      }
      for (int i_prod_dep = 0; i_prod_dep < _NUM_PROD_; i_prod_dep++) {
	if (_JAC_ID_(i_jac)<0 || forward_rate==0.0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] += (forward_rate) / state[_REACT_(i_phase*_NUM_REACT_+i_react_ind)] / 
		_mass_frac_TO_M_(_NUM_REACT_ + i_prod_dep) * state[_WATER_(i_phase)];
      }
    }

    // Add dependence on products for reactants and products (reverse reaction) 
    for (int i_prod_ind = 0; i_prod_ind < _NUM_PROD_; i_prod_ind++) {
      for (int i_react_dep = 0; i_react_dep < _NUM_REACT_; i_react_dep++) {
	if (_JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] += (reverse_rate) / state[_PROD_(i_phase*_NUM_PROD_+i_prod_ind)] / 
		_mass_frac_TO_M_(i_react_dep) * state[_WATER_(i_phase)];
      }
      for (int i_prod_dep = 0; i_prod_dep < _NUM_PROD_; i_prod_dep++) {
	if (_JAC_ID_(i_jac)<0 || reverse_rate==0.0) {i_jac++; continue;}
        J[_JAC_ID_(i_jac++)] += (-reverse_rate) / state[_PROD_(i_phase*_NUM_PROD_+i_prod_ind)] / 
		_mass_frac_TO_M_(_NUM_REACT_ + i_prod_dep) * state[_WATER_(i_phase)];
      }
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = 0; i_react_dep < _NUM_REACT_; i_react_dep++) {
      if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[_JAC_ID_(i_jac++)] += (reverse_rate-forward_rate) / _mass_frac_TO_M_(i_react_dep);
    }
    for (int i_prod_dep = 0; i_prod_dep < _NUM_PROD_; i_prod_dep++) {
      if (_JAC_ID_(i_jac)<0) {i_jac++; continue;}
      J[_JAC_ID_(i_jac++)] += (forward_rate-reverse_rate) / _mass_frac_TO_M_(_NUM_REACT_ + i_prod_dep);
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
#undef _C_
#undef _k_reverse_
#undef _k_forward_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _REACT_
#undef _PROD_
#undef _WATER_
#undef _ACTIVITY_COEFF_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _mass_frac_TO_M_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
