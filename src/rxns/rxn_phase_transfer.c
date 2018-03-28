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

#define _UNIV_GAS_CONST_ 8.314472

#define _del_H_ float_data[0]
#define _del_S_ float_data[1]
#define _Dg_ float_data[2]
#define _pre_c_rms_ float_data[3]
#define _A_ float_data[4]
#define _C_ float_data[5]
#define _c_rms_alpha_ float_data[6]
#define _equil_const_ float_data[7]
#define _NUM_AERO_PHASE_ int_data[0]
#define _GAS_SPEC_ int_data[1]
#define _NUM_INT_PROP_ 2
#define _NUM_FLOAT_PROP_ 8
#define _AERO_SPEC_(x) (int_data[_NUM_INT_PROP_ + x]-1)
#define _DERIV_ID_(x) int_data[_NUM_INT_PROP_ + _NUM_AERO_PHASE_ + x]
#define _JAC_ID_(x) int_data[_NUM_INT_PROP_ + 1 + 2*(_NUM_AERO_PHASE_) + x]
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+1+(5*_NUM_AERO_PHASE_))
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_phase_transfer_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  jac_struct[_GAS_SPEC_][_GAS_SPEC_] = true;
  for (int i_aero_spec = 0; i_aero_spec < _NUM_AERO_PHASE_; i_aero_spec++) {
    jac_struct[_AERO_SPEC_(i_aero_spec)][_REACT_] = true;
    jac_struct[_REACT_][_AERO_SPEC_(i_aero_spec)] = true;
    jac_struct[_AERO_SPEC_(i_aero_spec)][_AERO_SPEC_(i_aero_spec)] = true;
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
void * rxn_phase_transfer_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Update the time derivative ids
  _DERIV_ID_(0) = deriv_ids[_REACT_];
  for (int i=0; i < _NUM_AERO_PHASE_; i++)
	  _DERIV_ID_(i + 1) = deriv_ids[_AERO_SPEC_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  _JAC_ID_(i_jac++) = jac_ids[_REACT_][_REACT_];  
  for (int i_aero_spec = 0; i_aero_spec < _NUM_AERO_PHASE_; i_aero_spec++) {
      _JAC_ID_(i_jac++) = jac_ids[_AERO_SPEC_(i_aero_spec)][_REACT_];
      _JAC_ID_(i_jac++) = jac_ids[_REACT_][_AERO_SPEC_(i_aero_spec)];
      _JAC_ID_(i_jac++) = jac_ids[_AERO_SPEC_(i_aero_spec)][_AERO_SPEC_(i_aero_spec)];
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
void * rxn_phase_transfer_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the mass accomodation coefficient if the N* parameter
  // was provided, otherwise set it to 1.0
  realtype mass_acc = 1.0
  if (_del_H_!=0.0 || _del_S_!=0.0) {
    realtype del_G = _del_H_ - _TEMPERATURE_ * _del_S_; 
    mass_acc = exp(-del_G/(_UNIV_GAS_CONST_ * _TEMPERATURE_));
    mass_acc = mass_acc / (1.0 + mass_acc);
  }

  // Save c_rms * mass_acc for use in mass transfer rate calc
  _c_rms_alpha_ = _pre_c_rms_ * sqrt(_TEMPERATURE_) * mass_acc;

  // Calculate the Henry's Law equilibrium rate constant
  if (_B_==0.0) {
    _equil_const_ = _A_;
  } else {
    _equil_const_ = _A_ * exp(_B_ * (1.0 - _TEMPERATURE_/298.0);
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param state Pointer to the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_phase_transfer_calc_deriv_contrib(realtype *state, realtype *deriv,
		void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the reaction rate
  realtype rate = _RATE_CONSTANT_;
  for (int i_spec=0; i_spec<_NUM_REACT_; i_spec++) rate *= state[_REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<_NUM_REACT_; i_spec++, i_dep_var++) {
      if (_DERIV_ID_(i_dep_var) < 0) continue; 
      deriv[_DERIV_ID_(i_dep_var)] -= rate;
    }
    for (int i_spec=0; i_spec<_NUM_PROD_; i_spec++, i_dep_var++) {
      if (_DERIV_ID_(i_dep_var) < 0) continue; 
      deriv[_DERIV_ID_(i_dep_var)] += rate*_yield_(i_spec);
    }
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
void * rxn_phase_transfer_calc_jac_contrib(realtype *state, realtype *J,
		void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the reaction rate
  realtype rate = _RATE_CONSTANT_;
  for (int i_spec=0; i_spec<_NUM_REACT_; i_spec++) rate *= state[_REACT_(i_spec)];

  // Add contributions to the Jacobian
  if (rate!=ZERO) {
    int i_elem = 0;
    for (int i_ind=0; i_ind<_NUM_REACT_; i_ind++) {
      for (int i_dep=0; i_dep<_NUM_REACT_; i_dep++, i_elem++) {
	if (_JAC_ID_(i_elem) < 0) continue;
	J[_JAC_ID_(i_elem)] -= rate / state[_REACT_(i_ind)];
      }
      for (int i_dep=0; i_dep<_NUM_PROD_; i_dep++, i_elem++) {
	if (_JAC_ID_(i_elem) < 0) continue;
	J[_JAC_ID_(i_elem)] += _yield_(i_dep) * rate / state[_REACT_(i_ind)];
      }
    }
  }

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);

}

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_phase_transfer_skip(void *rxn_data)
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
void * rxn_phase_transfer_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nPhase Transfer reaction\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Return the reaction rate for the current conditions
 *
 * \param rxn_data Pointer to the reaction data
 * \param state Pointer to the state array
 * \param env Pointer to the environmental state array
 * \param rate Pointer to a double value to store the calculated rate
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_phase_transfer_get_rate(void *rxn_data, realtype *state, realtype *env, realtype *rate)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the reaction rate
  rxn_phase_transfer_update_env_state(env, rxn_data);
  *rate = _RATE_CONSTANT_;
  for (int i_spec=0; i_spec<_NUM_REACT_; i_spec++) *rate *= state[_REACT_(i_spec)];

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _A_
#undef _B_
#undef _C_
#undef _D_
#undef _E_
#undef _CONV_
#undef _RATE_CONSTANT_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _REACT_
#undef _PROD_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _yield_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
