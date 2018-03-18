/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * CMAQ_H2O2 reaction solver functions
 *
*/
/** \file
 * \brief CMAQ_H2O2 reaction solver functions
*/
#include "../rxn_solver.h"

// TODO Lookup environmental indicies during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define _NUM_REACT_ int_data[0]
#define _NUM_PROD_ int_data[1]
#define _k1_A_ float_data[0]
#define _k1_B_ float_data[1]
#define _k1_C_ float_data[2]
#define _k2_A_ float_data[3]
#define _k2_B_ float_data[4]
#define _k2_C_ float_data[5]
#define _CONV_ float_data[6]
#define _RATE_CONSTANT_ float_data[7]
#define _NUM_INT_PROP_ 2
#define _NUM_FLOAT_PROP_ 8
#define _REACT_(x) (int_data[_NUM_INT_PROP_ + x]-1)
#define _PROD_(x) (int_data[_NUM_INT_PROP_ + _NUM_REACT_ + x]-1)
#define _DERIV_ID_(x) int_data[_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x]
#define _JAC_ID_(x) int_data[_NUM_INT_PROP_ + 2*(_NUM_REACT_+_NUM_PROD_) + x]
#define _yield_(x) float_data[_NUM_FLOAT_PROP_ + x]
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_+(_NUM_REACT_+2)*(_NUM_REACT_+_NUM_PROD_))
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_+_NUM_PROD_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_CMAQ_H2O2_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  for (int i_ind = 0; i_ind < _NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < _NUM_REACT_; i_dep++) {
      jac_struct[_REACT_(i_dep)][_REACT_(i_ind)] = true;
    }
    for (int i_dep = 0; i_dep < _NUM_PROD_; i_dep++) {
      jac_struct[_PROD_(i_dep)][_REACT_(i_ind)] = true;
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
void * rxn_CMAQ_H2O2_update_ids(int *deriv_ids, int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Update the time derivative ids
  for (int i=0; i < _NUM_REACT_; i++)
	  _DERIV_ID_(i) = deriv_ids[_REACT_(i)];
  for (int i=0; i < _NUM_PROD_; i++)
	  _DERIV_ID_(i + _NUM_REACT_) = deriv_ids[_PROD_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  for (int i_ind = 0; i_ind < _NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < _NUM_REACT_; i_dep++) {
      _JAC_ID_(i_jac++) = jac_ids[_REACT_(i_dep)][_REACT_(i_ind)];
    }
    for (int i_dep = 0; i_dep < _NUM_PROD_; i_dep++) {
      _JAC_ID_(i_jac++) = jac_ids[_PROD_(i_dep)][_REACT_(i_ind)];
    }
  }
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * For CMAQ_H2O2 reaction this only involves recalculating the rate 
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_CMAQ_H2O2_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Calculate the rate constant in (#/cc)
  // k = k1 + [M]*k2
  realtype conv = _CONV_ * _PRESSURE_PA_ / _TEMPERATURE_K_;
  _RATE_CONSTANT_ = (_k1_A_ 
	  * (_k1_C_==ZERO ? ONE : SUNRexp(_k1_C_/_TEMPERATURE_K_))
	  * (_k1_B_==ZERO ? ONE : SUNRpowerR(_TEMPERATURE_K_/((realtype)300.0), _k1_B_))
	  + _k2_A_ // [M] is included in _k2_A_
	  * (_k2_C_==ZERO ? ONE : SUNRexp(_k2_C_/_TEMPERATURE_K_))
	  * (_k2_B_==ZERO ? ONE : SUNRpowerR(_TEMPERATURE_K_/((realtype)300.0), _k2_B_))
	  * conv
	  ) * SUNRpowerI(conv, _NUM_REACT_-1);

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
void * rxn_CMAQ_H2O2_calc_deriv_contrib(double *state, double *deriv,
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
void * rxn_CMAQ_H2O2_calc_jac_contrib(realtype *state, realtype *J,
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
void * rxn_CMAQ_H2O2_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the CMAQ_H2O2 reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_CMAQ_H2O2_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  printf("\n\nCMAQ_H2O2 reaction\n");
  for (int i=0; i<_INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _k1_A_
#undef _k1_B_
#undef _k1_C_
#undef _k2_A_
#undef _k2_B_
#undef _k2_C_
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
