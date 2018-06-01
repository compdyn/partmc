/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Arrhenius reaction solver functions
 *
*/
/** \file
 * \brief Arrhenius reaction solver functions
*/
#ifdef PMC_USE_SUNDIALS

#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ int_data[0]
#define NUM_PROD_ int_data[1]
#define A_ float_data[0]
#define B_ float_data[1]
#define C_ float_data[2]
#define D_ float_data[3]
#define E_ float_data[4]
#define CONV_ float_data[5]
#define RATE_CONSTANT_ float_data[6]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 7
#define REACT_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x]-1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+(NUM_REACT_+2)*(NUM_REACT_+NUM_PROD_))
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+NUM_PROD_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      jac_struct[REACT_(i_dep)][REACT_(i_ind)] = true;
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      jac_struct[PROD_(i_dep)][REACT_(i_ind)] = true;
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}  

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Update the time derivative ids
  for (int i=0; i < NUM_REACT_; i++)
	  DERIV_ID_(i) = deriv_ids[REACT_(i)];
  for (int i=0; i < NUM_PROD_; i++)
	  DERIV_ID_(i + NUM_REACT_) = deriv_ids[PROD_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      JAC_ID_(i_jac++) = jac_ids[REACT_(i_dep)][REACT_(i_ind)];
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      JAC_ID_(i_jac++) = jac_ids[PROD_(i_dep)][REACT_(i_ind)];
    }
  }
  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Arrhenius reaction this only involves recalculating the rate 
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_update_env_state(realtype *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the rate constant in (#/cc)
  // k = A*exp(C/T) * (T/D)^B * (1+E*P)
  RATE_CONSTANT_ = A_ * SUNRexp(C_/TEMPERATURE_K_)
	  * (B_==ZERO ? ONE : SUNRpowerR(TEMPERATURE_K_/D_, B_))
	  * (E_==ZERO ? ONE : (ONE + E_*PRESSURE_PA_))
          * SUNRpowerI(CONV_*PRESSURE_PA_/TEMPERATURE_K_, NUM_REACT_-1);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for arrhenius reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_pre_calc(ModelData *model_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being computed (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_calc_deriv_contrib(ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the reaction rate
  realtype rate = RATE_CONSTANT_;
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++) rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue; 
      deriv[DERIV_ID_(i_dep_var)] -= rate;
    }
    for (int i_spec=0; i_spec<NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue; 
      deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_calc_jac_contrib(ModelData *model_data, realtype *J,
          void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the reaction rate
  realtype rate = RATE_CONSTANT_;
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++) rate *= state[REACT_(i_spec)];

  // Add contributions to the Jacobian
  if (rate!=ZERO) {
    int i_elem = 0;
    for (int i_ind=0; i_ind<NUM_REACT_; i_ind++) {
      for (int i_dep=0; i_dep<NUM_REACT_; i_dep++, i_elem++) {
	if (JAC_ID_(i_elem) < 0) continue;
	J[JAC_ID_(i_elem)] -= rate / state[REACT_(i_ind)];
      }
      for (int i_dep=0; i_dep<NUM_PROD_; i_dep++, i_elem++) {
	if (JAC_ID_(i_elem) < 0) continue;
	J[JAC_ID_(i_elem)] += YIELD_(i_dep) * rate / state[REACT_(i_ind)];
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Arrhenius reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_arrhenius_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nArrhenius reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef NUM_REACT_
#undef NUM_PROD_
#undef A_
#undef B_
#undef C_
#undef D_
#undef E_
#undef CONV_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef PROD_
#undef DERIV_ID_
#undef JAC_ID_
#undef YIELD_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_

#endif
