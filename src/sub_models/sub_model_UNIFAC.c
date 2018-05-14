/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * UNIFAC activity coefficient calculation
 *
 */
/** \file
 * \brief UNIFAC activity coefficient calculation
 */
#ifdef PMC_USE_SUNDIALS

#include "../sub_model_solver.h"

// TODO Lookup environmental indicies during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define _NUM_SPEC_ int_data[0]
#define _NUM_INT_PROP_ 1
#define _NUM_FLOAT_PROP_ 0
#define _INT_DATA_SIZE_ (_NUM_INT_PROP_)
#define _FLOAT_DATA_SIZE_ (_NUM_FLOAT_PROP_)

// Update types (These must match values in sub_model_UNIFAC.F90)
// (none right now)

/** \brief Get the Jacobian elements used for a particular row of the matrix
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param jac_row Array of flags indicating whether an element in the rown is used
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_get_used_jac_elem(void *sub_model_data, bool *jac_row)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update stored ids for elements used within a row of the Jacobian matrix
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param jac_row An array of new ids for one row of the Jacobian matrix
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_update_ids(void *sub_model_data, int *jac_row)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update sub-model data for new environmental conditions
 *
 * \param env_data Pointer to the environmental state array
 * \param sub_model_data Pointer to the sub-model data
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_update_env_state(void *sub_model_data, realtype *env_data)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Perform the sub-model calculations for the current model state
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param model_data Pointer to the model data including the current state and
 *                   environmental conditions
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_calculate(void *sub_model_data, ModelData *model_data)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * Derivatives are assumed to be of the form \f$ frac{dy}{dt} = A*S \f$, where
 * \f$A\f$ is the value passed to this function as \b base_val and \f$S\f$ is
 * the sub-model parameter used in the calculation. The row of the Jacobian
 * should correspond to \f$frac{dy'}{dx}\f$, where for each element \f$x\f$,
 * on the row, this function will add \f$A*frac{dS}{dx}\f$.
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param base_val The derivative 
 */
void * sub_model_UNIFAC_add_jac_contrib(void *sub_model_data,
         realtype base_val, realtype *jac_row)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Print the sub model data
 *
 * \param sub_model_data Pointer to the sub model data
 * \return The sub_model_data pointer advanced by the size of the sub-model
 */
void * sub_model_UNIFAC_print(void *sub_model_data)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _NUM_SPEC_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_
#endif
