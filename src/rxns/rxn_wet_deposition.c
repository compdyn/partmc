/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Wet deposition reaction solver functions
 *
*/
/** \file
 * \brief Wet deposition reaction solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define RXN_ID_ (int_data[0])
#define NUM_SPEC_ (int_data[1])
#define BASE_RATE_ float_data[0]
#define SCALING_ float_data[1]
#define RATE_CONSTANT_ float_data[2]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 3
#define REACT_(s) (int_data[NUM_INT_PROP_+s]-1)
#define DERIV_ID_(s) int_data[NUM_INT_PROP_+NUM_SPEC_+s]
#define JAC_ID_(s) int_data[NUM_INT_PROP_+2*NUM_SPEC_+s]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+3*NUM_SPEC_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    jac_struct[REACT_(i_spec)][REACT_(i_spec)] = true;
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
void * rxn_wet_deposition_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {

    // Update the time derivative id
    DERIV_ID_(i_spec) = deriv_ids[REACT_(i_spec)];

    // Update the Jacobian id
    JAC_ID_(i_spec) = jac_ids[REACT_(i_spec)][REACT_(i_spec)];

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data
 *
 * Wet deposition reactions can have their base (pre-scaling) rate constants
 * updated from the host model based on the calculations of an external
 * module. The structure of the update data is:
 *
 *  - \b int rxn_id (Id of one or more wet deposition reactions set by the
 *       host model using the
 *       \c pmc_rxn_wet_deposition::rxn_wet_deposition_t::set_rxn_id
 *       function prior to initializing the solver.)
 *  - \b double rate_const (New pre-scaling rate constant.)
 *
 * \param update_data Pointer to the updated reaction data
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_update_data(void *update_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  int *rxn_id = (int*) update_data;
  double *base_rate = (double*) &(rxn_id[1]);

  // Set the base wet deposition rate constants for matching reactions
  if (*rxn_id==RXN_ID_ && RXN_ID_!=0)
          BASE_RATE_ = (double) *base_rate;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update reaction data for new environmental conditions
 *
 * For wet deposition reactions this only involves recalculating the rate
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_update_env_state(double *env_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the rate constant in (1/s)
  RATE_CONSTANT_ = SCALING_ * BASE_RATE_;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for wet deposition reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_pre_calc(ModelData *model_data, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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
#ifdef PMC_USE_SUNDIALS
void * rxn_wet_deposition_calc_deriv_contrib(double *state, ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  //realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Add contributions to the time derivative
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (DERIV_ID_(i_spec) >= 0 )
        deriv[DERIV_ID_(i_spec)] -= RATE_CONSTANT_ * state[REACT_(i_spec)];
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void * rxn_wet_deposition_calc_jac_contrib(double *state, ModelData *model_data, realtype *J,
          void *rxn_data, double time_step)
{
  //realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Add contributions to the Jacobian
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (JAC_ID_(i_spec) >= 0) J[JAC_ID_(i_spec)] -= RATE_CONSTANT_;
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_wet_deposition_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nWet deposition reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Create update data for new wet deposition rates
 *
 * \return Pointer to a new rate update data object
 */
void * rxn_wet_deposition_create_rate_update_data()
{
  int *update_data = (int*) malloc(sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for wet deposition update data\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param rxn_id Id of wet deposition reactions to update
 * \param base_rate New pre-scaling wet deposition rate
 */
void rxn_wet_deposition_set_rate_update_data(void *update_data, int rxn_id,
          double base_rate)
{
  int *new_rxn_id = (int*) update_data;
  double *new_base_rate = (double*) &(new_rxn_id[1]);
  *new_rxn_id = rxn_id;
  *new_base_rate = base_rate;
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef RXN_ID_
#undef NUM_SPEC_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef DERIV_ID_
#undef JAC_ID_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
