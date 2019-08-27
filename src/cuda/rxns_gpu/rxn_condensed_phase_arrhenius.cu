/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Condensed Phase Arrhenius reaction solver functions
 *
*/
/** \file
 * \brief Condensed Phase Arrhenius reaction solver functions
*/
extern "C"{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns_gpu.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Small number
#define SMALL_NUMBER_ 1.0e-30

#define NUM_REACT_ (int_data[0*n_rxn])
#define NUM_PROD_ (int_data[1*n_rxn])
#define NUM_AERO_PHASE_ (int_data[2*n_rxn])
#define A_ (float_data[0*n_rxn])
#define B_ (float_data[1*n_rxn])
#define C_ (float_data[2*n_rxn])
#define D_ (float_data[3*n_rxn])
#define E_ (float_data[4*n_rxn])
#define RATE_CONSTANT_ (float_data[5*n_rxn])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 6
#define REACT_(x) (int_data[(NUM_INT_PROP_ + x)*n_rxn]-1)
#define PROD_(x) (int_data[(NUM_INT_PROP_+NUM_REACT_*NUM_AERO_PHASE_+x)*n_rxn]-1)
#define WATER_(x) (int_data[(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_)*NUM_AERO_PHASE_+x)*n_rxn]-1)
#define DERIV_ID_(x) (int_data[(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+1)*NUM_AERO_PHASE_+x)*n_rxn])
#define JAC_ID_(x) (int_data[(NUM_INT_PROP_+(2*(NUM_REACT_+NUM_PROD_)+1)*NUM_AERO_PHASE_+x)*n_rxn])
#define YIELD_(x) (float_data[(NUM_FLOAT_PROP_+x)*n_rxn])
#define UGM3_TO_MOLM3_(x) (float_data[(NUM_FLOAT_PROP_+NUM_PROD_+x)*n_rxn])
#define INT_DATA_SIZE_ (NUM_INT_PROP_+((NUM_REACT_+NUM_PROD_)*(NUM_REACT_+3)+1)*NUM_AERO_PHASE_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+2*NUM_PROD_+NUM_REACT_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_condensed_phase_arrhenius_get_used_jac_elem(void *rxn_data,
          bool **jac_struct)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Loop over all the instances of the specified phase
  for (int i_phase = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {

    // Add dependence on reactants for reactants and products
    for (int i_react_ind = i_phase*NUM_REACT_;
              i_react_ind < (i_phase+1)*NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][REACT_(i_react_ind)] = true;
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][REACT_(i_react_ind)] = true;
    }

    // Add dependence on aerosol-phase water for reactants and products in
    // aqueous reactions
    if (WATER_(i_phase)>=0) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        jac_struct[REACT_(i_react_dep)][WATER_(i_phase)] = true;
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        jac_struct[PROD_(i_prod_dep)][WATER_(i_phase)] = true;
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
void * rxn_gpu_condensed_phase_arrhenius_update_ids(ModelData *model_data,
          int *deriv_ids, int **jac_ids, void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Update the time derivative ids
  for (int i_phase = 0, i_deriv = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {
    for (int i_react = 0; i_react < NUM_REACT_; i_react++)
      DERIV_ID_(i_deriv++) = deriv_ids[REACT_(i_phase*NUM_REACT_+i_react)];
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++)
      DERIV_ID_(i_deriv++) = deriv_ids[PROD_(i_phase*NUM_PROD_+i_prod)];
  }

  // Update the Jacobian ids
  for (int i_phase = 0, i_jac = 0; i_phase < NUM_AERO_PHASE_; i_phase++) {

    // Add dependence on reactants for reactants and products
    for (int i_react_ind = i_phase*NUM_REACT_;
              i_react_ind < (i_phase+1)*NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = i_phase*NUM_REACT_;
                i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][REACT_(i_react_ind)];
      for (int i_prod_dep = i_phase*NUM_PROD_;
                i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][REACT_(i_react_ind)];
    }

    // Add dependence on aerosol-phase water for reactants and products
    for (int i_react_dep = i_phase*NUM_REACT_;
              i_react_dep < (i_phase+1)*NUM_REACT_; i_react_dep++)
      if (WATER_(i_phase)>=0) {
        JAC_ID_(i_jac++) = jac_ids[REACT_(i_react_dep)][WATER_(i_phase)];
      } else {
        JAC_ID_(i_jac++) = -1;
      }
    for (int i_prod_dep = i_phase*NUM_PROD_;
              i_prod_dep < (i_phase+1)*NUM_PROD_; i_prod_dep++)
      if (WATER_(i_phase)>=0) {
        JAC_ID_(i_jac++) = jac_ids[PROD_(i_prod_dep)][WATER_(i_phase)];
      } else {
        JAC_ID_(i_jac++) = -1;
      }

  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

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
__device__ void rxn_gpu_condensed_phase_arrhenius_update_env_state(double *rate_constants,
   int n_rxn2,double *double_pointer_gpu, double *env_data,
          void *rxn_data)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the rate constant in (M or mol/m3)
  // k = A*exp(C/T) * (T/D)^B * (1+E*P)
  RATE_CONSTANT_ = A_ * exp(C_/TEMPERATURE_K_)
          * (B_==0.0 ? 1.0 : pow(TEMPERATURE_K_/D_, B_))
          * (E_==0.0 ? 1.0 : (1.0 + E_*PRESSURE_PA_));

  rate_constants[0] = RATE_CONSTANT_;

}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for condensed_phase_arrhenius reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_condensed_phase_arrhenius_pre_calc(ModelData *model_data,
          void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
__device__ void rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    double unit_conv = 1.0;
    if (WATER_(i_phase)>=0) {
      unit_conv = state[WATER_(i_phase)] * 1.0e-9; // convert from ug/m3->L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction
      // occurs
      if (unit_conv < SMALL_NUMBER_) {
        i_deriv += NUM_REACT_ + NUM_PROD_;
        continue;
      }
      unit_conv = 1.0/unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    //double rate = RATE_CONSTANT_;
    double rate = rate_constants[0];
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              UGM3_TO_MOLM3_(i_react) * unit_conv;
    }

    // Reactant change
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      //deriv[DERIV_ID_(i_deriv++)] -= rate /
	 //     (UGM3_TO_MOLM3_(i_react) * unit_conv);
      atomicAdd((double*)&(deriv[DERIV_ID_(i_deriv++)]), -(rate /
      (UGM3_TO_MOLM3_(i_react) * unit_conv)));
    }

    // Products change
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      //deriv[DERIV_ID_(i_deriv++)] += rate * YIELD_(i_prod) /
      //  (UGM3_TO_MOLM3_(NUM_REACT_+i_prod) * unit_conv);
      atomicAdd((double*)&(deriv[DERIV_ID_(i_deriv++)]),rate * YIELD_(i_prod) /
	      (UGM3_TO_MOLM3_(NUM_REACT_+i_prod) * unit_conv));
    }

  }

}
#endif


/** \brief Calculate contributions to the time derivative f(t,y) from this
 * reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
void rxn_cpu_condensed_phase_arrhenius_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0, i_deriv = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    double unit_conv = 1.0;
    if (WATER_(i_phase)>=0) {
      unit_conv = state[WATER_(i_phase)] * 1.0e-9; // convert from ug/m3->L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction
      // occurs
      if (unit_conv < SMALL_NUMBER_) {
        i_deriv += NUM_REACT_ + NUM_PROD_;
        continue;
      }
      unit_conv = 1.0/unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    //double rate = RATE_CONSTANT_;
  double rate = rate_constants[0];
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              UGM3_TO_MOLM3_(i_react) * unit_conv;
    }

    // Reactant change
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] -= rate /
	      (UGM3_TO_MOLM3_(i_react) * unit_conv);
    }

    // Products change
    for (int i_prod = 0; i_prod < NUM_PROD_; i_prod++) {
      if (DERIV_ID_(i_deriv)<0) {i_deriv++; continue;}
      deriv[DERIV_ID_(i_deriv++)] += rate * YIELD_(i_prod) /
        (UGM3_TO_MOLM3_(NUM_REACT_+i_prod) * unit_conv);
    }

  }

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step of the itegrator (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
__device__ void rxn_gpu_condensed_phase_arrhenius_calc_jac_contrib(double *rate_constants, double *state,
          double *J, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate Jacobian contributions for each aerosol phase
  for (int i_phase=0, i_jac = 0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // If this is an aqueous reaction, get the unit conversion from mol/m3 -> M
    double unit_conv = 1.0;
    if (WATER_(i_phase)>=0) {
      unit_conv = state[WATER_(i_phase)] * 1.0e-9; // convert from ug/m3->L/m3

      // For aqueous reactions, if no aerosol water is present, no reaction
      // occurs
      if (unit_conv < SMALL_NUMBER_) {
        i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + 1);
        continue;
      }
      unit_conv = 1.0/unit_conv;
    }

    // Calculate the reaction rate rate (M/s or mol/m3/s)
    //double rate = RATE_CONSTANT_;
  double rate = rate_constants[0];
    for (int i_react = 0; i_react < NUM_REACT_; i_react++) {
      rate *= state[REACT_(i_phase*NUM_REACT_+i_react)] *
              UGM3_TO_MOLM3_(i_react) * unit_conv;
    }

    // No Jac contributions to add if the rate is zero
    if (rate==0.0) {
      i_jac += (NUM_REACT_ + NUM_PROD_) * (NUM_REACT_ + 1);
      continue;
    }

    // Add dependence on reactants for reactants and products
    for (int i_react_ind = 0; i_react_ind < NUM_REACT_; i_react_ind++) {
      for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
	if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
        atomicAdd((double*)&(J[JAC_ID_(i_jac++)]), -rate /
          state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
          (UGM3_TO_MOLM3_(i_react_dep) * unit_conv));
      }
      for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
	if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
        atomicAdd((double*)&(J[JAC_ID_(i_jac++)]),
          rate * YIELD_(i_prod_dep) /
          state[REACT_(i_phase*NUM_REACT_+i_react_ind)] /
          (UGM3_TO_MOLM3_(NUM_REACT_+i_prod_dep) * unit_conv));
      }
    }

    // Add dependence on aerosol-phase water for reactants and products in
    // aqueous reactions
    for (int i_react_dep = 0; i_react_dep < NUM_REACT_; i_react_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd((double*)&(J[JAC_ID_(i_jac++)]),
          (NUM_REACT_-1) * rate / state[WATER_(i_phase)] /
          (UGM3_TO_MOLM3_(i_react_dep) * unit_conv));
    }
    for (int i_prod_dep = 0; i_prod_dep < NUM_PROD_; i_prod_dep++) {
      if (JAC_ID_(i_jac)<0) {i_jac++; continue;}
      atomicAdd((double*)&(J[JAC_ID_(i_jac++)]),
          -(NUM_REACT_-1) * rate * YIELD_(i_prod_dep) /
          state[WATER_(i_phase)] /
          (UGM3_TO_MOLM3_(NUM_REACT_+i_prod_dep) * unit_conv));
    }

  }

}
#endif

/** \brief Retrieve Int data size
 *
 * \param rxn_data Pointer to the reaction data
 * \return The data size of int array
 */
void * rxn_gpu_condensed_phase_arrhenius_get_float_pointer(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) float_data;
}

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_condensed_phase_arrhenius_skip(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Condensed Phase Arrhenius reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_condensed_phase_arrhenius_print(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nCondensed Phase Arrhenius reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef SMALL_NUMBER_

#undef NUM_REACT_
#undef NUM_PROD_
#undef NUM_AERO_PHASE_
#undef A_
#undef B_
#undef C_
#undef D_
#undef E_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef PROD_
#undef WATER_
#undef DERIV_ID_
#undef JAC_ID_
#undef YIELD_
#undef UGM3_TO_MOLM3_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
}