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
#include "../rxn_solver.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472

#define DELTA_H_ float_data[0]
#define DELTA_S_ float_data[1]
#define DIFF_COEFF_ float_data[2]
#define PRE_C_AVG_ float_data[3]
#define A_ float_data[4]
#define C_ float_data[5]
#define C_AVG_ALPHA_ float_data[6]
#define EQUIL_CONST_ float_data[7]
#define CONV_ float_data[8]
#define MW_ float_data[9]
#define UGM3_TO_PPM_ float_data[10]
#define SMALL_NUMBER_ float_data[11]
#define NUM_AERO_PHASE_ int_data[0]
#define GAS_SPEC_ (int_data[1]-1)
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 12
#define AERO_SPEC_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define AERO_WATER_(x) (int_data[NUM_INT_PROP_ + NUM_AERO_PHASE_ + x]-1)
#define AERO_PHASE_ID_(x) (int_data[NUM_INT_PROP_ + 2*NUM_AERO_PHASE_ + x]-1)
#define AERO_REP_ID_(x) (int_data[NUM_INT_PROP_ + 3*(NUM_AERO_PHASE_) + x]-1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + 4*(NUM_AERO_PHASE_) + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+2+(10*NUM_AERO_PHASE_))
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_)

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct 2D array of flags indicating potentially non-zero 
 *                   Jacobian elements
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_get_used_jac_elem(void *rxn_data,
          bool **jac_struct)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  jac_struct[GAS_SPEC_][GAS_SPEC_] = true;
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
    jac_struct[AERO_SPEC_(i_aero_phase)][GAS_SPEC_] = true;
    jac_struct[GAS_SPEC_][AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[AERO_SPEC_(i_aero_phase)][AERO_SPEC_(i_aero_phase)] = true;
    jac_struct[GAS_SPEC_][AERO_WATER_(i_aero_phase)] = true;
    jac_struct[AERO_SPEC_(i_aero_phase)][AERO_WATER_(i_aero_phase)] = true;
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
void * rxn_HL_phase_transfer_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Update the time derivative ids
  DERIV_ID_(0) = deriv_ids[GAS_SPEC_];
  for (int i=0; i < NUM_AERO_PHASE_; i++)
	  DERIV_ID_(i + 1) = deriv_ids[AERO_SPEC_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][GAS_SPEC_];  
  for (int i_aero_phase = 0; i_aero_phase < NUM_AERO_PHASE_; i_aero_phase++) {
      JAC_ID_(i_jac++) = jac_ids[AERO_SPEC_(i_aero_phase)][GAS_SPEC_];
      JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][AERO_SPEC_(i_aero_phase)];
      JAC_ID_(i_jac++) =
              jac_ids[AERO_SPEC_(i_aero_phase)][AERO_SPEC_(i_aero_phase)];
      JAC_ID_(i_jac++) = jac_ids[GAS_SPEC_][AERO_WATER_(i_aero_phase)];
      JAC_ID_(i_jac++) = 
              jac_ids[AERO_SPEC_(i_aero_phase)][AERO_WATER_(i_aero_phase)];
    }

  // Calculate a small number based on the integration tolerances to use
  // during solving. TODO find a better place to do this
  realtype *abs_tol = model_data->abs_tol;
  SMALL_NUMBER_ = ( abs_tol[GAS_SPEC_] > abs_tol[AERO_SPEC_(0)] ?
                    abs_tol[AERO_SPEC_(0)] / 10.0 : abs_tol[GAS_SPEC_] / 10.0 );

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
void * rxn_HL_phase_transfer_update_env_state(double *env_data,
          void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the mass accomodation coefficient if the N* parameter
  // was provided, otherwise set it to 1.0
  double mass_acc = 1.0;
  if (DELTA_H_!=0.0 || DELTA_S_!=0.0) {
    double del_G = DELTA_H_ - TEMPERATURE_K_ * DELTA_S_; 
    mass_acc = exp(-del_G/(UNIV_GAS_CONST_ * TEMPERATURE_K_));
    mass_acc = mass_acc / (1.0 + mass_acc);
  }

  // Save c_rms * mass_acc for use in mass transfer rate calc
  C_AVG_ALPHA_ = PRE_C_AVG_ * sqrt(TEMPERATURE_K_) * mass_acc;

  // Calculate the Henry's Law equilibrium rate constant in units of
  // (ug_x/ug_H2O/ppm) where x is the aerosol-phase species. (A was saved in
  // units of M/ppm.)
  if (C_==0.0) {
    EQUIL_CONST_ = A_ * MW_;
  } else {
    EQUIL_CONST_ = A_ * exp(C_ * (1.0/TEMPERATURE_K_ - 1.0/298.0)) * MW_;
  }

  // Calculate the conversion from ug/m^3 -> ppm
  UGM3_TO_PPM_ = CONV_ * TEMPERATURE_K_ / PRESSURE_PA_;

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for phase_transfer reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_pre_calc(ModelData *model_data, void *rxn_data)
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
void * rxn_HL_phase_transfer_calc_deriv_contrib(ModelData *model_data,
          realtype *deriv, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)
   
    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number concentration
                                                // (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase 
    // mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] / 
        (aero_conc_type==0?1.0:number_conc) < SMALL_NUMBER_) continue;

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) + 
              4.0*radius/(3.0*C_AVG_ALPHA_));
  
    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (UGM3_TO_PPM_ *
	    EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

    // Calculate gas-phase condensation rate (ppm/s)
    // (Slow down the condensation as gas-phase concentrations approach zero to
    //  help out the solver.)
    cond_rate *= ( state[GAS_SPEC_] > SMALL_NUMBER_ ?
                   state[GAS_SPEC_] - SMALL_NUMBER_ : 0.0 );

    // Calculate aerosol-phase evaporation rate (ug/m^3/s)
    // (Slow down the condensation as gas-phase concentrations approach zero to
    //  help out the solver.)
    evap_rate *= ( state[AERO_SPEC_(i_phase)] > SMALL_NUMBER_ ?
                   state[AERO_SPEC_(i_phase)] - SMALL_NUMBER_ : 0.0 );

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (DERIV_ID_(0)>=0) {
      if (aero_conc_type==0) {
        // Scale the changes to the gas-phase by the number of particles for
        // per-particle aerosol concentrations
        deriv[DERIV_ID_(0)] += number_conc *
                (evap_rate * UGM3_TO_PPM_ - cond_rate);
      } else {
        // No scaling for aerosol concentrations with total mass per aerosol
        // phase
        deriv[DERIV_ID_(0)] += evap_rate * UGM3_TO_PPM_ - cond_rate;
      }
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (DERIV_ID_(1+i_phase)>=0) deriv[DERIV_ID_(1+i_phase)] += 
            cond_rate / UGM3_TO_PPM_ - evap_rate;

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
void * rxn_HL_phase_transfer_calc_jac_contrib(ModelData *model_data,
          realtype *J, void *rxn_data, double time_step)
{
  realtype *state = model_data->state;
  realtype *env_data = model_data->env;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate derivative contributions for each aerosol phase
  for (int i_phase=0; i_phase<NUM_AERO_PHASE_; i_phase++) {

    // Get the particle effective radius (m)
    realtype radius;
    aero_rep_get_effective_radius(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &radius);			// particle effective radius (m)
   
    // Get the particle number concentration (#/cc)
    realtype number_conc;
    aero_rep_get_number_conc(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase),	// aerosol phase index
		  &number_conc);		// particle number concentration
                                                // (#/cc)

    // Check the aerosol concentration type (per-particle or total per-phase mass)
    int aero_conc_type = aero_rep_get_aero_conc_type(
		  model_data,			// model data
		  AERO_REP_ID_(i_phase),	// aerosol representation index
		  AERO_PHASE_ID_(i_phase));	// aerosol phase index

    // If no aerosol water is present, no transfer occurs
    if (state[AERO_WATER_(i_phase)] / 
        (aero_conc_type==0?1.0:number_conc) < SMALL_NUMBER_) continue;

    // If the radius or number concentration are zero, no transfer occurs
    if (radius <= ZERO || number_conc <= ZERO) continue;

    // Calculate the rate constant for diffusion limited mass transfer to the
    // aerosol phase (1/s)
    realtype cond_rate = 1.0/(radius*radius/(3.0*DIFF_COEFF_) + 
              4.0*radius/(3.0*C_AVG_ALPHA_));
  
    // Calculate the evaporation rate constant (1/s)
    realtype evap_rate = cond_rate / (UGM3_TO_PPM_ *
	    EQUIL_CONST_ * state[AERO_WATER_(i_phase)]);

    // Adjust rates for lower threshholds (see derivative calculation)
    cond_rate *= ( state[GAS_SPEC_] > SMALL_NUMBER_ ? 1.0 : 0.0 );
    evap_rate *= ( state[AERO_SPEC_(i_phase)] > SMALL_NUMBER_ ? 1.0 : 0.0 );

    // Change in the gas-phase is evaporation - condensation (ppm/s)
    if (aero_conc_type==0) {
      // Scale the changes to the gas-phase by the number of particles for
      // per-particle aerosol concentrations
      if (JAC_ID_(1+i_phase*5+1)>=0) 
	      J[JAC_ID_(1+i_phase*5+1)] += number_conc * evap_rate * 
                      UGM3_TO_PPM_;
      if (JAC_ID_(1+i_phase*5+3)>=0) 
	      J[JAC_ID_(1+i_phase*5+3)] += - number_conc * evap_rate * 
                      UGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)] / 
                      state[AERO_WATER_(i_phase)];
      if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= number_conc * cond_rate;
    } else {
      // No scaling for aerosol concentrations with total mass per aerosol phase
      if (JAC_ID_(1+i_phase*5+1)>=0) J[JAC_ID_(1+i_phase*5+1)] +=
              evap_rate * UGM3_TO_PPM_;
      if (JAC_ID_(1+i_phase*5+3)>=0) J[JAC_ID_(1+i_phase*5+3)] -= 
              evap_rate * UGM3_TO_PPM_ * state[AERO_SPEC_(i_phase)] /
              state[AERO_WATER_(i_phase)];
      if (JAC_ID_(0)>=0) J[JAC_ID_(0)] -= cond_rate;
    }

    // Change in the aerosol-phase species is condensation - evaporation
    // (ug/m^3/s)
    if (JAC_ID_(1+i_phase*5)>=0) J[JAC_ID_(1+i_phase*5)] +=
            cond_rate / UGM3_TO_PPM_;
    if (JAC_ID_(1+i_phase*5+2)>=0) J[JAC_ID_(1+i_phase*5+2)] -= evap_rate; 
    if (JAC_ID_(1+i_phase*5+4)>=0) J[JAC_ID_(1+i_phase*5+4)] -= - evap_rate * 
	    state[AERO_SPEC_(i_phase)] / state[AERO_WATER_(i_phase)];
  
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Advance the reaction data pointer to the next reaction
 * 
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_skip(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_HL_phase_transfer_print(void *rxn_data)
{
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nPhase Transfer reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++) 
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
 
  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef UNIV_GAS_CONST_

#undef DELTA_H_
#undef DELTA_S_
#undef DIFF_COEFF_
#undef PRE_C_AVG_
#undef A_
#undef C_
#undef C_AVG_ALPHA_
#undef EQUIL_CONST_
#undef CONV_
#undef MW_
#undef UGM3_TO_PPM_
#undef NUM_AERO_PHASE_
#undef GAS_SPEC_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef AERO_SPEC_
#undef AERO_WATER_
#undef AERO_PHASE_ID_
#undef AERO_REP_ID_
#undef DERIV_ID_
#undef JAC_ID_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
