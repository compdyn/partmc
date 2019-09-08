/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * ZSR Aerosol Water sub model solver functions
 *
*/
/** \file
 * \brief ZSR Aerosol Water sub model solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../sub_models.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define SMALL_NUMBER_ 1.0e-30

#define ACT_TYPE_JACOBSON 1
#define ACT_TYPE_EQSAM 2

#define NUM_PHASE_ (int_data[0])
#define GAS_WATER_ID_ (int_data[1]-1)
#define NUM_ION_PAIR_ (int_data[2])
#define INT_DATA_SIZE_ (int_data[3])
#define FLOAT_DATA_SIZE_ (int_data[4])
#define PPM_TO_RH_ (sub_model_env_data[0])
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 0
#define NUM_ENV_PARAM_ 1
#define PHASE_ID_(p) (int_data[NUM_INT_PROP_+p]-1)
#define PAIR_INT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PAIR_FLOAT_PARAM_LOC_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIR_+x]-1)
#define TYPE_(x) (int_data[PAIR_INT_PARAM_LOC_(x)])
#define JACOB_NUM_CATION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define JACOB_NUM_ANION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+2])
#define JACOB_CATION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+3])
#define JACOB_ANION_ID_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+4])
#define JACOB_NUM_Y_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+5])
#define JACOB_GAS_WATER_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+p]
#define JACOB_CATION_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+NUM_PHASE_+p]
#define JACOB_ANION_JAC_ID_(p,x) int_data[PAIR_INT_PARAM_LOC_(x)+6+2*NUM_PHASE_+p]
#define EQSAM_NUM_ION_(x) (int_data[PAIR_INT_PARAM_LOC_(x)+1])
#define EQSAM_GAS_WATER_JAC_ID_(p,x) (int_data[PAIR_INT_PARAM_LOC_(x)+2+p])
#define EQSAM_ION_ID_(x,y) (int_data[PAIR_INT_PARAM_LOC_(x)+2+NUM_PHASE_+y])
#define EQSAM_ION_JAC_ID_(p,x,y) int_data[PAIR_INT_PARAM_LOC_(x)+2+NUM_PHASE_+EQSAM_NUM_ION_(x)+y*NUM_PHASE_+p]
#define JACOB_low_RH_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define JACOB_CATION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define JACOB_ANION_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define JACOB_Y_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])
#define EQSAM_NW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)])
#define EQSAM_ZW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+1])
#define EQSAM_ION_PAIR_MW_(x) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+2])
#define EQSAM_ION_MW_(x,y) (float_data[PAIR_FLOAT_PARAM_LOC_(x)+3+y])

// Update types (These must match values in sub_model_UNIFAC.F90)
// (none right now)

/** \brief Flag Jacobian elements used by this sub model
 *
 * ZSR aerosol water sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_struct A matrix of flags for needed Jac elements
 */
void sub_model_ZSR_aerosol_water_get_used_jac_elem(int *sub_model_int_data,
    double *sub_model_float_data, bool **jac_struct)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Loop through the dependent species - aerosol water and set all Jacobian
  // elements for each
  for (int i_phase=0; i_phase < NUM_PHASE_; ++i_phase) {

    // Flag the gas-phase water species
    jac_struct[PHASE_ID_(i_phase)][GAS_WATER_ID_] = true;

    // Flag elements for each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; ++i_ion_pair) {

      // Flag aerosol elements by calculation type
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Flag the anion and cation Jacobian elements
          jac_struct[PHASE_ID_(i_phase)]
                    [PHASE_ID_(i_phase)+JACOB_CATION_ID_(i_ion_pair)] = true;
          jac_struct[PHASE_ID_(i_phase)]
                    [PHASE_ID_(i_phase)+JACOB_ANION_ID_(i_ion_pair)] = true;
          break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Flag the ion Jacobian elements
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
            jac_struct[PHASE_ID_(i_phase)]
                      [PHASE_ID_(i_phase)+
                          EQSAM_ION_ID_(i_ion_pair,i_ion)] = true;
          }
          break;
      }
    }
  }
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * ZSR aerosol water sub models are assumed to be at equilibrium
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param deriv_ids Indices for state array variables on the solver state array
 * \param jac_ids Indices for Jacobian elements in the sparse data array
 */
void sub_model_ZSR_aerosol_water_update_ids(int *sub_model_int_data,
    double *sub_model_float_data, int *deriv_ids, int **jac_ids)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Loop through the dependent species - aerosol water and set all Jacobian
  // elements for each
  for (int i_phase=0; i_phase < NUM_PHASE_; ++i_phase) {

    // Flag elements for each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; ++i_ion_pair) {

      // Flag aerosol elements by calculation type
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Save the gas-phase water species
          JACOB_GAS_WATER_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)][GAS_WATER_ID_];

          // Save the cation and anion Jacobian elements
          JACOB_CATION_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)]
                     [PHASE_ID_(i_phase)+JACOB_CATION_ID_(i_ion_pair)];
          JACOB_ANION_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)]
                     [PHASE_ID_(i_phase)+JACOB_ANION_ID_(i_ion_pair)];
          break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Save the gas-phase water species
          EQSAM_GAS_WATER_JAC_ID_(i_phase,i_ion_pair) =
              jac_ids[PHASE_ID_(i_phase)][GAS_WATER_ID_];

          // Save the ion Jacobian elements
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
            EQSAM_ION_JAC_ID_(i_phase,i_ion_pair,i_ion) =
                jac_ids[PHASE_ID_(i_phase)]
                       [PHASE_ID_(i_phase)+EQSAM_ION_ID_(i_ion_pair,i_ion)];
          }
          break;
      }
    }
  }
}

/** \brief Update sub model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 */
void sub_model_ZSR_aerosol_water_update_env_state(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate PPM_TO_RH_
  // From MOSAIC code - reference to Seinfeld & Pandis page 181
  // TODO Figure out how to have consistent RH<->ppm conversions
  double t_steam = 373.15; 		// steam temperature (K)
  double a = 1.0 - t_steam/TEMPERATURE_K_;

  a = (((-0.1299*a - 0.6445)*a - 1.976)*a + 13.3185)*a;
  double water_vp = 101325.0 * exp(a); 		// (Pa)

  PPM_TO_RH_ = PRESSURE_PA_ / water_vp / 1.0e6;	// (1/ppm)
}

/** \brief Do pre-derivative calculations
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data, including the state array
 */
void sub_model_ZSR_aerosol_water_calculate(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data)
{
  double *state = model_data->grid_cell_state;
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Calculate the water activity---i.e., relative humidity (0-1)
  double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];

  // Calculate the total aerosol water for each instance of the aerosol phase
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {
    double *water = &(state[PHASE_ID_(i_phase)]);
    *water = 0.0;

    // Get the contribution from each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; i_ion_pair++) {

      double molality;
      double j_aw, e_aw;
      double conc;

      // Determine which type of activity calculation should be used
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Determine whether to use the minimum RH in the calculation
          j_aw = a_w>JACOB_low_RH_(i_ion_pair) ? a_w :
                  JACOB_low_RH_(i_ion_pair);

          // Calculate the molality of the pure binary ion pair solution
	  molality = 0.0;
          for (int i_order=0; i_order<JACOB_NUM_Y_(i_ion_pair); i_order++)
		  molality += JACOB_Y_(i_ion_pair, i_order) *
                          pow(j_aw,i_order);
          molality *= molality; // (mol/kg)

	  // Calculate the water associated with this ion pair
          double cation = state[PHASE_ID_(i_phase) +
                  JACOB_CATION_ID_(i_ion_pair)] /
		  JACOB_NUM_CATION_(i_ion_pair) /
                  JACOB_CATION_MW_(i_ion_pair) / 1000.0; // (umol/m3)
          double anion = state[PHASE_ID_(i_phase) +
                  JACOB_ANION_ID_(i_ion_pair)] /
		  JACOB_NUM_ANION_(i_ion_pair) /
                  JACOB_ANION_MW_(i_ion_pair) / 1000.0; // (umol/m3)
	  conc = (cation>anion ? anion : cation);
          conc = (conc>0.0 ? conc : 0.0);
          *water += conc / molality * 1000.0; // (ug/m3)

	  break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Keep the water activity within the range specified in EQSAM
          e_aw = a_w > 0.99 ? 0.99 : a_w;
          e_aw = e_aw < 0.001 ? 0.001 : e_aw;

	  // Calculate the molality of the ion pair
	  molality = (EQSAM_NW_(i_ion_pair) * 55.51 * 18.01 /
                    EQSAM_ION_PAIR_MW_(i_ion_pair) / 1000.0 * (1.0/e_aw-1.0));
	  molality = pow(molality, EQSAM_ZW_(i_ion_pair)); // (mol/kg)

	  // Calculate the water associated with this ion pair
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
	    conc = state[PHASE_ID_(i_phase)+EQSAM_ION_ID_(i_ion_pair,i_ion)];
            conc = (conc>0.0 ? conc : 0.0);
            *water += conc / EQSAM_ION_MW_(i_ion_pair,i_ion) /
                    molality; // (ug/m3);
	  }

	  break;
      }
    }
  }
}

/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current time step [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_ZSR_aerosol_water_get_jac_contrib(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data, realtype *J, double time_step)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *state    = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate the water activity---i.e., relative humidity (0-1)
  double a_w = PPM_TO_RH_ * state[GAS_WATER_ID_];
  double d_aw_d_wg = PPM_TO_RH_;

  // Calculate the total aerosol water for each instance of the aerosol phase
  for (int i_phase=0; i_phase<NUM_PHASE_; i_phase++) {

    // Get the contribution from each ion pair
    for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; i_ion_pair++) {

      double molality, d_molal_d_wg;
      double conc, d_conc_d_ion;
      double e_aw, d_eaw_d_wg;
      double j_aw, d_jaw_d_wg;

      // Determine which type of activity calculation should be used
      switch (TYPE_(i_ion_pair)) {

	// Jacobson et al. (1996)
	case ACT_TYPE_JACOBSON :

          // Determine whether to use the minimum RH in the calculation
          j_aw = a_w>JACOB_low_RH_(i_ion_pair) ? a_w :
                  JACOB_low_RH_(i_ion_pair);
          d_jaw_d_wg = a_w>JACOB_low_RH_(i_ion_pair) ? d_aw_d_wg : 0.0;

          // Calculate the molality of the pure binary ion pair solution
	  molality = JACOB_Y_(i_ion_pair, 0);
          d_molal_d_wg = 0.0;
          for (int i_order=1; i_order<JACOB_NUM_Y_(i_ion_pair); i_order++) {
            molality += JACOB_Y_(i_ion_pair, i_order) *
                        pow(j_aw,i_order);
            d_molal_d_wg += JACOB_Y_(i_ion_pair, i_order) * i_order *
                            pow(j_aw,(i_order-1));
          }
          d_molal_d_wg *= d_jaw_d_wg;

	  // Calculate the water associated with this ion pair
          double cation = state[PHASE_ID_(i_phase) +
                  JACOB_CATION_ID_(i_ion_pair)] /
		  JACOB_NUM_CATION_(i_ion_pair) /
                  JACOB_CATION_MW_(i_ion_pair) / 1000.0; // (umol/m3)
          double d_cation_d_C = 1.0 /
		  JACOB_NUM_CATION_(i_ion_pair) /
                  JACOB_CATION_MW_(i_ion_pair) / 1000.0;
          double anion = state[PHASE_ID_(i_phase) +
                  JACOB_ANION_ID_(i_ion_pair)] /
		  JACOB_NUM_ANION_(i_ion_pair) /
                  JACOB_ANION_MW_(i_ion_pair) / 1000.0; // (umol/m3)
          double d_anion_d_A = 1.0 /
		  JACOB_NUM_ANION_(i_ion_pair) /
                  JACOB_ANION_MW_(i_ion_pair) / 1000.0;

          // Determine the limiting ion and set Jacobian elements
          if (cation>anion) {

            // Anion-limited conditions
            if (anion>0.0) {
              J[JACOB_GAS_WATER_JAC_ID_(i_phase,i_ion_pair)] +=
                  -2.0 * anion / pow(molality, 3) * 1000.0 * d_molal_d_wg;
              J[JACOB_ANION_JAC_ID_(i_phase,i_ion_pair)] +=
                  1.0 / pow(molality, 2) * 1000.0 * d_anion_d_A;
            }
          } else {

            // Cation-limited conditions
            if (cation>0.0) {
              J[JACOB_GAS_WATER_JAC_ID_(i_phase,i_ion_pair)] +=
                  -2.0 * cation / pow(molality, 3) * 1000.0 * d_molal_d_wg;
              J[JACOB_CATION_JAC_ID_(i_phase,i_ion_pair)] +=
                  1.0 / pow(molality, 2) * 1000.0 * d_cation_d_C;
            }
          }

	  break;

	// EQSAM (Metger et al., 2002)
	case ACT_TYPE_EQSAM :

          // Keep the water activity within the range specified in EQSAM
          e_aw = a_w > 0.99 ? 0.99 : a_w;
          e_aw = e_aw < 0.001 ? 0.001 : e_aw;
          d_eaw_d_wg = a_w > 0.99 ? 0.0 : d_aw_d_wg;
          d_eaw_d_wg = a_w < 0.001 ? 0.0 : d_eaw_d_wg;

	  // Calculate the molality of the ion pair
	  molality = (EQSAM_NW_(i_ion_pair) * 55.51 * 18.01 /
                    EQSAM_ION_PAIR_MW_(i_ion_pair) / 1000.0 * (1.0/e_aw-1.0));
          d_molal_d_wg = -EQSAM_NW_(i_ion_pair) * 55.51 * 18.01 /
                    EQSAM_ION_PAIR_MW_(i_ion_pair) / 1000.0 / pow(e_aw,2) *
                    d_eaw_d_wg;
          d_molal_d_wg = EQSAM_ZW_(i_ion_pair) *
                         pow(molality, EQSAM_ZW_(i_ion_pair)-1.0) *
                         d_molal_d_wg;
	  molality = pow(molality, EQSAM_ZW_(i_ion_pair)); // (mol/kg)

	  // Calculate the Jacobian contributions
	  for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); i_ion++) {
	    conc = state[PHASE_ID_(i_phase)+EQSAM_ION_ID_(i_ion_pair,i_ion)];
            conc = (conc>0.0 ? conc : 0.0);
            d_conc_d_ion = (conc>0.0 ? 1.0 : 0.0);

            // Gas-phase water contribution
            J[EQSAM_GAS_WATER_JAC_ID_(i_phase,i_ion_pair)] +=
                -1.0 * conc / EQSAM_ION_MW_(i_ion_pair,i_ion) /
                pow(molality,2) * d_molal_d_wg;

            // Ion contribution
            J[EQSAM_ION_JAC_ID_(i_phase,i_ion_pair,i_ion)] +=
                d_conc_d_ion / EQSAM_ION_MW_(i_ion_pair,i_ion) /
                molality;
	  }

	  break;
      }
    }
  }
}
#endif

/** \brief Print the ZSR Aerosol Water sub model parameters
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 */
void sub_model_ZSR_aerosol_water_print(int *sub_model_int_data,
    double *sub_model_float_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  printf("\n\nZSR aerosol water sub model\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);
  printf("\nNumber of phases: %d", NUM_PHASE_);
  printf("\nGas-phase water state id: %d", GAS_WATER_ID_);
  printf("\nNumber of ion pairs: %d", NUM_ION_PAIR_);
  printf("\n*** Phase state ids (index of water in each phase) ***");
  for (int i_phase=0; i_phase<NUM_PHASE_; ++i_phase) {
    printf("\n  phase %d: %d", i_phase, PHASE_ID_(i_phase));
  }
  printf("\n*** Ion-Pair info ***");
  for (int i_ion_pair=0; i_ion_pair<NUM_ION_PAIR_; ++i_ion_pair) {
    printf("\n  ION PAIR %d", i_ion_pair);
    switch (TYPE_(i_ion_pair)) {
      case (ACT_TYPE_JACOBSON) :
        printf("\n    *** JACOBSON ***");
        printf("\n    low RH: %le", JACOB_low_RH_(i_ion_pair));
        printf("\n    number of cations: %d number of anions: %d",
               JACOB_NUM_CATION_(i_ion_pair), JACOB_NUM_ANION_(i_ion_pair));
        printf("\n    cation id: %d anion id: %d",
               JACOB_CATION_ID_(i_ion_pair), JACOB_ANION_ID_(i_ion_pair));
        printf("\n    cation MW: %le anion MW: %le",
               JACOB_CATION_MW_(i_ion_pair), JACOB_ANION_MW_(i_ion_pair));
        printf("\n    number of Y parameters: %d:", JACOB_NUM_Y_(i_ion_pair));
        for (int i_Y=0; i_Y<JACOB_NUM_Y_(i_ion_pair); ++i_Y)
          printf(" Y(%d)=%le", i_Y, JACOB_Y_(i_ion_pair,i_Y));
        for (int i_phase=0; i_phase<NUM_PHASE_; ++i_phase) {
          printf("\n    PHASE %d:", i_phase);
          printf(" gas-phase water Jac id: %d",
                 JACOB_GAS_WATER_JAC_ID_(i_phase,i_ion_pair));
          printf(" cation Jac id: %d",
                 JACOB_CATION_JAC_ID_(i_phase,i_ion_pair));
          printf(" anion Jac id: %d",
                 JACOB_ANION_JAC_ID_(i_phase,i_ion_pair));
        }
        break;
      case (ACT_TYPE_EQSAM) :
        printf("\n    *** EQSAM ***");
        printf("\n    NW: %le ZW: %le ion pair MW: %le",
               EQSAM_NW_(i_ion_pair), EQSAM_ZW_(i_ion_pair),
               EQSAM_ION_PAIR_MW_(i_ion_pair));
        printf("\n    number of ions: %d", EQSAM_NUM_ION_(i_ion_pair));
        printf("\n    IONS");
        for (int i_ion=0; i_ion<EQSAM_NUM_ION_(i_ion_pair); ++i_ion) {
          printf("\n      ion id: %d", EQSAM_ION_ID_(i_ion_pair, i_ion));
          for (int i_phase=0; i_phase<NUM_PHASE_; ++i_phase) {
            printf("\n        phase: %d gas-phase water Jac id: %d "
                   "ion Jac id: %d", i_phase,
                   EQSAM_GAS_WATER_JAC_ID_(i_phase,i_ion_pair),
                   EQSAM_ION_JAC_ID_(i_phase,i_ion_pair,i_ion));
          }
        }
        break;
      default :
        printf("\n !!! INVALID TYPE SPECIFIED: %d", TYPE_(i_ion_pair));
        break;
    }
  }
}
