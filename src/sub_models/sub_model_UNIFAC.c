/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * UNIFAC activity coefficient calculation
 *
 */
/** \file
 * \brief UNIFAC activity coefficient calculation
 *
 * For more info see the pmc_sub_model_UNIFAC module
 *
 * Equation references are to Marcolli and Peter, ACP 5(2), 1501-1527, 2005.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../sub_models.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_UNIQUE_PHASE_ (int_data[0])
#define NUM_GROUP_ (int_data[1])
#define TOTAL_INT_PROP_ (int_data[2])
#define TOTAL_FLOAT_PROP_ (int_data[3])
#define NUM_INT_PROP_ 4
#define NUM_FLOAT_PROP_ 0
#define PHASE_INT_LOC_(p) (int_data[NUM_INT_PROP_+p]-1)
#define PHASE_FLOAT_LOC_(p) (int_data[NUM_INT_PROP_+NUM_UNIQUE_PHASE_+p]-1)
#define NUM_PHASE_INSTANCE_(p) (int_data[PHASE_INT_LOC_(p)])
#define NUM_SPEC_(p) (int_data[PHASE_INT_LOC_(p)+1])
#define PHASE_INST_FLOAT_LOC_(p,c) (int_data[PHASE_INT_LOC_(p)+2+c]-1)
#define PHASE_INST_ID_(p,c) (int_data[PHASE_INT_LOC_(p)+2+NUM_PHASE_INSTANCE_(p)+c]-1)
#define SPEC_ID_(p,i) (int_data[PHASE_INT_LOC_(p)+2+2*NUM_PHASE_INSTANCE_(p)+i])
#define GAMMA_ID_(p,i) (int_data[PHASE_INT_LOC_(p)+2+2*NUM_PHASE_INSTANCE_(p)+NUM_SPEC_(p)+i])
#define JAC_ID_(p,c,i) int_data[PHASE_INT_LOC_(p)+2+2*NUM_PHASE_INSTANCE_(p)+(c+2)*NUM_SPEC_(p)+i]
#define V_IK_(p,i,k) (int_data[PHASE_INT_LOC_(p)+2+2*NUM_PHASE_INSTANCE_(p)+(k+2+NUM_PHASE_INSTANCE_(p))*NUM_SPEC_(p)+i])

#define Q_K_(k) (float_data[k])
#define R_K_(k) (float_data[NUM_GROUP_+k])
#define THETA_M_(m) (float_data[2*NUM_GROUP_+m])
#define A_MN_(m,n) (float_data[(m+3)*NUM_GROUP_+n])
#define PSI_MN_(m,n) (float_data[(m+3+NUM_GROUP_)*NUM_GROUP_+n])
#define R_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+i])
#define Q_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+NUM_SPEC_(p)+i])
#define L_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+2*NUM_SPEC_(p)+i])
#define MW_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+3*NUM_SPEC_(p)+i])
#define X_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+4*NUM_SPEC_(p)+i])
#define LN_GAMMA_IK_(p,i,k) (float_data[PHASE_FLOAT_LOC_(p)+i*NUM_GROUP_+5*NUM_SPEC_(p)+k])
#define GAMMA_I_(p,c,i) (float_data[PHASE_INST_FLOAT_LOC_(p,c)+i])

#define INT_DATA_SIZE_ (TOTAL_INT_PROP_)
#define FLOAT_DATA_SIZE_ (TOTAL_FLOAT_PROP_)

// Update types (These must match values in sub_model_UNIFAC.F90)
// (none right now)

/** \brief Get the Jacobian elements used for a particular row of the matrix
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param jac_struct A matrix of flags for needed Jac elements
 */
void sub_model_UNIFAC_get_used_jac_elem(int *sub_model_int_data,
    double *sub_model_float_data, bool **jac_struct)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  for (int i_phase=0; i_phase < NUM_UNIQUE_PHASE_; ++i_phase)
    for (int i_inst=0; i_inst < NUM_PHASE_INSTANCE_(i_phase); ++i_inst)
      for (int i_spec=0; i_spec < NUM_SPEC_(i_phase); ++i_spec)
        jac_struct[PHASE_INST_ID_(i_phase, i_inst)+GAMMA_ID_(i_phase, i_spec)]
                  [PHASE_INST_ID_(i_phase, i_inst)+
                   SPEC_ID_(i_phase, i_spec)] = true;
}

/** \brief Update stored ids for elements used within a row of the Jacobian matrix
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param deriv_ids Indices for state array variables on the solver state array
 * \param jac_ids Indices for Jacobian elements in the sparse data array
 */
void sub_model_UNIFAC_update_ids(int *sub_model_int_data,
    double *sub_model_float_data, int *deriv_ids, int **jac_ids)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  for (int i_phase=0; i_phase < NUM_UNIQUE_PHASE_; ++i_phase)
    for (int i_inst=0; i_inst < NUM_PHASE_INSTANCE_(i_phase); ++i_inst)
      for (int i_spec=0; i_spec < NUM_SPEC_(i_phase); ++i_spec)
        JAC_ID_(i_phase, i_inst, i_spec) =
          jac_ids[PHASE_INST_ID_(i_phase, i_inst)+GAMMA_ID_(i_phase, i_spec)]
                 [PHASE_INST_ID_(i_phase, i_inst)+SPEC_ID_(i_phase, i_spec)];
}

/** \brief Update sub-model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param env_data Pointer to the environmental state array
 */
void sub_model_UNIFAC_update_env_state(int *sub_model_int_data,
    double *sub_model_float_data, double *env_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Update the interaction parameters
  for (int m=0; m<NUM_GROUP_; m++)
    for (int n=0; n<NUM_GROUP_; n++)
      PSI_MN_(m,n) = exp(-A_MN_(m,n)/TEMPERATURE_K_);

  // Calculate the pure liquid residual acitivity ceofficient ln(GAMMA_k^(i))
  // terms. Eq. 7 & 8
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {

      // Calculate the sum (Q_n * X_n) in the denominator of Eq. 9 for the
      // pure liquid
      double sum_Qn_Xn = 0.0;
      for (int i_group=0; i_group<NUM_GROUP_; i_group++)
        sum_Qn_Xn += Q_K_(i_group) * V_IK_(i_phase, i_spec, i_group);

      // Calculate THETA_m Eq. 9
      for (int m=0; m<NUM_GROUP_; m++)
        THETA_M_(m) = Q_K_(m) * V_IK_(i_phase, i_spec, m) / sum_Qn_Xn;

      // Calculate ln(GAMMA_k^(i))
      for (int k=0; k<NUM_GROUP_; k++) {
        double sum_m_A = 0.0; // ln( sum_m_A ) term in Eq. 8
        double sum_m_B = 0.0; // last term in Eq. 8
        for (int m=0; m<NUM_GROUP_; m++) {
          sum_m_A += THETA_M_(m) * PSI_MN_(m,k);
          double sum_n = 0.0;
          for (int n=0; n<NUM_GROUP_; n++)
            sum_n += THETA_M_(n) * PSI_MN_(n,m);
          sum_m_B += THETA_M_(m) * PSI_MN_(k,m) / sum_n;
        }
        LN_GAMMA_IK_(i_phase, i_spec, k) = Q_K_(k) *
                  (1.0 - log(sum_m_A) - sum_m_B);
      }
    }
  }
}

/** \brief Perform the sub-model calculations for the current model state
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param model_data Pointer to the model data including the current state and
 *                   environmental conditions
 */
void sub_model_UNIFAC_calculate(int *sub_model_int_data,
    double *sub_model_float_data, ModelData *model_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  // Loop through each instance of each phase to calculate activity
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_instance=0; i_instance<NUM_PHASE_INSTANCE_(i_phase);
        i_instance++) {

      // Get the total number of moles of species in this phase instance
      double total_umoles = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
        total_umoles += model_data->state[PHASE_INST_ID_(i_phase, i_instance)
          + SPEC_ID_(i_phase, i_spec)] / MW_I_(i_phase, i_spec);
      }

      // If there are no species present, skip the calculations
      if (total_umoles<=0.0) {
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          GAMMA_I_(i_phase, i_instance, i_spec) = 1.0;
        }
        continue;
      }

      // Update the mole fractions X_i
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
        X_I_(i_phase, i_spec) = model_data->state[
          PHASE_INST_ID_(i_phase, i_instance) + SPEC_ID_(i_phase, i_spec)]
          / MW_I_(i_phase, i_spec) / total_umoles;
      }

      // Calcualte the sum (Q_n * X_n) in the denominator of Eq. 9 for the
      // mixture
      double sum_Qn_Xn_mixture = 0.0;
      for (int n=0; n<NUM_GROUP_; n++) {
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          sum_Qn_Xn_mixture += Q_K_(n) * X_I_(i_phase, i_spec)
            * V_IK_(i_phase, i_spec, n);
        }
      }

      // Calculate the group surface area fractions THETA_m (Eq. 9)
      for (int m=0; m<NUM_GROUP_; m++) {
        THETA_M_(m) = 0.0;
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          THETA_M_(m) += Q_K_(m) * X_I_(i_phase, i_spec)
            * V_IK_(i_phase, i_spec, m);
        }
        THETA_M_(m) /= sum_Qn_Xn_mixture;
      }

      // Calculate activity coefficients for each species in the phase instance
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {

        // Calculate PHI_i (Eq. 4)
        double PHI_i = 0.0;
        for (int j_spec=0; j_spec<NUM_SPEC_(i_phase); j_spec++)
          PHI_i += R_I_(i_phase, j_spec) * X_I_(i_phase, j_spec);
        PHI_i = R_I_(i_phase, i_spec) * X_I_(i_phase, i_spec) / PHI_i;

        // Calculate THETA_i (Eq. 4)
        double THETA_i = 0.0;
        for (int j_spec=0; j_spec<NUM_SPEC_(i_phase); j_spec++)
          THETA_i += Q_I_(i_phase, j_spec) * X_I_(i_phase, j_spec);
        THETA_i = Q_I_(i_phase, i_spec) * X_I_(i_phase, i_spec) / THETA_i;

        // Calculate the combinatorial term ln(gamma_i^C) Eq. 3
        double lnGAMMA_I_C = 0.0;
        if (X_I_(i_phase, i_spec) > 0.0) {
          for (int j_spec=0; j_spec<NUM_SPEC_(i_phase); j_spec++)
            lnGAMMA_I_C += L_I_(i_phase, j_spec) * X_I_(i_phase, j_spec);
          lnGAMMA_I_C = log( PHI_i / X_I_(i_phase, i_spec) )
                       + 5.0 * Q_I_(i_phase, i_spec) * log( THETA_i / PHI_i )
                       + L_I_(i_phase, i_spec)
                       - PHI_i / X_I_(i_phase, i_spec) * lnGAMMA_I_C;
        }

        // Calculate the residual term ln(gamma_i^R) Eq. 7
        double lnGAMMA_I_R = 0.0;
        for (int k=0; k<NUM_GROUP_; k++) {
          double sum_m_A = 0.0; // ln( sum_m_A ) term in Eq. 8
          double sum_m_B = 0.0; // last term in Eq. 8
          for (int m=0; m<NUM_GROUP_; m++) {
            sum_m_A += THETA_M_(m) * PSI_MN_(m,k);
            double sum_n = 0.0;
            for (int n=0; n<NUM_GROUP_; n++)
              sum_n += THETA_M_(n) * PSI_MN_(n,m);
            sum_m_B += THETA_M_(m) * PSI_MN_(k,m) / sum_n;
          }
          // calculate ln(GAMMA_k) Eq. 8
          double ln_GAMMA_k = Q_K_(k) * (1.0 - log(sum_m_A) - sum_m_B);
          lnGAMMA_I_R += V_IK_(i_phase, i_spec, k)
                       * ( ln_GAMMA_k - LN_GAMMA_IK_(i_phase, i_spec, k));
        }

        // Calculate gamma_i Eq. 1 and convert to units of (m^3/ug)
        GAMMA_I_(i_phase, i_instance, i_spec) =
          exp( lnGAMMA_I_C + lnGAMMA_I_R ) / total_umoles
          / MW_I_(i_phase, i_spec);

        // Set the parameter on the model state
        model_data->state[PHASE_INST_ID_(i_phase, i_instance)+
                          GAMMA_ID_(i_phase, i_spec)] =
          GAMMA_I_(i_phase, i_instance, i_spec);
      }
    }
  }
}

// TODO finish adding J contributions
/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current time step [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_UNIFAC_get_jac_contrib(int *sub_model_int_data,
    double *sub_model_float_data, ModelData *model_data, realtype *J,
    double time_step)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
}
#endif

/** \brief Print the sub model data
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 */
void sub_model_UNIFAC_print(int *sub_model_int_data,
    double *sub_model_float_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  printf("\n\nUNIFAC sub model\n\n");
  printf("\nint_data");
  for (int i=0; i<INT_DATA_SIZE_; i++) printf(" %d", int_data[i]);
  printf("\nfloat_data");
  for (int i=0; i<FLOAT_DATA_SIZE_; i++) printf(" %le", float_data[i]);
  printf("\nNumber of unique phases %d", NUM_UNIQUE_PHASE_);
  printf("\nNumber of UNIFAC groups %d", NUM_GROUP_);
  printf("\n*** General group data ***");
  for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
    printf("\nGroup %d:", i_group);
    printf("\n  Q_K: %le R_K: %le Theta_M: %le", Q_K_(i_group), R_K_(i_group),
           THETA_M_(i_group));
    printf("\n A_MN (by group):");
    for (int j_group=0; j_group<NUM_GROUP_; ++j_group)
      printf(" %le", A_MN_(i_group, j_group));
    printf("\n Psi_MN (by group):");
    for (int j_group=0; j_group<NUM_GROUP_; ++j_group)
      printf(" %le", PSI_MN_(i_group, j_group));
  }
  printf("\n\n*** Phase data ***");
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; ++i_phase) {
    printf("\nPhase %d:", i_phase);
    printf("\n  Number of instances: %d", NUM_PHASE_INSTANCE_(i_phase));
    printf("\n  Number of species: %d", NUM_SPEC_(i_phase));
    printf("\n  Species ids:");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %d", SPEC_ID_(i_phase, i_spec));
    printf("\n  ** Instance data **");
    for (int i_inst=0; i_inst<NUM_PHASE_INSTANCE_(i_phase); ++i_inst) {
      printf("\n  Instance %d:", i_inst);
      printf("\n    Jacobian ids:");
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
        printf(" %d", JAC_ID_(i_phase, i_inst, i_spec));
      printf("\n    gamma_i (by species):");
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
        printf(" %le", GAMMA_I_(i_phase, i_inst, i_spec));

    }
    printf("\n  R_I (by species):");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %le", R_I_(i_phase, i_spec));
    printf("\n  Q_I (by species):");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %le", Q_I_(i_phase, i_spec));
    printf("\n  L_I (by species):");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %le", L_I_(i_phase, i_spec));
    printf("\n  MW_I (by species):");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %le", MW_I_(i_phase, i_spec));
    printf("\n  X_I (by species):");
    for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
      printf(" %le", X_I_(i_phase, i_spec));
    printf("\n  ** Phase-specific group data **");
    for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
      printf("\n  Group %d:", i_group);
      printf("\n    v_ik (by species):");
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
        printf(" %d", V_IK_(i_phase, i_spec, i_group));
      printf("\n    ln_gamma_ik (by species):");
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
        printf(" %le", LN_GAMMA_IK_(i_phase, i_spec, i_group));
    }
  }
}
