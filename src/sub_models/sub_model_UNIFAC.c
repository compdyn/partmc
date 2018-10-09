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
#define V_IK_(p,i,k) (int_data[PHASE_INT_LOC_(p)+2+2*NUM_PHASE_INSTANCE_(p)+(k+1)*NUM_SPEC_(p)+i])

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
 * \param sub_model_data Pointer to the sub-model data
 * \param jac_row Array of flags indicating whether an element in the rown is
 *                used
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_get_used_jac_elem(void *sub_model_data, bool *jac_row)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Get the id of a parameter in the condensed data block
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param identifiers For the UNIFAC model, the identifers are just the id
 *                    on the state array of the aerosol-phase species for
 *                    which the acitivty is needed.
 * \param parameter_id Parameter id for the requested activity coefficient if
 *                     found
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_get_parameter_id(void *sub_model_data,
          void *identifiers, int *parameter_id)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_instance=0; i_instance<NUM_PHASE_INSTANCE_(i_phase);
        i_instance++) {
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
        if (*((int*)identifiers) == SPEC_ID_(i_phase, i_spec)
            + PHASE_INST_ID_(i_phase, i_instance)) {
          *parameter_id = (int) (((int*)
                (&(GAMMA_I_(i_phase, i_instance, i_spec)))) - int_data);
          return (void*) &(float_data[FLOAT_DATA_SIZE_]);
        }
      }
    }
  }
  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Update sub-model data for new environmental conditions
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param env_data Pointer to the environmental state array
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_update_env_state(void *sub_model_data,
          double *env_data)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
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
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

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

      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

// TODO finish adding J contributions
/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * Derivatives are assumed to be of the form \f$\frac{dy}{dt} = A*S\f$, where
 * \f$A\f$ is the value passed to this function as \b base_val and \f$S\f$ is
 * the sub-model parameter used in the calculation. The row of the Jacobian
 * should correspond to \f$\frac{dy'}{dx}\f$, where for each element \f$x\f$,
 * on the row, this function will add \f$A*\frac{dS}{dx}\f$.
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param base_val The derivative
 * \param jac_row Pointer to the Jacobian row to modify
 */
void * sub_model_UNIFAC_add_jac_contrib(void *sub_model_data,
         double base_val, double *jac_row)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Skip through the sub model, only advancing the data pointer
 *
 * \param sub_model_data Pointer to the sub model data
 * \return The sub_model_data pointer advanced by the size of the sub-model
 */
void * sub_model_UNIFAC_skip(void *sub_model_data)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the sub model data
 *
 * \param sub_model_data Pointer to the sub model data
 * \return The sub_model_data pointer advanced by the size of the sub-model
 */
void * sub_model_UNIFAC_print(void *sub_model_data)
{
  int *int_data = (int*) sub_model_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nUNIFAC sub model\n\n");
  printf("\nint_data");
  for (int i=0; i<INT_DATA_SIZE_; i++) printf(" %d", int_data[i]);
  printf("\nfloat_data");
  for (int i=0; i<FLOAT_DATA_SIZE_; i++) printf(" %le", float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef NUM_UNIQUE_PHASE_
#undef NUM_GROUP_
#undef TOTAL_INT_PROP_
#undef TOTAL_FLOAT_PROP_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef PHASE_INT_LOC_
#undef PHASE_FLOAT_LOC_
#undef NUM_PHASE_INSTANCE_
#undef NUM_SPEC_
#undef PHASE_INST_FLOAT_LOC_
#undef PHASE_INST_ID_
#undef SPEC_ID_
#undef V_IK_

#undef Q_K_
#undef R_K_
#undef THETA_M_
#undef A_MN_
#undef PSI_MN_
#undef R_I_
#undef Q_I_
#undef L_I_
#undef MW_I_
#undef X_I_
#undef LN_GAMMA_IK_
#undef GAMMA_I_

#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
