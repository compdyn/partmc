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
      double total_group_moles = 0.0;
      for (int i_group=0; i_group<NUM_GROUP_; i_group++) {
        sum_Qn_Xn += Q_K_(i_group) * V_IK_(i_phase, i_spec, i_group);
        total_group_moles += V_IK_(i_phase, i_spec, i_group);
      }
      sum_Qn_Xn *= 1.0 / total_group_moles;

      // Calculate THETA_m Eq. 9
      for (int m=0; m<NUM_GROUP_; m++)
        THETA_M_(m) = Q_K_(m) * V_IK_(i_phase, i_spec, m) / sum_Qn_Xn /
                      total_group_moles;

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

  double *state = model_data->state;

  // Loop through each instance of each phase to calculate activity
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_instance=0; i_instance<NUM_PHASE_INSTANCE_(i_phase);
        i_instance++) {

      // Get the total number of moles of species in this phase instance
      double total_umoles = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
        total_umoles += state[PHASE_INST_ID_(i_phase, i_instance)
          + SPEC_ID_(i_phase, i_spec)] / MW_I_(i_phase, i_spec);
      }

      // If there are no species present, skip the calculations
      if (total_umoles<=0.0) {
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          GAMMA_I_(i_phase, i_instance, i_spec) = 1.0;
        }
        continue;
      }

      // Update the mole fractions x_i
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
        X_I_(i_phase, i_spec) = state[PHASE_INST_ID_(i_phase, i_instance) +
                                      SPEC_ID_(i_phase, i_spec)]
                                  / MW_I_(i_phase, i_spec) / total_umoles;
      }

      // Calculate the total number of moles of individual UNIFAC groups
      // and use this to create a conversion from total number of species
      // to total number of groups (for use in calculating group mole
      // fractions).
      double spec_umoles_to_group_umoles = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
        for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
          spec_umoles_to_group_umoles += V_IK_(i_phase, i_spec, i_group) *
              state[PHASE_INST_ID_(i_phase, i_instance) +
                    SPEC_ID_(i_phase, i_spec)] / MW_I_(i_phase, i_spec);
        }
      }
      spec_umoles_to_group_umoles *= 1.0/total_umoles;

      // Calcualte the sum (Q_n * X_n) in the denominator of Eq. 9 for the
      // mixture
      double sum_Qn_Xn_mixture = 0.0;
      for (int n=0; n<NUM_GROUP_; n++) {
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          sum_Qn_Xn_mixture += Q_K_(n) * X_I_(i_phase, i_spec)
            * V_IK_(i_phase, i_spec, n) / spec_umoles_to_group_umoles;
        }
      }

      // Calculate the group surface area fractions THETA_m (Eq. 9)
      for (int m=0; m<NUM_GROUP_; m++) {
        THETA_M_(m) = 0.0;
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {
          THETA_M_(m) += Q_K_(m) * X_I_(i_phase, i_spec)
            * V_IK_(i_phase, i_spec, m) / spec_umoles_to_group_umoles;
        }
        THETA_M_(m) /= sum_Qn_Xn_mixture;
      }

      // Calculate the denominators of Eq. 4
      double PHI_T = 0.0;
      for (int j_spec=0; j_spec<NUM_SPEC_(i_phase); j_spec++)
        PHI_T += R_I_(i_phase, j_spec) * X_I_(i_phase, j_spec);
      double THETA_T = 0.0;
      for (int j_spec=0; j_spec<NUM_SPEC_(i_phase); j_spec++)
        THETA_T += Q_I_(i_phase, j_spec) * X_I_(i_phase, j_spec);

      // Calculate activity coefficients for each species in the phase instance
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); i_spec++) {

        // Calculate PHI_i (Eq. 4)
        double PHI_i;
        PHI_i = R_I_(i_phase, i_spec) * X_I_(i_phase, i_spec) / PHI_T;

        // Calculate THETA_i (Eq. 4)
        double THETA_i;
        THETA_i = Q_I_(i_phase, i_spec) * X_I_(i_phase, i_spec) / THETA_T;

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
        state[PHASE_INST_ID_(i_phase, i_instance) +
              GAMMA_ID_(i_phase, i_spec)] =
          GAMMA_I_(i_phase, i_instance, i_spec);
      }
    }
  }
}

// TODO finish adding J contributions
/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * This derivation starts from equations (1)--(9) in \cite Marcolli2005.
 * The mole fraction of species \f$i\f$ is calculated as:
 * \f[
 *   \chi_i = \frac{m_i}{m_T},
 * \f]
 * where \f$m_T = \displaystyle\sum_j m_j\f$,
 * \f$m_i = \frac{c_i}{\mathrm{MW}_i}\f$
 * is the number in moles of species \f$i\f$,
 * \f$c_i\f$ is its mass concentration
 * (\f$\mathrm{ug} \: \mathrm{m}^{-3}\f$; the state variable), and
 * \f$\mathrm{MW}_i\f$ is its molecular weight
 * (\f$\mathrm{ug} \: \mathrm{mol}^{-1}\f$). Thus,
 * \f[
 * \begin{align*}
 *   \frac{\partial \chi_i}{\partial c_i} & =
 *      \frac{(m_T-m_i)}{\mathrm{MW}_i m_T^2}, \textrm{ and} \\
 *   \frac{\partial \chi_j}{\partial c_i} & =
 *      -\frac{m_i}{\mathrm{MW}_i m_T^2} \quad \text{for } i\neq j. \\
 * \end{align*}
 * \f]
 *
 * The partial derivative of \f$\Phi_i\f$ (Eq. 4) with respect to \f$c_i\f$
 * is derived as follows:
 * \f[
 * \begin{align*}
 *   \frac{\partial r_i x_i}{\partial c_i} & =
 *       r_i \frac{(m_T-m_i)}{\mathrm{MW}_i m_T^2}, \\
 *   \frac{\partial r_j x_j}{\partial c_i} & =
 *       - r_i \frac{m_i}{\mathrm{MW}_i m_T^2}
 *       \quad \text{for } i\neq j, \\
 *   \frac{\partial \displaystyle\sum_j r_j x_j}{\partial c_i} & =
 *       \frac{r_i}{\mathrm{MW}_i m_T} -
 *       \displaystyle\sum_j\frac{r_j m_j}{\mathrm{MW}_j m_T^2}, \\
 *   \rho & \equiv \displaystyle\sum_j\frac{r_j m_j}{\mathrm{MW}_j m_T^2},
 *   \quad \sigma \equiv \displaystyle\sum_k r_k x_k, \\
 *   \frac{\partial \Phi_j}{\partial c_i} & = \frac{
 *       \left[r_i \frac{(m^{\prime}_T-m_i)}{\mathrm{MW}_i m_T^2}
 *        \sigma -
 *        r_i x_i \left(\frac{r_i}{\mathrm{MW}_i m_T^2} - \rho
 *        \right) \right]} {\sigma^2}, \qquad
 *        m^{\prime}_T =
 *        \begin{cases}
 *          m_T  & \quad \text{if } i=j \\
 *          0    & \quad \text{if } i\neq j
 *        \end{cases}
 * \end{align*}
 * \f]
 *
 * Similarly, the partial derivative of \f$\Theta_i\f$ (Eq. 4) with respect
 * to \f$c_i\f$ is:
 * \f[
 * \begin{align*}
 *   \tau & \equiv \displaystyle\sum_j\frac{q_j m_j}{\mathrm{MW}_j m_T^2},
 *   \quad \omega \equiv \displaystyle\sum_k q_k x_k, \\
 *   \frac{\partial \Theta_j}{\partial c_i} & = \frac{
 *       \left[q_i \frac{(m^{\prime}_T-m_i)}{\mathrm{MW}_i m_T^2}
 *        \omega -
 *        q_i x_i \left(\frac{q_i}{\mathrm{MW}_i m_T^2} - \tau
 *        \right) \right]} {\omega^2}, \qquad
 *        m^{\prime}_T =
 *        \begin{cases}
 *          m_T  & \quad \text{if } i=j \\
 *          0    & \quad \text{if } i\neq j
 *        \end{cases}
 * \end{align*}
 * \f]
 *
 * From Eqs 5 and 6,
 * \f[
 *   \frac{\partial r_i}{\partial c_i} =
 *   \frac{\partial q_i}{\partial c_i} =
 *   \frac{\partial l_i}{\partial c_i} = 0.
 * \f]
 * For the last term in Eq. 3,
 * \f[
 * \begin{align*}
 *   \frac{\partial \displaystyle\sum_j x_j l_j}{\partial c_i} = &
 *       \frac{l_i}{\mathrm{MW}_i m_T} -
 *       \displaystyle\sum_j\frac{l_j m_j}{\mathrm{MW}_j m_T^2}, \\
 *   \kappa \equiv & \displaystyle\sum_j\frac{l_j m_j}{\mathrm{MW}_j m_T^2},
 *   \quad \mu \equiv \displaystyle\sum_j x_j l_j, \\
 *   \frac{\partial \mu}{\partial c_i} = & \frac{l_i}{\mathrm{MW}_i m_T}
 *       - \kappa.
 * \end{align*}
 * \f]
 * Thus, the partial derivative of the full combinatorial term (Eq. 3) is:
 * \f[
 * \begin{align*}
 *   \frac{\partial \ln{\gamma^C_j}}{\partial c_i} = &
 *       \frac{x_j}{\Phi_j}\frac{\left( \frac{\partial \Phi_j}{\partial c_i}
 *           x_j - \Phi_j\frac{\partial x_j}{\partial c_i} \right)}{x_j^2} \\
 *     & \quad + \frac{z}{2}q_j\frac{\Phi_j}{\Theta_j}\frac{\left(
 *           \frac{\partial \Theta_j}{\partial c_i} \Phi_j -
 *           \Theta_j\frac{\partial \Phi_j}{\partial c_i}\right)}{\Phi_j^2} \\
 *     & \quad - \mu\frac{\partial \mu}{\partial c_i}
 *         \frac{\left(\frac{\partial \Phi_j}{\partial c_i}x_j -
 *              \Phi_j\frac{\partial x_j}{\partial c_i} \right)}{x_j^2}.
 * \end{align*}
 * \f]
 * After some rearranging, this becomes:
 * \f[
 * \begin{align*}
 *   \frac{\partial \ln{\gamma^C_j}}{\partial c_i} = &
 *       \frac{1}{\Phi_j}\frac{\partial \Phi_j}{\partial c_i}
 *       - \frac{1}{x_j}\frac{\partial x_j}{\partial c_i} \\
 *     & \quad + \frac{z}{2}q_j\left(
 *       \frac{1}{\Theta_j}\frac{\partial \Theta_j}{\partial c_i}
 *       - \frac{1}{\Phi_j}\frac{\partial \Phi_j}{\partial c_i}\right) \\
 *     & \quad - \frac{\mu}{x_j}\frac{\partial \mu}{\partial c_i} \left(
 *       \frac{\partial \Phi_j}{\partial c_i} - \frac{\Phi_j}{x_j}
 *       \frac{\partial x_j}{\partial c_i} \right)
 * \end{align*}
 * \f]
 * As \f$\rho\f$, \f$\sigma\f$, \f$\tau\f$ \f$\omega\f$, \f$\kappa\f$,
 * and \f$\mu\f$ are independent of species \f$i\f$,
 * these can be calculated outside the loop over the independent species.
 * Moving to the residual term (Eq. 7), the mole fraction (\f$X_k\f$)
 * of group \f$k\f$ in the mixture is related to the mole fraction
 * (\f$x_i\f$) of species \f$i\f$ according to:
 * \f[
 *   X_k = \displaystyle\sum_i v_k^{(i)} x_i c_{xX},
 * \f]
 * where \f$c_{xX}\f$ is a conversion factor accounting for the difference
 * in total species and total group number concentrations:
 * \f[
 *   c_{xX} = \frac{m_T}{\displaystyle\sum_i v_k^{(i)} m_i}
 * \f]
 * Partial derivatives of the group interaction terms in Eq 9,
 * \f$\Theta_m\f$ and \f$\Psi_{mn}\f$, with respect to \f$c_i\f$ are
 * derived as follows:
 * \f[
 * \begin{align*}
 *   \frac{\partial X_k}{\partial c_i} = &
 *      \frac{v_k^{(i)}c_{xX}}{\mathrm{MW}_i m_T} -
 *      \displaystyle\sum_j \frac{v_k^{(j)} m_j c_{xX}}
 *      {\mathrm{MW}_j m_T^2}, \\
 *   \pi_k \equiv & \displaystyle\sum_j \frac{v_k^{(j)} m_j c_{xX}}
 *       {\mathrm{MW}_j m_T^2}, \\
 *   \frac{\partial \displaystyle\sum_n Q_n X_n}{\partial c_i} = &
 *       \displaystyle\sum_n Q_n \frac{v_n^{(i)} c_{xX}}
 *       {\mathrm{MW}_i m_T}
 *       - \displaystyle\sum_n\displaystyle\sum_j Q_n
 *       \frac{v_n^{(j)} m_j c_{xX}}{\mathrm{MW}_j m_T^2}, \\
 *   \Pi \equiv & \displaystyle\sum_n Q_n \pi_n =
 *       \displaystyle\sum_n\displaystyle\sum_j Q_n
 *       \frac{v_n^{(j)} m_j c_{xX}}{\mathrm{MW}_j m_T^2}, \\
 *   \frac{\partial \Theta_m}{\partial c_i} = & \left[
 *       \left(Q_m\frac{v_m^{(i)} c_{xX}}{\mathrm{MW}_i m_T^2} -
 *          Q_m\pi_m \right)
 *       \displaystyle\sum_n Q_n X_n
 *       - \left(\frac{v_m^{(i)} c_{xX}}{\mathrm{MW}_i m_T^2} - \Pi \right)
 *       Q_m X_m \right] / \left(\displaystyle\sum_n Q_n X_n \right)^2,
 *       \textrm{ and} \\
 *   \frac{\partial \Psi_{mn}}{\partial c_i} = & 0
 * \end{align*}
 * \f]
 * The partial derivative of the group residual activity coefficient (Eq. 8)
 * with respect to \f$c_i\f$ is thus:
 * \f[
 * \begin{align*}
 *   \frac{\partial \Gamma_k}{\partial c_i} = & Q_k \left[
 *     - \frac{1}{\displaystyle\sum_m\Theta_m\Psi_{mk}}
 *       \displaystyle\sum_m\frac{\partial\Theta_m}{\partial c_i}\Psi_{mk} \\
 *     \quad - \displaystyle\sum_m\left(
 *       \frac{\partial\Theta_m}{\partial c_i}\Psi_{km}
 *         \displaystyle\sum_n\Theta_n\Psi_{nm}
 *     - \Theta_m\Psi_{km}\displaystyle\sum_n
 *         \frac{\partial\Theta_n}{\partial c_i}
 *         \Psi_{nm}\right)/\left(
 *           \displaystyle\sum_n\Theta_n\Psi_{nm}\right)^2
 *     \right]
 * \end{align*}
 * \f]
 * After some rearranging, this becomes:
 * \f[
 *   \frac{\partial \Gamma_k}{\partial c_i} = - Q_k \displaystyle\sum_m \left(
 *     \frac{1}{\Theta_m}\frac{\partial\Theta_m}{\partial c_i}
 *     + \frac{\Psi_{km}}{\displaystyle\sum_n\Theta_n\Psi_{nm}}
 *       \frac{\partial\Theta_m}{\partial c_i}
 *     + \frac{\Theta_m\Psi_{km}}
 *            {\left(\displaystyle\sum_n\Theta_n\Psi_{nm}\right)^2}
 *       \displaystyle\sum_n\frac{\partial\Theta_n}{\partial c_i}\Psi_{nm}
 *   \right)
 * \f]
 * The left side of the three bracketed terms in the above equation are
 * independent of species \f$i\f$ and can be calculated outside of the loop
 * over the independent species. The partial derivative of the full residual
 * term with respect to \f$c_i\f$ is thus:
 * \f[
 * \begin{align*}
 *  \frac{\partial \ln{\Gamma_k^{(j)}}}{\partial c_i} = & 0 \\
 *  \frac{\partial \ln{\gamma_j^R}}{\partial c_i} = &
 *    \displaystyle\sum_k \frac{v_k^{(j)}}{\Gamma_k}
 *       \frac{\partial \Gamma_k}{\partial c_i}
 * \end{align*}
 * \f]
 * The overall equation for the partial derivative of \f$\gamma_j\f$ with
 * respect to species \f$i\f$ is:
 * \f[
 *   \frac{\partial\gamma_j}{\partial c_i} = \gamma_j \left(
 *     \frac{\partial\ln{\gamma_j^C}}{\partial c_i}
 *     + \frac{\partial\ln{\gamma_j^R}}{\partial c_i} \right)
 * \f]
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

  // Calculate the model parameters
  sub_model_UNIFAC_calculate(sub_model_int_data, sub_model_float_data,
                             model_data);

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
