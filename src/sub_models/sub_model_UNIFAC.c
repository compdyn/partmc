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
#define PHASE_ENV_LOC_(p) (int_data[NUM_INT_PROP_+2*NUM_UNIQUE_PHASE_+p]-1)
#define NUM_PHASE_INSTANCE_(p) (int_data[PHASE_INT_LOC_(p)])
#define NUM_SPEC_(p) (int_data[PHASE_INT_LOC_(p)+1])
#define PHASE_INST_ID_(p,c) (int_data[PHASE_INT_LOC_(p)+2+c]-1)
#define SPEC_ID_(p,i) (int_data[PHASE_INT_LOC_(p)+2+NUM_PHASE_INSTANCE_(p)+i])
#define GAMMA_ID_(p,i) (int_data[PHASE_INT_LOC_(p)+2+NUM_PHASE_INSTANCE_(p)+NUM_SPEC_(p)+i])
#define JAC_ID_(p,c,j,i) int_data[PHASE_INT_LOC_(p)+2+NUM_PHASE_INSTANCE_(p)+(2+c*NUM_SPEC_(p)+j)*NUM_SPEC_(p)+i]
#define V_IK_(p,i,k) (int_data[PHASE_INT_LOC_(p)+2+NUM_PHASE_INSTANCE_(p)+(k+2+NUM_PHASE_INSTANCE_(p)*NUM_SPEC_(p))*NUM_SPEC_(p)+i])

#define Q_K_(k) (float_data[k])
#define R_K_(k) (float_data[NUM_GROUP_+k])
#define X_K_(m) (float_data[2*NUM_GROUP_+m])
#define DTHETA_M_DC_I_(m) (float_data[3*NUM_GROUP_+m])
#define XI_M_(m) (float_data[4*NUM_GROUP_+m])
#define LN_GAMMA_K_(m) (float_data[5*NUM_GROUP_+m])
#define A_MN_(m,n) (float_data[(m+6)*NUM_GROUP_+n])
#define R_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+i])
#define Q_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+NUM_SPEC_(p)+i])
#define L_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+2*NUM_SPEC_(p)+i])
#define MW_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+3*NUM_SPEC_(p)+i])
#define X_I_(p,i) (float_data[PHASE_FLOAT_LOC_(p)+4*NUM_SPEC_(p)+i])

#define THETA_M_(m) (sub_model_env_data[m])
#define PSI_MN_(m,n) (sub_model_env_data[(m+1)*NUM_GROUP_+n])
#define LN_GAMMA_IK_(p,i,k) (sub_model_env_data[PHASE_ENV_LOC_(p)+i*NUM_GROUP_+k])

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
      for (int i_gamma=0; i_gamma < NUM_SPEC_(i_phase); ++i_gamma)
        for (int i_spec=0; i_spec < NUM_SPEC_(i_phase); ++i_spec)
          jac_struct[PHASE_INST_ID_(i_phase, i_inst)+GAMMA_ID_(i_phase, i_gamma)]
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
      for (int i_gamma=0; i_gamma < NUM_SPEC_(i_phase); ++i_gamma)
        for (int i_spec=0; i_spec < NUM_SPEC_(i_phase); ++i_spec)
          JAC_ID_(i_phase, i_inst, i_gamma, i_spec) =
            jac_ids[PHASE_INST_ID_(i_phase, i_inst)+GAMMA_ID_(i_phase, i_gamma)]
                   [PHASE_INST_ID_(i_phase, i_inst)+SPEC_ID_(i_phase, i_spec)];
}

/** \brief Update sub-model data for new environmental conditions
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 */
void sub_model_UNIFAC_update_env_state(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *env_data = model_data->grid_cell_env;

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
 * These calculations are based on equations (1)--(9) in \cite Marcolli2005.
 * Variable names used in the derivation of the Jacobian equations are
 * adopted here for consistency. Specifically, these are:
 * \f[
 *   \begin{align*}
 *    \sigma & \equiv \displaystyle\sum_k r_k x_k, \\
 *    \tau & \equiv \displaystyle\sum_k q_k x_k, \\
 *    \mu & \equiv \displaystyle\sum_j x_j l_j, \\
 *    c_{xX} & = \frac{1}{\displaystyle\sum_j
 *                   \displaystyle\sum_i v_j^{(i)} x_i}, \textrm{ and} \\
 *    \Pi & \equiv \displaystyle\sum_n Q_n X_n. \\
 *   \end{align*}
 * \f]
 * Using these variables and subscript \f$j\f$ for the species whose activity
 * is being calculated for easy comparison with the Jacobian calculations,
 * equations (1)--(9) in \cite Marcolli2005 become:
 * \f[
 *   \begin{align*}
 *    \ln{\gamma_j} & = \ln{\gamma_j^C} + \ln{\gamma_j^R} & (1) \\
 *    \alpha_w & = \gamma_w x_w & (2) \\
 *    \ln{\gamma_j^C} & = \ln{\frac{\Phi_j}{x_j}}
 *                       + \frac{z}{2}q_j\ln{\frac{\Theta_j}{\Phi_j}}
 *                       + l_j
 *                       - \frac{\Phi_j}{x_j}\mu & (3) \\
 *    \Phi_j & = \frac{r_j x_j}{\sigma}; \qquad
 *    \Theta_j = \frac{q_j x_j}{\tau} & (4) \\
 *    l_j & = \frac{z}{2}\left(r_j - q_j\right) - r_j + 1 & (5) \\
 *    r_j & = \displaystyle\sum_k v_k^{(j)}R_k; \qquad
 *    q_j   = \displaystyle\sum_k v_k^{(j)}Q_k & (6) \\
 *    \ln{\gamma_j^R} & = \displaystyle\sum_k v_k^{(j)} \left[
 *                        \ln{\Gamma_k} - \ln{\Gamma_k^{(j)}}\right] & (7) \\
 *    \ln{\Gamma_k} & = Q_k \left[ 1 - \ln{\Xi_k} -
 *                      \displaystyle\sum_m\frac{\Theta_m\Phi_{km}}{\Xi_m}
 *                      \right] & (8) \\
 *    \Theta_m & = \frac{Q_m X_m}{\Pi}; \qquad
 *    \Psi_{mn} = \mathrm{exp}\left[\frac{-a_{mn}}{T}\right] & (9) \\
 *   \end{align*}
 * \f]
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data including the current state and
 *                   environmental conditions
 */
void sub_model_UNIFAC_calculate(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;

  double *state = model_data->grid_cell_state;

  // Loop through each instance of each phase to calculate activity
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; ++i_phase) {
    for (int i_instance=0; i_instance<NUM_PHASE_INSTANCE_(i_phase);
        ++i_instance) {

      // Get the total number of moles of species in this phase instance
      double m_T = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
        m_T += state[PHASE_INST_ID_(i_phase, i_instance)
          + SPEC_ID_(i_phase, i_spec)] / MW_I_(i_phase, i_spec);
      }

      // If there are no species present, skip the calculations
      if (m_T<=0.0) {
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
          state[PHASE_INST_ID_(i_phase, i_instance) +
                GAMMA_ID_(i_phase, i_spec)] = 0.0;
        }
        continue;
      }

      // Pre-calculate mole fractions and variables that do not depend
      // on species j

      double sigma = 0.0;
      double tau   = 0.0;
      double mu    = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {

        // Mole fractions x_i
        double x_i = state[PHASE_INST_ID_(i_phase, i_instance) +
                           SPEC_ID_(i_phase, i_spec)]
                     / MW_I_(i_phase, i_spec) / m_T;
        X_I_(i_phase, i_spec) = x_i;

        // Combinatorial variables
        sigma += R_I_(i_phase, i_spec) * x_i;
        tau   += Q_I_(i_phase, i_spec) * x_i;
        mu    += L_I_(i_phase, i_spec) * x_i;
      }

      // Residual variables
      double c_xX = 0.0;
      double Pi   = 0.0;
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        X_K_(i_group)  = 0.0;
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
          double v_ik_x_i = V_IK_(i_phase, i_spec, i_group)
                            * X_I_(i_phase, i_spec);
          c_xX           += v_ik_x_i;
          X_K_(i_group)  += v_ik_x_i;
        }
      }
      c_xX = 1.0 / c_xX;
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        X_K_(i_group)    *= c_xX;
        Pi               += Q_K_(i_group) * X_K_(i_group);
      }
      // Group surface area fractions Theta_m (Eq. 9)
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group)
        THETA_M_(i_group) = Q_K_(i_group) * X_K_(i_group) / Pi;

      // Xi_m
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        XI_M_(i_group) = 0.0;
        for (int j_group=0; j_group<NUM_GROUP_; ++j_group) {
          XI_M_(i_group) += THETA_M_(j_group) * PSI_MN_(j_group, i_group);
        }
      }

      // Group residual activity coefficients ln(Gamma_k) (Eq. 8)
      for (int k=0; k<NUM_GROUP_; ++k) {
        LN_GAMMA_K_(k) = 1.0 - log( XI_M_(k) );
        for (int m=0; m<NUM_GROUP_; ++m)
          LN_GAMMA_K_(k) -= THETA_M_(m) * PSI_MN_(k,m) / XI_M_(m);
        LN_GAMMA_K_(k) *= Q_K_(k);
      }

      // Calculate activity coefficients for each species in the phase instance
      for (int j=0; j<NUM_SPEC_(i_phase); ++j) {

        // Mole fraction for species j
        double x_j = X_I_(i_phase, j);

        // Phi_j (Eq. 4)
        double Phi_j;
        Phi_j = R_I_(i_phase, j) * x_j / sigma;

        // Theta_j (Eq. 4)
        double Theta_j;
        Theta_j = Q_I_(i_phase, j) * x_j / tau;

        // Calculate the combinatorial term ln(gamma_i^C) Eq. 3
        double lngammaC_j = 0.0;
        if (X_I_(i_phase, j) > 0.0) {
          lngammaC_j = log( Phi_j / x_j )
                      + 5.0 * Q_I_(i_phase, j) * log( Theta_j / Phi_j )
                      + L_I_(i_phase, j)
                      - Phi_j / x_j * mu;
        }

        // Full residual term
        double lngammaR_j = 0.0;
        for (int k=0; k<NUM_GROUP_; ++k) {
          lngammaR_j += V_IK_(i_phase, j, k) *
                        ( LN_GAMMA_K_(k) - LN_GAMMA_IK_(i_phase, j, k) );
        }

        // Calculate gamma_i Eq. 1
        state[PHASE_INST_ID_(i_phase, i_instance) +
              GAMMA_ID_(i_phase, j)] =
          exp( lngammaC_j + lngammaR_j ) * X_I_(i_phase, j);
      }
    }
  }
}

/** \brief Add contributions to the Jacobian from derivates calculated using the output of this sub model
 *
 * This derivation starts from equations (1)--(9) in \cite Marcolli2005.
 * The mole fraction of species \f$i\f$ is calculated as:
 * \f[
 *    x_i = \frac{m_i}{m_T},
 * \f]
 *  where \f$m_T = \displaystyle\sum_j m_j\f$,
 *  \f$m_i = \frac{c_i}{\mathrm{MW}_i}\f$
 *  is the number in moles of species \f$i\f$,
 *  \f$c_i\f$ is its mass concentration
 *  (\f$\mathrm{ug} \: \mathrm{m}^{-3}\f$; the state variable), and
 *  \f$\mathrm{MW}_i\f$ is its molecular weight
 *  (\f$\mathrm{ug} \: \mathrm{mol}^{-1}\f$). Thus,
 * \f[
 *   \begin{align*}
 *    \frac{\partial x_i}{\partial c_i} & =
 *       \frac{(m_T-m_i)}{\mathrm{MW}_i m_T^2}, \textrm{ and} \\
 *    \frac{\partial x_j}{\partial c_i} & =
 *       -\frac{m_j}{\mathrm{MW}_i m_T^2} \quad \text{for } i\neq j. \\
 *   \end{align*}
 * \f]
 *  The partial derivative of \f$\Phi_j\f$ (Eq. 4) with respect to \f$c_i\f$
 *  is derived as follows:
 * \f[
 *   \begin{align*}
 *     \frac{\partial r_i x_i}{\partial c_i} & =
 *        r_i \frac{(m_T-m_i)}{\mathrm{MW}_i m_T^2}, \\
 *     \frac{\partial r_j x_j}{\partial c_i} & =
 *        - r_j \frac{m_j}{\mathrm{MW}_i m_T^2}
 *        \quad \text{for } i\neq j, \\
 *     \frac{\partial \displaystyle\sum_j r_j x_j}{\partial c_i} & =
 *        \frac{1}{\mathrm{MW}_i m_T}\left[r_i -
 *        \displaystyle\sum_j r_j x_j\right], \\
 *     \sigma & \equiv \displaystyle\sum_k r_k x_k, \\
 *     \frac{\partial \Phi_j}{\partial c_i} & =
 *     \frac{r_j}{\mathrm{MW}_i m_T \sigma^2}
 *     \left(\sigma^{\prime} - x_j r_i \right)
 *     , \qquad
 *         \sigma^{\prime} =
 *         \begin{cases}
 *           \sigma  & \quad \text{if } i=j \\
 *           0    & \quad \text{if } i\neq j
 *         \end{cases} \\
 *   \end{align*}
 * \f]
 *  Similarly, the partial derivative of \f$\Theta_j\f$ (Eq. 4) with respect
 *  to \f$c_i\f$ is:
 * \f[
 * \begin{align*}
 *    \tau & \equiv \displaystyle\sum_k q_k x_k, \\
 *    \frac{\partial \Theta_j}{\partial c_i} & =
 *     \frac{q_j}{\mathrm{MW}_i m_T \tau^2}
 *     \left(\tau^{\prime} - x_j q_i \right)
 *     , \qquad
 *         \tau^{\prime} =
 *         \begin{cases}
 *           \tau  & \quad \text{if } i=j \\
 *           0    & \quad \text{if } i\neq j
 *         \end{cases} \\
 * \end{align*}
 * \f]
 *  From Eqs 5 and 6,
 * \f[
 *    \frac{\partial r_j}{\partial c_i} =
 *    \frac{\partial q_j}{\partial c_i} =
 *    \frac{\partial l_j}{\partial c_i} = 0.
 * \f]
 *  For the last term in Eq. 3,
 * \f[
 * \begin{align*}
 *    \mu & \equiv \displaystyle\sum_j x_j l_j, \\
 *    \frac{\partial \mu}{\partial c_i} & = \frac{1}{\mathrm{MW}_i m_T}
 *        \left( l_i - \mu \right). \\
 * \end{align*}
 * \f]
 *  The partial derivative of the full combinatorial term (Eq. 3) is:
 * \f[
 * \begin{align*}
 *    \frac{\partial \ln{\gamma^C_j}}{\partial c_i} & =
 *        \frac{x_j}{\Phi_j}\frac{\left( \frac{\partial \Phi_j}{\partial c_i}
 *            x_j - \Phi_j\frac{\partial x_j}{\partial c_i} \right)}{x_j^2} \\
 *      & \quad + \frac{z}{2}q_j\frac{\Phi_j}{\Theta_j}\frac{\left(
 *            \frac{\partial \Theta_j}{\partial c_i} \Phi_j -
 *            \Theta_j\frac{\partial \Phi_j}{\partial c_i}\right)}{\Phi_j^2} \\
 *      & \quad - \frac{\left[ \left( \frac{\partial \Phi_j}{\partial c_i} \mu +
 *             \Phi_j \frac{\partial \mu}{\partial c_i} \right) x_j
 *             - \Phi_j \mu \frac{\partial x_j}{\partial c_i} \right]}{x_j^2}.
 * \end{align*}
 * \f]
 *  After some rearranging, this becomes:
 * \f[
 * \begin{align*}
 *    \frac{\partial \ln{\gamma^C_j}}{\partial c_i} & =
 *        \frac{1}{\Phi_j}\frac{\partial \Phi_j}{\partial c_i}
 *        - \frac{1}{x_j}\frac{\partial x_j}{\partial c_i} \\
 *      & \quad + \frac{z}{2}q_j\left(
 *        \frac{1}{\Theta_j}\frac{\partial \Theta_j}{\partial c_i}
 *        - \frac{1}{\Phi_j}\frac{\partial \Phi_j}{\partial c_i}\right) \\
 *      & \quad - \frac{\partial \Phi_j}{\partial c_i}\frac{\mu}{x_j}
 *      - \frac{\Phi_j}{x_j}\frac{\partial \mu}{\partial c_i}
 *      + \frac{\Phi_j \mu}{x_j^2}\frac{\partial x_j}{\partial c_i}
 * \end{align*}
 * \f]
 *  As \f$\sigma\f$, \f$\tau\f$,
 *  and \f$\mu\f$ are independent of species \f$i\f$,
 *  these can be calculated outside the loop over the independent species.
 *  Moving to the residual term (Eq. 7), the mole fraction (\f$X_p\f$)
 *  of group \f$p\f$ in the mixture is related to the mole fraction
 *  (\f$x_j\f$) of species \f$j\f$ according to:
 * \f[
 *    X_p = c_{xX} \omega_p, \textrm{ where }
 *    \omega_p \equiv \displaystyle\sum_j v_p^{(j)} x_j, \\
 * \f]
 *  and \f$c_{xX}\f$ is a conversion factor accounting for the difference
 *  in total species and total group number concentrations:
 * \f[
 *    c_{xX} = \frac{1}{\displaystyle\sum_p \omega_p}.
 * \f]
 *  Partial derivatives of the group interaction terms in Eq 9,
 *  \f$\Theta_m\f$ and \f$\Psi_{mn}\f$, with respect to \f$c_i\f$ are
 *  derived as follows:
 * \f[
 *   \begin{align*}
 *     \frac{\partial c_{xX}}{\partial c_i} & = - c_{xX}^2
 *       \displaystyle\sum_p \displaystyle\sum_j v_p^{(j)}
 *       \frac{\partial x_j}{\partial c_i}
 *       = - \frac{c_{xX}^2}{\mathrm{MW}_im_T} \left(
 *             \displaystyle\sum_p v_p^{(i)} -
 *             \displaystyle\sum_p \omega_p \right), \\
 *       & = \frac{c_{xX}}{\mathrm{MW}_im_T} \left(
 *               1 - c_{xX} \displaystyle\sum_p v_p^{(i)} \right), \\
 *     \frac{\partial X_n}{\partial c_i} & =
 *       c_{xX} \displaystyle\sum_j v_n^{(j)}
 *            \frac{\partial x_j}{\partial c_i}
 *       + \frac{\partial c_{xX}}{\partial c_i} \omega_n, \\
 *       & = \frac{c_{xX}}{\mathrm{MW}_im_T}
 *                  \left( v_n^{(i)} - \omega_n \right)
 *             + \frac{c_{xX}}{\mathrm{MW}_im_T}
 *                  \left( \omega_n - c_{xX}\omega_n
 *                      \displaystyle\sum_p v_p^{(i)} \right), \\
 *        & = \frac{c_{xX}}{\mathrm{MW}_im_T} \left(
 *                  v_n^{(i)} - X_n
 *                      \displaystyle\sum_p v_p^{(i)} \right), \\
 *     \Pi & \equiv \displaystyle\sum_n Q_n X_n, \\
 *     \frac{\partial \Pi}{\partial c_i} & =
 *        \displaystyle\sum_n Q_n
 *        \frac{c_{xX}}{\mathrm{MW}_i m_T}\left(
 *            v_n^{(i)} - X_n \displaystyle\sum_p v_p^{(i)} \right), \\
 *     & = \frac{c_{xX}}{\mathrm{MW}_i m_T} \left(
 *            \displaystyle\sum_n Q_n v_n^{(i)} -
 *            \Pi \displaystyle\sum_p v_p^{(i)} \right), \\
 *     \frac{\partial \Theta_m}{\partial c_i} & = \frac{
 *     \left(Q_m\frac{\partial X_m}{\partial c_i} \Pi -
 *         Q_m X_m \frac{\partial \Pi}{\partial c_i} \right)}{\Pi^2} , \\
 *     & = \frac{c_{xX} Q_m}{\mathrm{MW}_im_T\Pi^2} \left[
 *              \left( \Pi v_m^{(i)} - \Pi X_m
 *                  \displaystyle\sum_p v_p^{(i)} \right) -
 *              \left( X_m \displaystyle\sum_n Q_n v_n^{(i)} -
 *                  X_m \Pi \displaystyle\sum_p v_p^{(i)} \right) \right], \\
 *     & = \frac{c_{xX} Q_m}{\mathrm{MW}_im_T\Pi^2} \left[
 *            \Pi v_m^{(i)} - X_m\displaystyle\sum_n Q_n v_n^{(i)}\right],
 *        \textrm{ and} \\
 *    \frac{\partial \Psi_{mn}}{\partial c_i} & = 0
 *   \end{align*}
 * \f]
 *  The partial derivative of the group residual activity coefficient (Eq. 8)
 *  with respect to \f$c_i\f$ is:
 * \f[
 *   \begin{align*}
 *    \frac{\partial \ln{\Gamma_k}}{\partial c_i} & = Q_k \left[
 *      - \frac{1}{\displaystyle\sum_m\Theta_m\Psi_{mk}}
 *        \displaystyle\sum_m\frac{\partial\Theta_m}{\partial c_i}\Psi_{mk} \\
 *      \quad - \displaystyle\sum_m\left(
 *        \frac{\partial\Theta_m}{\partial c_i}\Psi_{km}
 *          \displaystyle\sum_n\Theta_n\Psi_{nm}
 *      - \Theta_m\Psi_{km}\displaystyle\sum_n
 *          \frac{\partial\Theta_n}{\partial c_i}
 *          \Psi_{nm}\right)/\left(
 *            \displaystyle\sum_n\Theta_n\Psi_{nm}\right)^2
 *      \right]
 *   \end{align*}
 * \f]
 *  After some rearranging, this becomes:
 * \f[
 *   \begin{align*}
 *    \Xi_m & \equiv \displaystyle\sum_n\Theta_n\Psi_{nm}, \\
 *    \frac{\partial \ln{\Gamma_k}}{\partial c_i} & =
 *      - Q_k \displaystyle\sum_m \left(
 *      \frac{\Psi_{mk}}{\Xi_k}\frac{\partial\Theta_m}{\partial c_i}
 *      + \frac{\Psi_{km}}{\Xi_m}
 *        \frac{\partial\Theta_m}{\partial c_i}
 *      - \frac{\Theta_m\Psi_{km}}
 *             {\Xi_m^2}
 *        \displaystyle\sum_n\frac{\partial\Theta_n}{\partial c_i}\Psi_{nm}
 *    \right)
 *   \end{align*}
 * \f]
 *  The left side of the three bracketed terms in the above equation are
 *  independent of species \f$i\f$ and can be calculated outside of the loop
 *  over the independent species. The partial derivative of the full residual
 *  term with respect to \f$c_i\f$ is:
 * \f[
 *   \begin{align*}
 *   \frac{\partial \ln{\Gamma_k^{(j)}}}{\partial c_i} & = 0 \\
 *   \frac{\partial \ln{\gamma_j^R}}{\partial c_i} & =
 *     \displaystyle\sum_k v_k^{(j)}
 *        \frac{\partial \ln{\Gamma_k}}{\partial c_i}
 *   \end{align*}
 * \f]
 *  The overall equation for the partial derivative of \f$\gamma_j\f$ with
 *  respect to species \f$i\f$ is:
 * \f[
 *    \frac{\partial\gamma_j}{\partial c_i} = \gamma_j \left(
 *      \frac{\partial\ln{\gamma_j^C}}{\partial c_i}
 *      + \frac{\partial\ln{\gamma_j^R}}{\partial c_i} \right)
 * \f]
 *
 * \param sub_model_int_data Pointer to the sub model integer data
 * \param sub_model_float_data Pointer to the sub model floating-point data
 * \param sub_model_env_data Pointer to the sub model environment-dependent data
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current time step [s]
 */
#ifdef PMC_USE_SUNDIALS
void sub_model_UNIFAC_get_jac_contrib(int *sub_model_int_data,
    double *sub_model_float_data, double *sub_model_env_data,
    ModelData *model_data, realtype *J, double time_step)
{
  int *int_data = sub_model_int_data;
  double *float_data = sub_model_float_data;
  double *state    = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Loop through each instance of each phase to calculate activity
  for (int i_phase=0; i_phase<NUM_UNIQUE_PHASE_; ++i_phase) {
    for (int i_instance=0; i_instance<NUM_PHASE_INSTANCE_(i_phase);
        ++i_instance) {

      // Get the total number of moles of species in this phase instance
      // (in umoles)
      double m_T = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
        m_T += state[PHASE_INST_ID_(i_phase, i_instance)
          + SPEC_ID_(i_phase, i_spec)] / MW_I_(i_phase, i_spec);
      }

      // If there are no species present, skip the calculations
      if (m_T<=0.0) {
        continue;
      }

      // Pre-calculate mole fractions and variables that do not depend on
      // dependent species j or independent species i

      double sigma = 0.0;
      double tau   = 0.0;
      double mu    = 0.0;
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {

        // Mole fractions x_i for this phase instance
        double x_i = state[PHASE_INST_ID_(i_phase, i_instance) +
                           SPEC_ID_(i_phase, i_spec)]
                     / MW_I_(i_phase, i_spec) / m_T;
        X_I_(i_phase, i_spec) = x_i;

        // Combinatorial variables
        sigma += R_I_(i_phase, i_spec) * x_i;
        tau   += Q_I_(i_phase, i_spec) * x_i;
        mu    += L_I_(i_phase, i_spec) * x_i;
      }

      // Residual variables
      double c_xX = 0.0;
      double Pi   = 0.0;
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        X_K_(i_group)  = 0.0;
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec) {
          double v_ik_x_i = V_IK_(i_phase, i_spec, i_group)
                            * X_I_(i_phase, i_spec);
          c_xX           += v_ik_x_i;
          X_K_(i_group)  += v_ik_x_i;
        }
      }
      c_xX = 1.0 / c_xX;
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        X_K_(i_group)    *= c_xX;
        Pi               += Q_K_(i_group) * X_K_(i_group);
      }

      // Group surface area fractions Theta_m (Eq. 9)
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group)
        THETA_M_(i_group) = Q_K_(i_group) * X_K_(i_group) / Pi;

      // Xi_m
      for (int i_group=0; i_group<NUM_GROUP_; ++i_group) {
        XI_M_(i_group) = 0;
        for (int j_group=0; j_group<NUM_GROUP_; ++j_group)
          XI_M_(i_group) += THETA_M_(j_group) * PSI_MN_(j_group, i_group);
      }

      // Group residual activity coefficients ln(Gamma_k) (Eq. 8)
      for (int k=0; k<NUM_GROUP_; ++k) {
        LN_GAMMA_K_(k) = 1.0 - log( XI_M_(k) );
        for (int m=0; m<NUM_GROUP_; ++m)
          LN_GAMMA_K_(k) -= THETA_M_(m) * PSI_MN_(k,m) / XI_M_(m);
        LN_GAMMA_K_(k) *= Q_K_(k);
      }

      // Loop over dependent species j
      for (int j=0; j<NUM_SPEC_(i_phase); ++j) {

        // Mole fraction for species j
        double x_j = X_I_(i_phase, j);

        // Phi_j (Eq. 4)
        double Phi_j;
        Phi_j = R_I_(i_phase, j) * x_j / sigma;

        // Theta_j (Eq. 4)
        double Theta_j;
        Theta_j = Q_I_(i_phase, j) * x_j / tau;

        // Full combinatorial term
        double lngammaC_j = 0.0;
        if (X_I_(i_phase, j) > 0.0) {
          lngammaC_j = log( Phi_j / x_j )
                      + 5.0 * Q_I_(i_phase, j) * log( Theta_j / Phi_j )
                      + L_I_(i_phase, j)
                      - Phi_j / x_j * mu;
        }

        // Full residual term
        double lngammaR_j = 0.0;
        for (int k=0; k<NUM_GROUP_; ++k) {
          lngammaR_j += V_IK_(i_phase, j, k) *
                        ( LN_GAMMA_K_(k) - LN_GAMMA_IK_(i_phase, j, k) );
        }

        // Activity coefficient gamma_j (in m3/ug)
        double gamma_j = exp(lngammaC_j + lngammaR_j);

        // Loop over independent species i
        double dx_j_dc_i, dPhi_j_dc_i, dTheta_j_dc_i; // combinatorial
        double dmu_dc_i, dlngammaC_j_dc_i;            //   terms
        for (int i=0; i<NUM_SPEC_(i_phase); ++i) {

          // combinatorial partial derivatives
          if (i==j) {
            dx_j_dc_i     = 1.0   - x_j;
            dPhi_j_dc_i   = sigma - x_j * R_I_(i_phase, i);
            dTheta_j_dc_i = tau   - x_j * Q_I_(i_phase, i);
          } else {
            dx_j_dc_i     = - x_j;
            dPhi_j_dc_i   = - x_j * R_I_(i_phase, i);
            dTheta_j_dc_i = - x_j * Q_I_(i_phase, i);
          }
          double MWimT_inv = 1.0 / (MW_I_(i_phase, i) * m_T);
          dx_j_dc_i     *= MWimT_inv;
          dPhi_j_dc_i   *= R_I_(i_phase, j) * MWimT_inv / sigma / sigma;
          dTheta_j_dc_i *= Q_I_(i_phase, j) * MWimT_inv / tau   / tau;

          dmu_dc_i = MWimT_inv * (L_I_(i_phase, i) - mu);

          // partial derivative of full combinatorial term ln(gammaC_j)
          dlngammaC_j_dc_i = dPhi_j_dc_i / Phi_j - dx_j_dc_i / x_j
                            + 5.0 * Q_I_(i_phase, j) *
                                (dTheta_j_dc_i / Theta_j - dPhi_j_dc_i / Phi_j)
                            - dPhi_j_dc_i * mu / x_j
                            - Phi_j / x_j * dmu_dc_i
                            + Phi_j * mu / x_j / x_j * dx_j_dc_i;


          // Residual partial derivatives

          // partial derivative of group surface area fractions Theta_m
          double sumQnvn_i = 0.0;
          for (int n=0; n<NUM_GROUP_; ++n)
            sumQnvn_i += Q_K_(n) * V_IK_(i_phase, i, n);
          for (int m=0; m<NUM_GROUP_; ++m) {
            DTHETA_M_DC_I_(m) = Q_K_(m) * c_xX * MWimT_inv / Pi / Pi *
                                ( Pi * V_IK_(i_phase, i, m)
                                  - X_K_(m) * sumQnvn_i );
          }


          // partial derivative of full residual terms ln(gammaR_j)
          double dlngammaR_j_dc_i = 0.0;
          for (int k=0; k<NUM_GROUP_; ++k) {

            // partial derivative of group residual activity coefficients
            // ln(Gamma_k)
            double dlnGamma_k_dc_i = 0.0;
            for (int m=0; m<NUM_GROUP_; ++m) {
              double sum_n = 0.0;
              for (int n=0; n<NUM_GROUP_; ++n)
                sum_n += DTHETA_M_DC_I_(n) * PSI_MN_(n,m);
              dlnGamma_k_dc_i += DTHETA_M_DC_I_(m) *
                                 ( PSI_MN_(m,k) / XI_M_(k)
                                   + PSI_MN_(k,m) / XI_M_(m) )
                                 - THETA_M_(m) * PSI_MN_(k,m) /
                                   pow( XI_M_(m), 2 ) * sum_n;
            }
            dlnGamma_k_dc_i *= -Q_K_(k);

            dlngammaR_j_dc_i += V_IK_(i_phase, j, k)
                                * dlnGamma_k_dc_i;
          }

          // partial derivative of activity coefficient gamma_j
          double dgamma_j_dc_i = gamma_j * x_j *
                                   (dlngammaC_j_dc_i + dlngammaR_j_dc_i)
                                 + gamma_j * dx_j_dc_i;

          // Add the partial derivative contribution to the
          // Jacobian
          J[JAC_ID_(i_phase, i_instance, j, i)] += dgamma_j_dc_i;

        }
      }
    }
  }

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
    printf("\n  Q_K: %le R_K: %le X_k: %le",
           Q_K_(i_group), R_K_(i_group), X_K_(i_group));
    printf("\n  Xi_M: %le ln(Gamma_k): %le", XI_M_(i_group),
           LN_GAMMA_K_(i_group));
    printf("\n  dTheta_n / dc_i): %le", DTHETA_M_DC_I_(i_group));
    printf("\n A_MN (by group):");
    for (int j_group=0; j_group<NUM_GROUP_; ++j_group)
      printf(" %le", A_MN_(i_group, j_group));
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
      printf("\n    State id: %d", PHASE_INST_ID_(i_phase, i_inst));
      printf("\n    Gamma ids:");
      for (int i_gamma=0; i_gamma<NUM_SPEC_(i_phase); ++i_gamma)
        printf(" %d", PHASE_INST_ID_(i_phase, i_inst) +
                      GAMMA_ID_(i_phase, i_gamma));
      printf("\n    Species ids:");
      for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
        printf(" %d", PHASE_INST_ID_(i_phase, i_inst) +
                      SPEC_ID_(i_phase, i_spec));
      printf("\n    Jacobian ids:");
      for (int i_gamma=0; i_gamma<NUM_SPEC_(i_phase); ++i_gamma)
        for (int i_spec=0; i_spec<NUM_SPEC_(i_phase); ++i_spec)
          printf(" J[%d][%d] =  %d", i_gamma, i_spec,
                 JAC_ID_(i_phase, i_inst, i_gamma, i_spec));

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
    }
  }
}
