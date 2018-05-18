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
#ifdef PMC_USE_SUNDIALS

#include "../sub_model_solver.h"

// TODO Lookup environmental indicies during initialization
#define _TEMPERATURE_K_ env_data[0]
#define _PRESSURE_PA_ env_data[1]

#define _NUM_UNIQUE_PHASE_ (int_data[0])
#define _NUM_GROUP_ (int_data[1])
#define _TOTAL_INT_PROP_ (int_data[2])
#define _TOTAL_FLOAT_PROP_ (int_data[3])
#define _NUM_INT_PROP_ 4
#define _NUM_FLOAT_PROP_ 0
#define _PHASE_INT_LOC_(p) (int_data[_NUM_INT_PROP_+p]-1)
#define _PHASE_FLOAT_LOC_(p) (int_data[_NUM_INT_PROP_+_NUM_UNIQUE_PHASE_+p]-1)
#define _NUM_PHASE_INSTANCE_(p) (int_data[_PHASE_INT_LOC_(p)])
#define _NUM_SPEC_(p) (int_data[_PHASE_INT_LOC_(p)+1])
#define _PHASE_INST_FLOAT_LOC_(p,c) (int_data[_PHASE_INT_LOC_(p)+2+c]-1)
#define _PHASE_INST_ID_(p,c) (int_data[_PHASE_INT_LOC_(p)+2+_NUM_PHASE_INSTANCE_(p)+c]-1)
#define _SPEC_ID_(p,i) (int_data[_PHASE_INT_LOC_(p)+2+2*_NUM_PHASE_INSTANCE_(p)+i])
#define _v_ik_(p,i,k) (int_data[_PHASE_INT_LOC_(p)+2+2*_NUM_PHASE_INSTANCE_(p)+(k+1)*_NUM_SPEC_(p)+i])

#define _Q_k_(k) (float_data[k])
#define _R_k_(k) (float_data[_NUM_GROUP_+k])
#define _THETA_m_(m) (float_data[2*_NUM_GROUP_+m])
#define _a_mn_(m,n) (float_data[(m+3)*_NUM_GROUP_+n])
#define _PSI_mn_(m,n) (float_data[(m+3+_NUM_GROUP_)*_NUM_GROUP_+n])
#define _r_i_(p,i) (float_data[_PHASE_FLOAT_LOC_(p)+i])
#define _q_i_(p,i) (float_data[_PHASE_FLOAT_LOC_(p)+_NUM_SPEC_(p)+i])
#define _l_i_(p,i) (float_data[_PHASE_FLOAT_LOC_(p)+2*_NUM_SPEC_(p)+i])
#define _MW_i_(p,i) (float_data[_PHASE_FLOAT_LOC_(p)+3*_NUM_SPEC_(p)+i])
#define _X_i_(p,i) (float_data[_PHASE_FLOAT_LOC_(p)+4*_NUM_SPEC_(p)+i])
#define _ln_GAMMA_ik_(p,i,k) (float_data[_PHASE_FLOAT_LOC_(p)+i*_NUM_GROUP_+5*_NUM_SPEC_(p)+k])
#define _gamma_i_(p,c,i) (float_data[_PHASE_INST_FLOAT_LOC_(p,c)+i])

#define _INT_DATA_SIZE_ (_TOTAL_INT_PROP_)
#define _FLOAT_DATA_SIZE_ (_TOTAL_FLOAT_PROP_)

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

/** \brief Get the id of a parameter in the condensed data block
 *
 * \param sub_mode_data Pointer to the sub-model data
 * \param identifiers For the UNIFAC model, the identifers are just the id
 *                    on the state array of the aerosol-phase species for
 *                    which the acitivty is needed.
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_get_parameter_id(void *sub_model_data, void *identifiers,
    int *parameter_id)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  for (int i_phase=0; i_phase<_NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_instance=0; i_instance<_NUM_PHASE_INSTANCE_(i_phase);
        i_instance++) {
      for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
        if (*((int*)identifiers) == _SPEC_ID_(i_phase, i_spec) 
            + _PHASE_INST_ID_(i_phase, i_instance)) {
          *parameter_id = (int) (((int*) 
                (&(_gamma_i_(i_phase, i_instance, i_spec)))) - int_data);
          return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
        }
      }
    }
  }
  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

/** \brief Update sub-model data for new environmental conditions
 *
 * \param sub_model_data Pointer to the sub-model data
 * \param env_data Pointer to the environmental state array
 * \return The sub_model_data pointer advanced by the size of the sub model
 */
void * sub_model_UNIFAC_update_env_state(void *sub_model_data, realtype *env_data)
{
  int *int_data = (int*) sub_model_data;
  realtype *float_data = (realtype*) &(int_data[_INT_DATA_SIZE_]);

  // Update the interaction parameters
  for (int m=0; m<_NUM_GROUP_; m++)
    for (int n=0; n<_NUM_GROUP_; n++)
      _PSI_mn_(m,n) = exp(-_a_mn_(m,n)/_TEMPERATURE_K_);

  // Calculate the pure liquid residual acitivity ceofficient ln(GAMMA_k^(i))
  // terms. Eq. 7 & 8
  for (int i_phase=0; i_phase<_NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
      
      // Calculate the sum (Q_n * X_n) in the denominator of Eq. 9 for the
      // pure liquid
      realtype sum_Qn_Xn = 0.0;
      for (int i_group=0; i_group<_NUM_GROUP_; i_group++)
        sum_Qn_Xn += _Q_k_(i_group) * _v_ik_(i_phase, i_spec, i_group);

      // Calculate THETA_m Eq. 9
      for (int m=0; m<_NUM_GROUP_; m++)
        _THETA_m_(m) = _Q_k_(m) * _v_ik_(i_phase, i_spec, m) / sum_Qn_Xn;

      // Calculate ln(GAMMA_k^(i))
      for (int k=0; k<_NUM_GROUP_; k++) {
        realtype sum_m_A = 0.0; // ln( sum_m_A ) term in Eq. 8
        realtype sum_m_B = 0.0; // last term in Eq. 8
        for (int m=0; m<_NUM_GROUP_; m++) {
          sum_m_A += _THETA_m_(m) * _PSI_mn_(m,k);
          realtype sum_n = 0.0;
          for (int n=0; n<_NUM_GROUP_; n++)
            sum_n += _THETA_m_(n) * _PSI_mn_(n,m);
          sum_m_B += _THETA_m_(m) * _PSI_mn_(k,m) / sum_n;
        }
        _ln_GAMMA_ik_(i_phase, i_spec, k) = _Q_k_(k) * 
                  (1.0 - log(sum_m_A) - sum_m_B);
      }
    } 
  }

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

  // Loop through each instance of each phase to calculate activity
  for (int i_phase=0; i_phase<_NUM_UNIQUE_PHASE_; i_phase++) {
    for (int i_instance=0; i_instance<_NUM_PHASE_INSTANCE_(i_phase);
        i_instance++) {

      // Get the total number of moles of species in this phase instance
      realtype total_umoles = 0.0;
      for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
        total_umoles += model_data->state[_PHASE_INST_ID_(i_phase, i_instance)
          + _SPEC_ID_(i_phase, i_spec)] / _MW_i_(i_phase, i_spec);
      }

      // If there are no species present, skip the calculations
      if (total_umoles<=0.0) {
        for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
          _gamma_i_(i_phase, i_instance, i_spec) = 1.0;
        }
        continue;
      }
      
      // Update the mole fractions X_i
      for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
        _X_i_(i_phase, i_spec) = model_data->state[
          _PHASE_INST_ID_(i_phase, i_instance) + _SPEC_ID_(i_phase, i_spec)]
          / _MW_i_(i_phase, i_spec) / total_umoles;
      }

      // Calcualte the sum (Q_n * X_n) in the denominator of Eq. 9 for the
      // mixture
      realtype sum_Qn_Xn_mixture = 0.0;
      for (int n=0; n<_NUM_GROUP_; n++) {
        for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
          sum_Qn_Xn_mixture += _Q_k_(n) * _X_i_(i_phase, i_spec) 
            * _v_ik_(i_phase, i_spec, n);
        }
      }

      // Calculate the group surface area fractions THETA_m (Eq. 9)
      for (int m=0; m<_NUM_GROUP_; m++) {
        _THETA_m_(m) = 0.0;
        for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {
          _THETA_m_(m) += _Q_k_(m) * _X_i_(i_phase, i_spec)
            * _v_ik_(i_phase, i_spec, m);
        }
        _THETA_m_(m) /= sum_Qn_Xn_mixture;
      }

      // Calculate activity coefficients for each species in the phase instance
      for (int i_spec=0; i_spec<_NUM_SPEC_(i_phase); i_spec++) {

        // Calculate PHI_i (Eq. 4)
        realtype PHI_i = 0.0;
        for (int j_spec=0; j_spec<_NUM_SPEC_(i_phase); j_spec++)
          PHI_i += _r_i_(i_phase, j_spec) * _X_i_(i_phase, j_spec);
        PHI_i = _r_i_(i_phase, i_spec) * _X_i_(i_phase, i_spec) / PHI_i;

        // Calculate THETA_i (Eq. 4)
        realtype THETA_i = 0.0;
        for (int j_spec=0; j_spec<_NUM_SPEC_(i_phase); j_spec++)
          THETA_i += _q_i_(i_phase, j_spec) * _X_i_(i_phase, j_spec);
        THETA_i = _q_i_(i_phase, i_spec) * _X_i_(i_phase, i_spec) / THETA_i;

        // Calculate the combinatorial term ln(gamma_i^C) Eq. 3
        realtype ln_gamma_i_C = 0.0;
        if (_X_i_(i_phase, i_spec) > 0.0) {
          for (int j_spec=0; j_spec<_NUM_SPEC_(i_phase); j_spec++)
            ln_gamma_i_C += _l_i_(i_phase, j_spec) * _X_i_(i_phase, j_spec);
          ln_gamma_i_C = log( PHI_i / _X_i_(i_phase, i_spec) )
                       + 5.0 * _q_i_(i_phase, i_spec) * log( THETA_i / PHI_i )
                       + _l_i_(i_phase, i_spec)
                       - PHI_i / _X_i_(i_phase, i_spec) * ln_gamma_i_C;
        }

        // Calculate the residual term ln(gamma_i^R) Eq. 7
        realtype ln_gamma_i_R = 0.0;
        for (int k=0; k<_NUM_GROUP_; k++) {
          realtype sum_m_A = 0.0; // ln( sum_m_A ) term in Eq. 8
          realtype sum_m_B = 0.0; // last term in Eq. 8
          for (int m=0; m<_NUM_GROUP_; m++) {
            sum_m_A += _THETA_m_(m) * _PSI_mn_(m,k);
            realtype sum_n = 0.0;
            for (int n=0; n<_NUM_GROUP_; n++)
              sum_n += _THETA_m_(n) * _PSI_mn_(n,m);
            sum_m_B += _THETA_m_(m) * _PSI_mn_(k,m) / sum_n;
          }
          // calculate ln(GAMMA_k) Eq. 8
          realtype ln_GAMMA_k = _Q_k_(k) * (1.0 - log(sum_m_A) - sum_m_B);
          ln_gamma_i_R += _v_ik_(i_phase, i_spec, k) 
                       * ( ln_GAMMA_k - _ln_GAMMA_ik_(i_phase, i_spec, k));
        }

        // Calculate gamma_i Eq. 1 and convert to units of (m^3/ug)
        _gamma_i_(i_phase, i_instance, i_spec) = 
          exp( ln_gamma_i_C + ln_gamma_i_R ) / total_umoles 
          / _MW_i_(i_phase, i_spec);
 
      }
    }
  }

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

/** \brief Skip through the sub model, only advancing the data pointer
 *
 * \param sub_model_data Pointer to the sub model data
 * \return The sub_model_data pointer advanced by the size of the sub-model
 */
void * sub_model_UNIFAC_skip(void *sub_model_data)
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

  printf("\n\nUNIFAC sub model\n\n");
  printf("\nint_data");
  for (int i=0; i<_INT_DATA_SIZE_; i++) printf(" %d", int_data[i]);
  printf("\nfloat_data");
  for (int i=0; i<_FLOAT_DATA_SIZE_; i++) printf(" %le", float_data[i]);

  return (void*) &(float_data[_FLOAT_DATA_SIZE_]);
}

#undef _TEMPERATURE_K_
#undef _PRESSURE_PA_

#undef _NUM_UNIQUE_PHASE_
#undef _NUM_GROUP_
#undef _TOTAL_INT_PROP_
#undef _TOTAL_FLOAT_PROP_
#undef _NUM_INT_PROP_
#undef _NUM_FLOAT_PROP_
#undef _PHASE_INT_LOC_
#undef _PHASE_FLOAT_LOC_
#undef _NUM_PHASE_INSTANCE_
#undef _NUM_SPEC_
#undef _PHASE_INST_FLOAT_LOC_
#undef _PHASE_INST_ID_
#undef _SPEC_ID_
#undef _v_ik_

#undef _Q_k_
#undef _R_k_
#undef _THETA_m_
#undef _a_mn_
#undef _PSI_mn_
#undef _r_i_
#undef _q_i_
#undef _l_i_
#undef _MW_i_
#undef _X_i_
#undef _ln_GAMMA_ik_
#undef _gamma_i_

#undef _INT_DATA_SIZE_
#undef _FLOAT_DATA_SIZE_

#endif
