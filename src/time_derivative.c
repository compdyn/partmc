/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Functions of the time derivative structure
 *
 */
/** \file
 * \brief Functions of the time derivative structure
 */
#include "time_derivative.h"

int time_derivative_initialize(TimeDerivative *time_deriv, int num_spec) {
  if (num_spec <= 0) return 0;

  time_deriv->production_rates =
      (long double *)malloc(num_spec * sizeof(long double));
  if (time_deriv->production_rates == NULL) return 0;

  time_deriv->loss_rates =
      (long double *)malloc(num_spec * sizeof(long double));
  if (time_deriv->loss_rates == NULL) {
    free(time_deriv->production_rates);
    return 0;
  }

  time_deriv->num_spec = num_spec;

  return 1;
}

void time_derivative_reset(TimeDerivative *time_deriv) {
  for (int i_spec = 0; i_spec < time_deriv->num_spec; ++i_spec) {
    time_deriv->production_rates[i_spec] = 0.0;
    time_deriv->loss_rates[i_spec] = 0.0;
  }
}

void time_derivative_output(TimeDerivative *time_deriv, double *dest_array) {
  long double *prod_rate = time_deriv->production_rates;
  long double *loss_rate = time_deriv->loss_rates;
  long double rate, loss_est;
  for (int i_spec = 0; i_spec < time_deriv->num_spec; ++i_spec) {
    rate = *prod_rate - *loss_rate;
    if (rate == 0.0) {
      loss_est = 1.0;
    } else {
      loss_est = fabsl(rate / (*prod_rate + *loss_rate));
      loss_est /= (loss_est + MAX_PRECISION_LOSS);
    }
    *(dest_array++) = loss_est * rate;
    ++prod_rate;
    ++loss_rate;
  }
}

void time_derivative_add_value(TimeDerivative *time_deriv, int spec_id,
                               long double rate_contribution) {
  if (rate_contribution > 0.0) {
    time_deriv->production_rates[spec_id] += rate_contribution;
  } else {
    time_deriv->loss_rates[spec_id] += -rate_contribution;
  }
}

void time_derivative_free(TimeDerivative *time_deriv) {
  free(time_deriv->production_rates);
  free(time_deriv->loss_rates);
}
