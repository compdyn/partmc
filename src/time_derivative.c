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
#include <math.h>
#include <stdio.h>

int time_derivative_initialize(TimeDerivative *time_deriv,
                               unsigned int num_spec) {
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

#ifdef PMC_DEBUG
  time_deriv->last_max_loss_precision = 0.0;
#endif

  return 1;
}

void time_derivative_reset(TimeDerivative time_deriv) {
  for (unsigned int i_spec = 0; i_spec < time_deriv.num_spec; ++i_spec) {
    time_deriv.production_rates[i_spec] = 0.0;
    time_deriv.loss_rates[i_spec] = 0.0;
  }
}

void time_derivative_output(TimeDerivative time_deriv, double *dest_array,
                            double *deriv_est, unsigned int output_precision) {
  long double *r_p = time_deriv.production_rates;
  long double *r_l = time_deriv.loss_rates;

#ifdef PMC_DEBUG
  time_deriv.last_max_loss_precision = 1.0;
#endif

  for (unsigned int i_spec = 0; i_spec < time_deriv.num_spec; ++i_spec) {
    double prec_loss = 1.0;
    if (*r_p + *r_l != 0.0) {
      if (deriv_est) {
        long double scale_fact;
        scale_fact =
            1.0 / (*r_p + *r_l) /
            (1.0 / (*r_p + *r_l) + MAX_PRECISION_LOSS / fabsl(*r_p - *r_l));
        *dest_array =
            scale_fact * (*r_p - *r_l) + (1.0 - scale_fact) * (*deriv_est);
      } else {
        *dest_array = *r_p - *r_l;
      }
#ifdef PMC_DEBUG
      if (*r_p != 0.0 && *r_l != 0.0) {
        prec_loss = *r_p > *r_l ? 1.0 - *r_l / *r_p : 1.0 - *r_p / *r_l;
        if (prec_loss < time_deriv.last_max_loss_precision)
          time_deriv.last_max_loss_precision = prec_loss;
      }
#endif
    } else {
      *dest_array = 0.0;
    }
    ++r_p;
    ++r_l;
    ++dest_array;
    if (deriv_est) ++deriv_est;
#ifdef PMC_DEBUG
    if (output_precision == 1) {
      printf("\nspec %d prec_loss %le", i_spec, -log(prec_loss) / log(2.0));
    }
#endif
  }
}

void time_derivative_add_value(TimeDerivative time_deriv, unsigned int spec_id,
                               long double rate_contribution) {
  if (rate_contribution > 0.0) {
    time_deriv.production_rates[spec_id] += rate_contribution;
  } else {
    time_deriv.loss_rates[spec_id] += -rate_contribution;
  }
}

#ifdef PMC_DEBUG
double time_derivative_max_loss_precision(TimeDerivative time_deriv) {
  return -log(time_deriv.last_max_loss_precision) / log(2.0);
}
#endif

void time_derivative_free(TimeDerivative time_deriv) {
  free(time_deriv.production_rates);
  free(time_deriv.loss_rates);
}
