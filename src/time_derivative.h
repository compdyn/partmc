/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header for the time derivative structure and related functions
 *
 */
/** \file
 * \brief Header for the time derivative structure and related functions
 */
#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include <math.h>
#include <stdlib.h>

// Threshhold for precisition loss in rate calculations
#define MAX_PRECISION_LOSS 1.0e-14

/* Time derivative for solver species */
typedef struct {
  unsigned int num_spec;          // Number of species in the derivative
  long double *production_rates;  // Production rates for all species
  long double *loss_rates;        // Loss rates for all species
#ifdef PMC_DEBUG
  double last_max_loss_precision;  // Maximum loss of precision at last output
#endif
} TimeDerivative;

/** \brief Initialize the derivative
 *
 * \param time_deriv Pointer to the TimeDerivative object
 * \param num_spec Number of species to include in the derivative
 * \return Flag indicating whether the derivative was sucessfully initialized
 *         (0 = false; 1 = true)
 */
int time_derivative_initialize(TimeDerivative *time_deriv,
                               unsigned int num_spec);

/** \brief Reset the derivative
 *
 * \param time_deriv TimeDerivative object
 */
void time_derivative_reset(TimeDerivative time_deriv);

/** \brief Output the current derivative array
 *
 * \param time_deriv TimeDerivative object
 * \param dest_array Pointer to the destination array
 * \param deriv_est Pointer to an estimate of the derivative array (optional)
 * \param output_precision Output the estimated loss of precision for each
 * species if output_precision == 1
 */
void time_derivative_output(TimeDerivative time_deriv, double *dest_array,
                            double *deriv_est, unsigned int output_precision);

/** \brief Add a contribution to the time derivative
 *
 * \param time_deriv TimeDerivative object
 * \param spec_id Index of the species to update rates for
 * \param rate_contribution Value to add to the time derivative for speces
 * spec_id
 */
void time_derivative_add_value(TimeDerivative time_deriv, unsigned int spec_id,
                               long double rate_contribution);

#ifdef PMC_DEBUG
/** \brief Maximum loss of precision at the last output of the derivative
 *         in bits
 *
 * \param time_deriv TimeDerivative object
 * \return maximum loss of precision
 */
double time_derivative_max_loss_precision(TimeDerivative time_deriv);
#endif

/** \brief Free memory associated with a TimeDerivative
 *
 * \param time_deriv TimeDerivative object
 */
void time_derivative_free(TimeDerivative time_deriv);

#endif
