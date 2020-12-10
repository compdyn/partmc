/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Jacobian functions
 *
 */
/** \file
 * \brief Jacobian functions
 */
#include "Jacobian.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "time_derivative.h"

#define SMALL_NUMBER 1e-90

int jacobian_initialize(Jacobian *jac, unsigned int num_spec,
                        unsigned int **jac_struct) {
  unsigned int num_elem = 0;
  for (unsigned int i_col = 0; i_col < num_spec; ++i_col)
    for (unsigned int i_row = 0; i_row < num_spec; ++i_row)
      if (jac_struct[i_col][i_row] == 1) ++num_elem;

  jac->num_spec = num_spec;
  jac->num_elem = num_elem;
  jac->col_ptrs = (unsigned int *)malloc((num_spec + 1) * sizeof(unsigned int));
  if (!jac->col_ptrs) return 0;
  jac->row_ids = (unsigned int *)malloc(num_elem * sizeof(unsigned int));
  if (!jac->row_ids) {
    free(jac->col_ptrs);
    return 0;
  }
  jac->production_partials =
      (long double *)malloc(num_elem * sizeof(long double));
  if (!jac->production_partials) {
    free(jac->col_ptrs);
    free(jac->row_ids);
    return 0;
  }
  jac->loss_partials = (long double *)malloc(num_elem * sizeof(long double));
  if (!jac->loss_partials) {
    free(jac->col_ptrs);
    free(jac->row_ids);
    free(jac->production_partials);
    return 0;
  }

  unsigned int i_elem = 0;
  unsigned int i_col = 0;
  for (; i_col < num_spec; ++i_col) {
    jac->col_ptrs[i_col] = i_elem;
    for (unsigned int i_row = 0; i_row < num_spec; ++i_row) {
      if (jac_struct[i_row][i_col] == 1) jac->row_ids[i_elem++] = i_row;
    }
  }
  jac->col_ptrs[i_col] = i_elem;
  return 1;
}

unsigned int jacobian_get_element_id(Jacobian jac, unsigned int col_id,
                                     unsigned int row_id) {
  if (col_id >= jac.num_spec || col_id < 0) {
    printf(
        "\nError: Bad Jacobian column id: %u. Expected value between 0 and "
        "%u\n",
        col_id, jac.num_spec);
    exit(EXIT_FAILURE);
  }

  for (unsigned int i_elem = jac.col_ptrs[col_id];
       i_elem < jac.col_ptrs[col_id + 1]; ++i_elem) {
    if (jac.row_ids[i_elem] == row_id) return i_elem;
  }

  printf("\nError: Invalid Jacobian element specified: %u %u\n", col_id,
         row_id);
  exit(EXIT_FAILURE);
  return 0;
}

void jacobian_reset(Jacobian jac) {
  for (unsigned int i_elem = 0; i_elem < jac.num_elem; ++i_elem) {
    jac.production_partials[i_elem] = 0.0;
    jac.loss_partials[i_elem] = 0.0;
  }
}

void jacobian_output(Jacobian jac, double *dest_array) {
  for (unsigned int i_col = 0; i_col < jac.num_spec; ++i_col) {
    for (unsigned int i_elem = jac.col_ptrs[i_col];
         i_elem < jac.col_ptrs[i_col + 1]; ++i_elem) {
      long double drf_dy = jac.production_partials[i_elem];
      long double drr_dy = jac.loss_partials[i_elem];
      dest_array[i_elem] = drf_dy - drr_dy;
    }
  }
}

void jacobian_add_value(Jacobian jac, unsigned int elem_id,
                        unsigned int prod_or_loss,
                        long double jac_contribution) {
  if (prod_or_loss == JACOBIAN_PRODUCTION)
    jac.production_partials[elem_id] += jac_contribution;
  if (prod_or_loss == JACOBIAN_LOSS)
    jac.loss_partials[elem_id] += jac_contribution;
}

void jacobian_free(Jacobian jac) {
  free(jac.col_ptrs);
  free(jac.row_ids);
  free(jac.production_partials);
  free(jac.loss_partials);
}
