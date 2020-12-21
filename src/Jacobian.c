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

#define BUFFER_SIZE 10
#define SMALL_NUMBER 1e-90

int jacobian_initialize_empty(Jacobian *jac, unsigned int num_spec) {
  jac->num_spec = num_spec;
  jac->num_elem = 0;
  jac->elements = (JacobianColumnElements *)malloc(
      num_spec * sizeof(JacobianColumnElements));
  if (!jac->elements) {
    jacobian_free(jac);
    return 0;
  }
  for (unsigned int i_col = 0; i_col < num_spec; ++i_col) {
    jac->elements[i_col].array_size = BUFFER_SIZE;
    jac->elements[i_col].number_of_elements = 0;
    jac->elements[i_col].row_ids =
        (unsigned int *)malloc(BUFFER_SIZE * sizeof(unsigned int));
    if (!jac->elements[i_col].row_ids) {
      jacobian_free(jac);
      return 0;
    }
  }
  jac->col_ptrs = NULL;
  jac->row_ids = NULL;
  jac->production_partials = NULL;
  jac->loss_partials = NULL;
  return 1;
}

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
  jac->elements = NULL;
  return 1;
}

// Add buffer space for Jacobian column elements (returns 1 on success, 0
// otherwise)
int jacobian_column_elements_add_space(JacobianColumnElements *column) {
  unsigned int *temp_ids;
  column->array_size += BUFFER_SIZE;
  temp_ids = (unsigned int *)malloc(column->array_size * sizeof(unsigned int));
  if (!temp_ids) return 0;
  for (unsigned int i_elem = 0; i_elem < column->number_of_elements; ++i_elem) {
    temp_ids[i_elem] = column->row_ids[i_elem];
  }
  free(column->row_ids);
  column->row_ids = temp_ids;
  return 1;
}

void jacobian_register_element(Jacobian *jac, unsigned int dep_id,
                               unsigned int ind_id) {
  if (!jac->elements) {
    printf(
        "\n\nERROR - Trying to register elements in a Jacobian that has "
        "already been built.\n\n");
    exit(EXIT_FAILURE);
  }
  JacobianColumnElements *column = &(jac->elements[ind_id]);
  for (unsigned int i_elem = 0; i_elem < column->number_of_elements; ++i_elem) {
    if (column->row_ids[i_elem] == dep_id) return;
  }
  if (column->array_size == column->number_of_elements) {
    jacobian_column_elements_add_space(column);
  }
  jac->elements[ind_id].row_ids[jac->elements[ind_id].number_of_elements] =
      dep_id;
  ++(jac->elements[ind_id].number_of_elements);
}

int compare_ids(const void *a, const void *b) {
  return (*(unsigned int *)a - *(unsigned int *)b);
}

unsigned int jacobian_build_matrix(Jacobian *jac) {
  if (!jac->elements) {
    printf(
        "\n\nERROR - Trying to build a Jacobian that has already been "
        "built.\n\n");
    exit(EXIT_FAILURE);
  }
  jac->num_elem = 0;
  for (unsigned int i_col = 0; i_col < jac->num_spec; ++i_col) {
    JacobianColumnElements *column = &(jac->elements[i_col]);
    qsort(column->row_ids, column->number_of_elements, sizeof(unsigned int),
          compare_ids);
    jac->num_elem += column->number_of_elements;
  }
  jac->col_ptrs =
      (unsigned int *)malloc((jac->num_spec + 1) * sizeof(unsigned int));
  if (!jac->col_ptrs) {
    jacobian_free(jac);
    return 0;
  }
  jac->row_ids = (unsigned int *)malloc(jac->num_elem * sizeof(unsigned int));
  if (!jac->row_ids) {
    jacobian_free(jac);
    return 0;
  }
  unsigned int i_elem = 0;
  for (unsigned int i_col = 0; i_col < jac->num_spec; ++i_col) {
    jac->col_ptrs[i_col] = i_elem;
    JacobianColumnElements *column = &(jac->elements[i_col]);
    for (unsigned int j_elem = 0; j_elem < column->number_of_elements;
         ++j_elem) {
      jac->row_ids[i_elem++] = column->row_ids[j_elem];
    }
  }
  jac->col_ptrs[jac->num_spec] = i_elem;
  if (i_elem != jac->num_elem) {
    printf("\n\n ERROR - Internal error building Jacobain matrix %d %d\n\n",
           i_elem, jac->num_elem);
    exit(EXIT_FAILURE);
  }
  jac->production_partials =
      (long double *)malloc(jac->num_elem * sizeof(long double));
  if (!jac->production_partials) {
    jacobian_free(jac);
    return 0;
  }
  jac->loss_partials =
      (long double *)malloc(jac->num_elem * sizeof(long double));
  if (!jac->loss_partials) {
    jacobian_free(jac);
    return 0;
  }
  for (unsigned int i_col = 0; i_col < jac->num_spec; ++i_col) {
    jacobian_column_elements_free(&(jac->elements[i_col]));
  }
  free(jac->elements);
  jac->elements = NULL;
  jacobian_reset(*jac);
  return 1;
}

unsigned int jacobian_number_of_elements(Jacobian jac) { return jac.num_elem; }

unsigned int jacobian_column_pointer_value(Jacobian jac, unsigned int col_id) {
  return jac.col_ptrs[col_id];
}

unsigned int jacobian_row_index(Jacobian jac, unsigned int elem_id) {
  return jac.row_ids[elem_id];
}

unsigned int jacobian_get_element_id(Jacobian jac, unsigned int dep_id,
                                     unsigned int ind_id) {
  if (ind_id >= jac.num_spec || ind_id < 0) {
    printf(
        "\nError: Bad Jacobian column id: %u. Expected value between 0 and "
        "%u\n",
        ind_id, jac.num_spec);
    exit(EXIT_FAILURE);
  }
  for (unsigned int i_elem = jac.col_ptrs[ind_id];
       i_elem < jac.col_ptrs[ind_id + 1]; ++i_elem) {
    if (jac.row_ids[i_elem] == dep_id) return i_elem;
  }
  return -1;
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

void jacobian_print(Jacobian jac) {
  printf("\n   *** Jacobian ***");
  printf("\nnumber of variables: %d", jac.num_spec);
  printf("\nnumber of non-zero Jacobian elements: %d", jac.num_elem);
  if (!jac.col_ptrs && !jac.row_ids && !jac.production_partials &&
      !jac.loss_partials && jac.elements) {
    printf("\nstatus: building Jacobian");
    for (unsigned int i_col = 0; i_col < jac.num_spec; ++i_col) {
      for (unsigned int i_elem = 0;
           i_elem < jac.elements[i_col].number_of_elements; ++i_elem) {
        printf("\n  col = %6d row = %6d", i_col,
               jac.elements[i_col].row_ids[i_elem]);
      }
    }
  } else if (jac.col_ptrs && jac.row_ids && jac.production_partials &&
             jac.loss_partials) {
    printf("\nstatus: Jacobian built");
    for (unsigned int i_col = 0; i_col < jac.num_spec; ++i_col) {
      for (unsigned int i_elem = jac.col_ptrs[i_col];
           i_elem < jac.col_ptrs[i_col + 1]; ++i_elem) {
        printf("\n  col = %6d row = %6d production = %Le loss = %Le", i_col,
               jac.row_ids[i_elem], jac.production_partials[i_elem],
               jac.loss_partials[i_elem]);
      }
    }
  } else {
    printf("\nstatus: invalid state");
  }
  printf("\n  *** end Jacobian ***");
}

void jacobian_column_elements_free(JacobianColumnElements *column) {
  if (column->row_ids) {
    free(column->row_ids);
    column->row_ids = NULL;
  }
}

void jacobian_free(Jacobian *jac) {
  if (jac->col_ptrs) {
    free(jac->col_ptrs);
    jac->col_ptrs = NULL;
  }
  if (jac->row_ids) {
    free(jac->row_ids);
    jac->row_ids = NULL;
  }
  if (jac->production_partials) {
    free(jac->production_partials);
    jac->production_partials = NULL;
  }
  if (jac->loss_partials) {
    free(jac->loss_partials);
    jac->loss_partials = NULL;
  }
  if (jac->elements) {
    for (unsigned int i_col = 0; i_col < jac->num_spec; ++i_col) {
      jacobian_column_elements_free(&(jac->elements[i_col]));
    }
    free(jac->elements);
    jac->elements = NULL;
  }
}
