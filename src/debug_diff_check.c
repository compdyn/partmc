/* Copyright (C) 2020 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Difference checker for int_data, float_data, env_data (until we can get
 * rid of these entirely in version 2.0!)
 */
/** \file
 * \brief model element data difference checker - NOT THREAD SAFE!
 */
#include "debug_diff_check.h"

/* Starting index for data for each model element (plus 1+the last value) */
typedef struct {
  int n_elements;  // number of model elements

  int *int_data;    // integer data
  int *float_data;  // float data
  int *env_data;    // environmental data

} ModelElementDataIndices;

/* Model element data pointers */
typedef struct {
  int *int_data;       // start of the integer data
  double *float_data;  // start of the float data
  double *env_data;    // start of the environmental data for each grid cell

} ModelElementDataPointers;

/* Model element data and indices */
typedef struct {
  int n_cells;  // number of grid cells

  ModelElementDataPointers ptrs;
  ModelElementDataIndices indices;

} ModelElement;

/* Model data */
typedef struct {
  ModelElement reactions;
  ModelElement aero_reps;
  ModelElement sub_models;

} Model;

/* Difference checker data */
typedef struct {
  Model current;     // current model data used by solver
  Model last_check;  // model data values at last check

} DifferenceCheckerData;

DifferenceCheckerData diff_data[2];
int num_solvers = 0;

// Allocate index arrays for a specified number of elements
void allocate_index_arrays(ModelElement *model_element, int num_elements) {
  // Set the number of elements and allocate the arrays
  model_element->indices.n_elements = num_elements;
  model_element->indices.int_data =
      (int *)malloc(sizeof(int) * (num_elements + 1));
  if (model_element->indices.int_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }
  model_element->indices.float_data =
      (int *)malloc(sizeof(int) * (num_elements + 1));
  if (model_element->indices.float_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }
  model_element->indices.env_data =
      (int *)malloc(sizeof(int) * (num_elements + 1));
  if (model_element->indices.env_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }
}

// Attach to a set of data
void attach_to_data(ModelElement *model_element, int num_elements,
                    int num_cells, int *int_data, double *float_data,
                    double *env_data, int *int_indices, int *float_indices,
                    int *env_indices) {
  // allocate the index arrays
  allocate_index_arrays(model_element, num_elements);

  // copy the indices
  for (int i = 0; i <= num_elements; ++i) {
    model_element->indices.int_data[i] = int_indices[i];
    model_element->indices.float_data[i] = float_indices[i];
    model_element->indices.env_data[i] = env_indices[i];
  }

  // set the number of grid cells
  model_element->n_cells = num_cells;

  // set data pointers
  model_element->ptrs.int_data = int_data;
  model_element->ptrs.float_data = float_data;
  model_element->ptrs.env_data = env_data;
}

// Copy a set of data
void copy_data(ModelElement from, ModelElement *to) {
  // allocate the index arrays
  allocate_index_arrays(to, from.indices.n_elements);

  // copy the indices
  for (int i = 0; i <= from.indices.n_elements; ++i) {
    to->indices.int_data[i] = from.indices.int_data[i];
    to->indices.float_data[i] = from.indices.float_data[i];
    to->indices.env_data[i] = from.indices.env_data[i];
  }

  // set the number of grid cells
  to->n_cells = from.n_cells;

  // allocate the data arrays
  to->ptrs.int_data =
      (int *)malloc(sizeof(int) * to->indices.int_data[to->indices.n_elements]);
  if (to->ptrs.int_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }
  to->ptrs.float_data = (double *)malloc(
      sizeof(double) * to->indices.float_data[to->indices.n_elements]);
  if (to->ptrs.float_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }
  to->ptrs.env_data = (double *)malloc(
      sizeof(double) * to->indices.env_data[to->indices.n_elements] *
      to->n_cells);
  if (to->ptrs.env_data == NULL) {
    printf("\n\nERROR allocating space for diff checker");
    exit(EXIT_FAILURE);
  }

  // copy the data
  for (int i = 0; i < to->indices.int_data[to->indices.n_elements]; ++i)
    to->ptrs.int_data[i] = from.ptrs.int_data[i];
  for (int i = 0; i < to->indices.float_data[to->indices.n_elements]; ++i)
    to->ptrs.float_data[i] = from.ptrs.float_data[i];
  for (int i = 0;
       i < to->indices.env_data[to->indices.n_elements] * to->n_cells; ++i)
    to->ptrs.env_data[i] = from.ptrs.env_data[i];
}

// Initialize the difference checker data
void diff_check_init(ModelData model_data) {
  DifferenceCheckerData *dd;
  dd = &diff_data[num_solvers++];

  printf("\nInitializing solver diff checker %d n_rxn = %d", num_solvers,
         model_data.n_rxn);

  // Point to the actual model data
  attach_to_data(&dd->current.reactions, model_data.n_rxn, model_data.n_cells,
                 model_data.rxn_int_data, model_data.rxn_float_data,
                 model_data.rxn_env_data, model_data.rxn_int_indices,
                 model_data.rxn_float_indices, model_data.rxn_env_idx);
  attach_to_data(
      &dd->current.aero_reps, model_data.n_aero_rep, model_data.n_cells,
      model_data.aero_rep_int_data, model_data.aero_rep_float_data,
      model_data.aero_rep_env_data, model_data.aero_rep_int_indices,
      model_data.aero_rep_float_indices, model_data.aero_rep_env_idx);
  attach_to_data(
      &dd->current.sub_models, model_data.n_sub_model, model_data.n_cells,
      model_data.sub_model_int_data, model_data.sub_model_float_data,
      model_data.sub_model_env_data, model_data.sub_model_int_indices,
      model_data.sub_model_float_indices, model_data.sub_model_env_idx);

  // Create a set of data to hold the previous values
  copy_data(dd->current.reactions, &dd->last_check.reactions);
  copy_data(dd->current.aero_reps, &dd->last_check.aero_reps);
  copy_data(dd->current.sub_models, &dd->last_check.sub_models);
}

// Compare and update two model elements
int compare_and_update(ModelElement current, ModelElement *last_check,
                       char *element_type, bool do_compare) {
  int num_diff = 0;

  if (current.indices.n_elements != last_check->indices.n_elements) {
    printf("\n  %s: n_element difference: current = %d, last check = %d",
           element_type, current.indices.n_elements,
           last_check->indices.n_elements);
    ++num_diff;
  }
  if (current.n_cells != last_check->n_cells) {
    printf("\n  %s: n_cells difference: current = %d, last check = %d",
           element_type, current.n_cells, last_check->n_cells);
    ++num_diff;
  }
  for (int e = 0; e < last_check->indices.n_elements; ++e) {
    // int data
    if (current.indices.int_data[e] != last_check->indices.int_data[e]) {
      printf(
          "\n  %s[%d]: start int data difference: current = %d, last check = "
          "%d",
          element_type, e, current.indices.int_data[e],
          last_check->indices.int_data[e]);
      ++num_diff;
    }
    if (current.indices.int_data[e + 1] !=
        last_check->indices.int_data[e + 1]) {
      printf(
          "\n  %s[%d]: end int data difference: current = %d, last check = %d",
          element_type, e, current.indices.int_data[e + 1],
          last_check->indices.int_data[e + 1]);
      ++num_diff;
    }
    int count = 0;
    for (int i = last_check->indices.int_data[e];
         i < last_check->indices.int_data[e + 1]; ++i) {
      if (current.ptrs.int_data[i] != last_check->ptrs.int_data[i] &&
          do_compare) {
        printf(
            "\n  %s[%d]: int datum %d (at index %d) difference: current = %d, "
            "last check = %d",
            element_type, e, count, i, current.ptrs.int_data[i],
            last_check->ptrs.int_data[i]);
        ++num_diff;
      }
      last_check->ptrs.int_data[i] = current.ptrs.int_data[i];
      ++count;
    }

    // float data
    if (current.indices.float_data[e] != last_check->indices.float_data[e]) {
      printf(
          "\n  %s[%d]: start float data difference: current = %d, last check = "
          "%d",
          element_type, e, current.indices.float_data[e],
          last_check->indices.float_data[e]);
      ++num_diff;
    }
    if (current.indices.float_data[e + 1] !=
        last_check->indices.float_data[e + 1]) {
      printf(
          "\n  %s[%d]: end float data difference: current = %d, last check = "
          "%d",
          element_type, e, current.indices.float_data[e + 1],
          last_check->indices.float_data[e + 1]);
      ++num_diff;
    }
    count = 0;
    for (int i = last_check->indices.float_data[e];
         i < last_check->indices.float_data[e + 1]; ++i) {
      if (current.ptrs.float_data[i] != last_check->ptrs.float_data[i] &&
          do_compare) {
        printf(
            "\n  %s[%d]: float datum %d (at index %d) difference: current = "
            "%lg, last check = %lg",
            element_type, e, count, i, current.ptrs.float_data[i],
            last_check->ptrs.float_data[i]);
        ++num_diff;
      }
      last_check->ptrs.float_data[i] = current.ptrs.float_data[i];
      ++count;
    }

    // env data
    if (current.indices.float_data[e] != last_check->indices.float_data[e]) {
      printf(
          "\n  %s[%d]: start float data difference: current = %d, last check = "
          "%d",
          element_type, e, current.indices.float_data[e],
          last_check->indices.float_data[e]);
      ++num_diff;
    }
    if (current.indices.float_data[e + 1] !=
        last_check->indices.float_data[e + 1]) {
      printf(
          "\n  %s[%d]: end float data difference: current = %d, last check = "
          "%d",
          element_type, e, current.indices.float_data[e + 1],
          last_check->indices.float_data[e + 1]);
      ++num_diff;
    }
    for (int i_cell = 0; i_cell < current.n_cells; ++i_cell) {
      double *curr_env_ptr =
          current.ptrs.env_data +
          i_cell * current.indices.env_data[current.indices.n_elements];
      double *last_env_ptr =
          last_check->ptrs.env_data +
          i_cell * last_check->indices.env_data[last_check->indices.n_elements];
      count = 0;
      for (int i = last_check->indices.env_data[e];
           i < last_check->indices.env_data[e + 1]; ++i) {
        if (curr_env_ptr[i] != last_env_ptr[i] && do_compare) {
          printf(
              "\n  %s[%d]: env datum %d (at index %d) for cell %d difference: "
              "current = %lg, last check = %lg",
              element_type, e, count, i, i_cell, curr_env_ptr[i],
              last_env_ptr[i]);
          ++num_diff;
        }
        last_env_ptr[i] = curr_env_ptr[i];
        ++count;
      }
    }

  }  // loop on elements

  return num_diff;
}

// Do a model data difference check
void diff_check(char *message) {
  int num_diff = 0;

  for (int i = 0; i < num_solvers; ++i) {
    DifferenceCheckerData *dd;
    dd = &diff_data[i];
    printf("\nchecking solver %d", i + 1);
    num_diff += compare_and_update(dd->current.reactions,
                                   &dd->last_check.reactions, "reaction", true);
    num_diff += compare_and_update(dd->current.aero_reps,
                                   &dd->last_check.aero_reps, "aero rep", true);
    num_diff += compare_and_update(
        dd->current.sub_models, &dd->last_check.sub_models, "sub model", true);
  }

  if (num_diff == 0) {
    printf("\nNO DIFFERENCES: %s\n", message);
  } else {
    printf("\n FOUND %i DIFFERENCES: %s\n", num_diff, message);
  }
  fflush(stdout);
}

// Do a model data update with no compare (except data dimensions)
void diff_check_update_only(char *message) {
  int num_diff = 0;

  for (int i = 0; i < num_solvers; ++i) {
    DifferenceCheckerData *dd;
    dd = &diff_data[i];
    num_diff += compare_and_update(
        dd->current.reactions, &dd->last_check.reactions, "reaction", false);
    num_diff += compare_and_update(
        dd->current.aero_reps, &dd->last_check.aero_reps, "aero rep", false);
    num_diff += compare_and_update(
        dd->current.sub_models, &dd->last_check.sub_models, "sub model", false);
  }

  if (num_diff == 0) {
    printf("\nDIFFERENCE CHECKER UPDATE: %s\n", message);
  } else {
    printf("\n DIFFERENCE CHECKER UPDATE WITH %i STRUCTURAL DIFFERENCES: %s\n",
           num_diff, message);
  }
  fflush(stdout);
}
