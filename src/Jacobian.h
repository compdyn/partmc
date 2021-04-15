/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header for the Jacobian structure and related functions
 *
 */
/** \file
 * \brief Header for the Jacobian structure and related functions
 */
#ifndef JACOBIAN_H_
#define JACOBIAN_H_

#include <math.h>
#include <stdlib.h>

// Flags for specifying production or loss elements
#define JACOBIAN_PRODUCTION 0
#define JACOBIAN_LOSS 1

/* Registered elements for a column in the Jacobian */
typedef struct {
  unsigned int array_size;  // Size of the array of flagged elements
  unsigned int
      number_of_elements;  // Number of registered elements in the column
  unsigned int
      *row_ids;  // Array of row ids for each registered element in the column
} JacobianColumnElements;

/* Jacobian for solver species */
typedef struct {
  unsigned int num_spec;   // Number of species
  unsigned int num_elem;   // Number of potentially non-zero Jacobian elements
  unsigned int *col_ptrs;  // Index of start/end of each column in data array
  unsigned int *row_ids;   // Row id of each Jacobian element in data array
  long double
      *production_partials;    // Data array for productions rate partial derivs
  long double *loss_partials;  // Data array for loss rate partial derivs
  JacobianColumnElements *elements;  // Jacobian elements flagged for inclusion
} Jacobian;

/** \brief Initialize the Jacobian
 *
 * Sets up a sparse matrix (num_spec x num_spec) containing zero elements.
 * Elements can be added using the \c jacobian_register_element function.
 *
 * \param jac Jacobian object
 * \param num_spec Number of species
 * \return Flag indicating whether the Jacobian was successfully initialized
 *         (0 = false; 1 = true)
 */
int jacobian_initialize_empty(Jacobian *jac, unsigned int num_spec);

/** \brief Initialize the Jacobian
 *
 * \param jac Pointer to the Jacobian object
 * \param num_spec Number of species
 * \param jac_struct Dense matrix of flags indicating whether an element is
 *                   (1) potentially non-zero or (0) not.
 * \return Flag indicating whether the derivative was successfully initialized
 *         (0 = false; 1 = true)
 */
int jacobian_initialize(Jacobian *jac, unsigned int num_spec,
                        unsigned int **jac_struct);

/** \brief Adds an element to the sparse matrix
 *
 * \param jac Jacobian object
 * \param dep_id Dependent species index
 * \param ind_id Independent species index
 */
void jacobian_register_element(Jacobian *jac, unsigned int dep_id,
                               unsigned int ind_id);

/** \brief Builds the sparse matrix with the registered elements
 *
 * \return 1 on success, 0 otherwise
 */
unsigned int jacobian_build_matrix(Jacobian *jac);

/** \brief Returns the number of elements in the Jacobian
 *
 * \param jac Jacobian object
 * \return Number of Jacobian elements
 */
unsigned int jacobian_number_of_elements(Jacobian jac);

/** \brief Returns the value of a column pointer
 *
 * \param jac Jacobian object
 * \param col_id Column index (0...number of columns)
 * \return Column pointer value
 */
unsigned int jacobian_column_pointer_value(Jacobian jac, unsigned int col_id);

/** \brief Returns the row for a given Jacobian element
 *
 * \param jac Jacobian object
 * \param elem_id Jacobian element index (0...number of elements-1)
 * \return Row index for given element
 */
unsigned int jacobian_row_index(Jacobian jac, unsigned int elem_id);

/** \brief Get an element id in the Jacobian data arrays
 *
 * If the element is not included in the sparse matrix, -1 is returned.
 *
 * \param jac Jacobian object
 * \param dep_id Dependent species index
 * \param ind_id Independent species index
 * \return Index of Jacobian element in the data array
 */
unsigned int jacobian_get_element_id(Jacobian jac, unsigned int dep_id,
                                     unsigned int ind_id);

/** \brief Reset the Jacobian
 *
 * \param jac Jacobian matrix
 */
void jacobian_reset(Jacobian jac);

/** \brief Output the Jacobian
 *
 * \param jac Jacobian object
 * \param dest_array Pointer to the array to save Jacobian data to
 */
void jacobian_output(Jacobian jac, double *dest_array);

/** \brief Add a contribution to the Jacobian
 *
 * \param jac Jacobian object
 * \param elem_id Index of the element to update in the data array
 * \param prod_or_loss Flag indicating whether to update the (0) production or
 *                          (1) loss elements
 * \param jac_contribution Value to add to the Jacobian element
 *                         (contributions to loss elements should be positive if
 *                         the contribution increases the loss)
 */
void jacobian_add_value(Jacobian jac, unsigned int elem_id,
                        unsigned int prod_or_loss,
                        long double jac_contribution);

/** \brief Prints the Jacobian structure
 *
 * \param jac Jacobian object
 */
void jacobian_print(Jacobian jac);

/** \brief Free memory associated with a JacobianColumnElements
 *
 * \param column Jacobian column elements
 */
void jacobian_column_elements_free(JacobianColumnElements *column);

/** \brief Free memory associated with a Jacobian
 *
 * \param jac Jacobian object
 */
void jacobian_free(Jacobian *jac);

#endif
