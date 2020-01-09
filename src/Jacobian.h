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

/* Jacobian for solver species */
typedef struct {
  unsigned int num_spec;   // Number of species
  unsigned int num_elem;   // Number of potentially non-zero Jacobian elements
  unsigned int *col_ptrs;  // Index of start/end of each column in data array
  unsigned int *row_ids;   // Row id of each Jacobian element in data array
  long double
      *production_partials;    // Data array for productions rate partial derivs
  long double *loss_partials;  // Data array for loss rate partial derivs
} Jacobian;

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

/** \brief Get an element id in the Jacobian data arrays
 *
 * \param jac Jacobian object
 * \param col_id Column index
 * \param row_id Row index
 * \return Index of Jacobian element in the data array
 */
unsigned int jacobian_get_element_id(Jacobian jac, unsigned int col_id,
                                     unsigned int row_id);

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

/** \brief Free memory associated with a Jacobian
 *
 * \param jac Jacobian object
 */
void jacobian_free(Jacobian jac);

#endif
