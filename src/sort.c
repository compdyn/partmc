/* Copyright (C) 2011 Matthew West
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 */

/** \file
 * \brief Wrapper routines for C qsort.
 */

#include <stdlib.h>

/** \brief Helper function for integer_sort_c()
 */
int pair_compare(const void *a, const void *b)
{
        int a_val = *((int*)a);
        int b_val = *((int*)b);

        if (a_val < b_val) {
                return -1;
        } else if (a_val > b_val) {
                return 1;
        } else {
                return 0;
        }
}

/** \brief Sort the given data array and return the permutation.
 *
 * On return the \c data array is sorted and the \c perm array
 * contains the permutation, so that <tt>new_data[i] =
 * data[perm[i]]</tt>, where \c data is the original data and \c
 * new_data is the sorted data.
 *
 * \param n The length of \c data and \c perm.
 * \param data The data array (sorted on return).
 * \param perm The permutation on return: <tt>new_data[i] =
 * data[perm[i]]</tt>.
 */
int integer_sort_c(int n, int *data, int *perm)
{
        int *data_perm;
        int i;

        data_perm = (int*)malloc(sizeof(int) * 2 * n);
        for (i = 0; i < n; i++) {
                data_perm[2 * i] = data[i];
                data_perm[2 * i + 1] = i + 1;
        }
        qsort(data_perm, n, 2 * sizeof(int), pair_compare);
        for (i = 0; i < n; i++) {
                data[i] = data_perm[2 * i];
                perm[i] = data_perm[2 * i + 1];
        }
        free(data_perm);
}
