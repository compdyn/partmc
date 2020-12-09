/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Debug and stats functions
 *
 */
// todo fix this...move all to this folder and compile correctly...
// todo: i think is not necessary include all time stdio.h.. reorder includes to left only camp_common with the common includes

// todo gprof to csv: https://stackoverflow.com/questions/28872400/convert-simple-ascii-table-to-csv
//Nvidia remote profiling: http://docs.nvidia.com/cuda/nsight-eclipse-edition-getting-started-guide/index.html#remote-development


#include "camp_debug_2.h"
#include <stdio.h>
#include <stdlib.h>

/** \brief Print derivative array
 *
 * \param deriv Derivative array
 */
static void print_derivative_2(N_Vector deriv) {
  // printf(" deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
    printf(" deriv: % -le", NV_DATA_S(deriv)[i]);
    printf(" index: %d \n", i);
  }
}