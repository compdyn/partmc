/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 */
/** \file
 * \brief common functions for c tests
 */
#ifndef PMC_TEST_COMMON_H
#define PMC_TEST_COMMON_H

#include "../src/camp_common.h"

// Tolerances
#define PMC_TEST_COMMON_ABS_TOL 1.0e-30
#define PMC_TEST_COMMON_REL_TOL 1.0e-10

// Assert function
#define ASSERT_MSG(y, z) pmc_assert(__func__, __LINE__, y, z);
#define ASSERT(y) pmc_assert(__func__, __LINE__, y, "Unknown error");
#define ASSERT_CLOSE_MSG(x, y, z) pmc_assert_close(__func__, __LINE__, x, y, z);
#define ASSERT_CLOSE(x, y) pmc_assert_close(__func__, __LINE__, x, y, "Unknown error");

// Assert function def
int pmc_assert(const char *func, const int line, bool eval,
               const char *message) {
  if (eval) {
    return 0;
  }
  printf("\n[ERROR] line %4d in %s(): %s", line, func, message);
  return 1;
}

// Assert close function def
int pmc_assert_close(const char *func, const int line, double val1, double val2, const char *message) {
  bool eval = val1 == val2 ? true :
              fabs(val1 - val2) <= PMC_TEST_COMMON_ABS_TOL ||
              2.0 * fabs(val1 - val2) / fabs(val1 + val2) <= PMC_TEST_COMMON_REL_TOL;
  return pmc_assert( func, line, eval, message);
}


#endif
