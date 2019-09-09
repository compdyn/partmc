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

// Assert function
#define ASSERT_MSG(y, z) pmc_assert(__func__, __LINE__, y, z);
#define ASSERT(y) pmc_assert(__func__, __LINE__, y, "Unknown error");

// Assert function def
int pmc_assert(const char *func, const int line, bool eval,
               const char *message) {
  if (eval) {
    return 0;
  }
  printf("\n[ERROR] line %4d in %s(): %s", line, func, message);
  return 1;
}

#endif
