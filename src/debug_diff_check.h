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
#ifndef DEBUG_DIFF_CHECK_H_
#define DEBUG_DIFF_CHECK_H_

#include "camp_common.h"

// Initialize the difference checker data
void diff_check_init(ModelData model_data);

// Do a model data difference check
void diff_check(char* message);

// Update the last checked state without checking for differences
void diff_check_update_only(char* message);

#endif  // DEBUG_DIFF_CHECK_H_
