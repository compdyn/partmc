/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Single particle aerosol representation functions
 *
 */
/** \file
 * \brief Single particle aerosol representation functions
 */
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_phase_solver_gpu.h"
#include "../aeros_gpu/aero_reps_gpu.h"
#include "../camp_gpu_solver.h"

  /*

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define UPDATE_RADIUS 0
#define UPDATE_NUMBER 1

#define NUM_PHASE_ int_data[0]
#define AERO_REP_ID_ int_data[1]
#define RADIUS_ aero_rep_env_data[0]
#define NUMBER_CONC_ aero_rep_env_data[1]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 0
#define NUM_ENV_PARAM_ 2
#define PHASE_STATE_ID_(x) (int_data[NUM_INT_PROP_+x]-1)
#define PHASE_MODEL_DATA_ID_(x) (int_data[NUM_INT_PROP_+NUM_PHASE_+x]-1)
#define PHASE_NUM_JAC_ELEM_(x) int_data[NUM_INT_PROP_+2*NUM_PHASE_+x]
#define PHASE_MASS_(x) (float_data[NUM_FLOAT_PROP_+x])
#define PHASE_AVG_MW_(x) (float_data[NUM_FLOAT_PROP_+NUM_PHASE_+x])

*/


}
