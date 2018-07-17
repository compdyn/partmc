/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for GPU solver functions for Arrhenius reactions
 *
*/
/** \file
 * \brief Header file for GPU solver functions for Arrhenius reactions
*/
#ifndef RXN_GPU_ARRHENIUS_H_
#define RXN_GPU_ARRHENIUS_H_
#include "../phlex_gpu_solver.h"

void * rxn_gpu_arrhenius_calc_deriv_contrib( void *rxn_data, 
         ModelDeviceData mdd, realtype *deriv );

#endif
