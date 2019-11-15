/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for reaction functions
 *
 * TODO Automatically generate rxn_solver.c and rxn_solver.h code
 * maybe using cmake?
 *
 */
/** \file
 * \brief Header file for reaction solver functions
 */
#ifndef RXNS2_H_
#define RXNS2_H_
#include "camp_common.h"

// arrhenius
#ifdef PMC_USE_SUNDIALS
void rxn2_arrhenius_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
                                      int *rxn_int_data, double *rxn_float_data,
                                      double *rxn_env_data, realtype time_step);

#endif


#endif
