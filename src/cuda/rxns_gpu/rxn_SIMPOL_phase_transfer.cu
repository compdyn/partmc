/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Phase Transfer reaction solver functions
 *
*/
/** \file
 * \brief Phase Transfer reaction solver functions
*/
extern "C"{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../aeros_gpu/aero_rep_solver_gpu.h"
#include "../aeros_gpu/sub_model_solver_gpu.h"

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]
 //todo fix
// Universal gas constant (J/mol/K)
#define UNIV_GAS_CONST_ 8.314472
// Small number for ignoring low concentrations
#define VERY_SMALL_NUMBER_ 1.0e-30

// Jacobian set indices
#define JAC_GAS 0
#define JAC_AERO 1

#define DELTA_H_ float_data[0*n_rxn]
#define DELTA_S_ float_data[1*n_rxn]
#define DIFF_COEFF_ float_data[2*n_rxn]
#define PRE_C_AVG_ float_data[3*n_rxn]
#define B1_ float_data[4*n_rxn]
#define B2_ float_data[5*n_rxn]
#define B3_ float_data[6*n_rxn]
#define B4_ float_data[7*n_rxn]
#define CONV_ float_data[8*n_rxn]
#define MW_ float_data[9*n_rxn]
#define SMALL_NUMBER_ float_data[10*n_rxn]
#define NUM_AERO_PHASE_ int_data[0*n_rxn]
#define GAS_SPEC_ (int_data[1*n_rxn]-1)
#define C_AVG_ALPHA_ rxn_env_data[0]
#define EQUIL_CONST_ rxn_env_data[1]
#define UGM3_TO_PPM_ rxn_env_data[2]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 11
#define NUM_ENV_PARAM_ 3
#define AERO_SPEC_(x) (int_data[(NUM_INT_PROP_ + x)*n_rxn]-1)
#define AERO_ACT_ID_(x) (int_data[(NUM_INT_PROP_ + NUM_AERO_PHASE_ + x)*n_rxn])
#define AERO_PHASE_ID_(x) (int_data[(NUM_INT_PROP_ + 2*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define AERO_REP_ID_(x) (int_data[(NUM_INT_PROP_ + 3*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define DERIV_ID_(x) (int_data[(NUM_INT_PROP_ + 4*(NUM_AERO_PHASE_) + x)*n_rxn])
#define GAS_ACT_JAC_ID_(x) int_data[(NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x)*n_rxn]
#define AERO_ACT_JAC_ID_(x) int_data[(NUM_INT_PROP_ + 1 + 6*(NUM_AERO_PHASE_) + x)*n_rxn]
#define JAC_ID_(x) (int_data[(NUM_INT_PROP_ + 1 + 5*(NUM_AERO_PHASE_) + x)*n_rxn])
#define PHASE_INT_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 10*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define PHASE_FLOAT_LOC_(x) (int_data[(NUM_INT_PROP_ + 2 + 11*(NUM_AERO_PHASE_) + x)*n_rxn]-1)
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[PHASE_INT_LOC_(x)*n_rxn])
#define PHASE_JAC_ID_(x,s,e) int_data[(PHASE_INT_LOC_(x)+1+s*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define EFF_RAD_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+e]
#define NUM_CONC_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define MASS_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+2*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]
#define MW_JAC_ELEM_(x,e) float_data[(PHASE_FLOAT_LOC_(x)+3*NUM_AERO_PHASE_JAC_ELEM_(x)+e)*n_rxn]

}