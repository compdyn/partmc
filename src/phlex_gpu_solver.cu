/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * GPU solver functions
 *
*/
/** \file
 * \brief GPU solver functions
*/
#include "aero_rep_solver.h"
#include "cuda_util.h"
#include "phlex_gpu_solver.h"
#include "rxn_gpu_solver.h"
#include "rxn_solver.h"
#include "sub_model_solver.h"

// TODO Should this be input data?
#define NUM_BLOCKS_  5
#define NUM_THREADS_ 100

// State variable types (Must match parameters defined in pmc_chem_spec_data module)
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

/** Create a new set of GPU solver data and add it to a ModelData object
 *
 * \param model_data Pointer to the model data to use to build GPU solver data
 */
void phlex_gpu_solver_new( ModelData * model_data )
{
  // Get the device data object
  ModelDeviceData * model_dev_data = 
          (ModelDeviceData*) (model_data->model_dev_data);

  // Set the nubmer of blocks and threads
  (*model_dev_data).num_blocks  = NUM_BLOCKS_;
  (*model_dev_data).num_threads = NUM_THREADS_;

  // Get a device pointer to the working state array
  HANDLE_ERROR( cudaHostRegister( model_data->state, 
                                  model_data->n_state_var * sizeof(int),
                                  0
                                ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                                  (void**) &(model_dev_data->dev_state), 
                                  model_data->state,
                                  0
                                ) );

  // Set up the working derivative array
  HANDLE_ERROR( cudaHostRegister( (void**) &(model_data->deriv),
                               model_data->deriv_size * sizeof(PMC_C_FLOAT),
                               0
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_deriv),
                               model_data->deriv,
                               0
                             ) );
  
  // Set up the working Jacobian data array
  HANDLE_ERROR( cudaHostRegister( (void**) &(model_data->jac),
                               model_data->jac_size * sizeof(PMC_C_FLOAT),
                               0
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_jac),
                               model_data->jac,
                               0
                             ) );

  // Initialize the reaction data
  rxn_gpu_solver_new( 
          ( RxnDeviceData* )( ( *model_dev_data ).rxn_dev_data ),
          model_data->rxn_data );

}

/** \brief Compute the time derivative f(t,y)
 *
 * \param t Current model time (s)
 * \param y Dependent variable array
 * \param deriv Time derivative vector f(t,y) to calculate
 * \param solver_data Pointer to the solver data
 * \return Status code
 */
int phlex_gpu_solver_f(realtype t, N_Vector y, N_Vector deriv, void *solver_data)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  ModelDeviceData *mdd = ( ModelDeviceData* ) ( md->model_dev_data );
  realtype time_step;

  // Update the state array with the current dependent variable values
  // Signal a recoverable error (positive return value) for negative 
  // concentrations.
  for (int i_spec=0, i_dep_var=0; i_spec<md->n_state_var; i_spec++) {
    if (md->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      if (NV_DATA_S(y)[i_dep_var] < 0.0) return 1;
      md->state[i_spec] = (PMC_C_FLOAT) (NV_DATA_S(y)[i_dep_var++]);
    }
  }

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Run pre-derivative calculations
  rxn_pre_calc(md);

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);
  
  // Calculate the time derivative f(t,y) for GPU rxns
  dim3 dimGrid( mdd->num_blocks );
  dim3 dimBlock( mdd->num_threads );
  rxn_gpu_calc_deriv<<< dimGrid, dimBlock >>>(*mdd, (PMC_C_FLOAT) time_step);

  // Copy working derivative array to solver derivative
  // TODO eliminate copying
  for (int i_spec=0; i_spec<NV_LENGTH_S(deriv); i_spec++)
          NV_DATA_S(deriv)[i_spec] = (realtype) mdd->host_deriv[i_spec];

  // Calculate the remaining time derivatives f(t,y)
  rxn_calc_deriv(md, deriv, (PMC_C_FLOAT) time_step);
  
  return (0);

}

void phlex_gpu_solver_free(ModelDeviceData * model_device_data)
{
  free( model_device_data );
}
