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
#include "cuda_util.h"
extern "C" {
#include "aero_rep_solver.h"
#include "phlex_gpu_solver.h"
#include "rxn_gpu_solver.h"
#include "rxn_solver.h"
#include "sub_model_solver.h"
}

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
extern "C"
void phlex_gpu_solver_new( ModelData * model_data )
{
  // Allocate a new ModelDeviceData object
  model_data->model_dev_data = (ModelDeviceData*) 
                                  malloc( sizeof(ModelDeviceData) );

  // Get the device data object
  ModelDeviceData * model_dev_data = 
          (ModelDeviceData*) (model_data->model_dev_data);

  // Set the nubmer of blocks and threads
  (*model_dev_data).num_blocks  = NUM_BLOCKS_;
  (*model_dev_data).num_threads = NUM_THREADS_;

  // Set up the working state array
  HANDLE_ERROR( cudaHostAlloc( (void**) &(model_dev_data->host_state), 
                               model_data->n_state_var * sizeof(PMC_C_FLOAT),
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_state), 
                               (void*) model_dev_data->host_state,
                               0
                             ) );

  // Set up the working environmental array
  HANDLE_ERROR( cudaHostAlloc( (void**) &(model_dev_data->host_env), 
                               model_data->n_env_var * sizeof(PMC_C_FLOAT),
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_env), 
                               (void*) model_dev_data->host_env,
                               0
                             ) );

  // Set up the working derivative array
  HANDLE_ERROR( cudaHostAlloc( (void**) &(model_dev_data->host_deriv),
                               model_data->deriv_size * sizeof(PMC_SOLVER_C_FLOAT),
                               cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_deriv),
                               (void*) model_dev_data->host_deriv,
                               0
                             ) );
  (*model_dev_data).deriv_size = model_data->deriv_size;

  // Set up the working Jacobian data array
  HANDLE_ERROR( cudaHostAlloc( (void**) &(model_dev_data->host_jac),
                               ( model_data->jac_size > 0 ? model_data->jac_size : 1 ) 
                                    * sizeof(PMC_SOLVER_C_FLOAT),
                               cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( 
                               (void**) &(model_dev_data->dev_jac),
                               (void*) model_dev_data->host_jac,
                               0
                             ) );
  (*model_dev_data).jac_size = model_data->jac_size;

  // Initialize the reaction data
  rxn_gpu_solver_new( model_dev_data, model_data->rxn_data );

}

/** \brief Update the environmental state
 *
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
extern "C"
void phlex_gpu_solver_update_env_state(ModelData *model_data, PMC_C_FLOAT *env)
{
  ModelDeviceData *mdd = ( ModelDeviceData* ) ( model_data->model_dev_data );

  // Update the environmental state for GPU functions
  for( int i_var = 0; i_var < model_data->n_env_var; i_var++ )
    mdd->host_env[ i_var ] = env[ i_var ];

  // Update the environmental state for reactions with GPU functions
  dim3 dimGrid( mdd->num_blocks );
  dim3 dimBlock( mdd->num_threads );
  rxn_gpu_update_env_state<<< dimGrid, dimBlock >>>( *mdd );
  cudaDeviceSynchronize();

  // Update the remaining reactions
  rxn_update_env_state( model_data, env );

}

/** \brief Compute the time derivative f(t,y)
 *
 * \param t Current model time (s)
 * \param y Dependent variable array
 * \param deriv Time derivative vector f(t,y) to calculate
 * \param solver_data Pointer to the solver data
 * \return Status code
 */
extern "C"
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
      mdd->host_state[i_spec] = md->state[i_spec] = (PMC_C_FLOAT) (NV_DATA_S(y)[i_dep_var]);
      mdd->host_deriv[i_dep_var++] = 0.0;     
    } else {
      mdd->host_state[i_spec] = md->state[i_spec];
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
  cudaDeviceSynchronize();

  // Calculate the remaining time derivatives f(t,y)
  rxn_calc_deriv(md, mdd->host_deriv, (PMC_C_FLOAT) time_step);
  
  // Copy working derivative array to solver derivative
  for (int i_spec=0; i_spec<NV_LENGTH_S(deriv); i_spec++)
          NV_DATA_S(deriv)[i_spec] = (realtype) mdd->host_deriv[i_spec];

  return (0);

}

/** \brief Compute the Jacobian
 *
 * \param t Current model time (s)
 * \param y Dependent variable array
 * \param deriv Time derivative vector f(t,y)
 * \param J Jacobian to calculate
 * \param solver_data Pointer to the solver data
 * \param tmp1 Unused vector
 * \param tmp2 Unused vector
 * \param tmp3 Unused vector
 * \return Status code
 */
int phlex_gpu_solver_Jac(realtype t, N_Vector y, N_Vector deriv, SUNMatrix J,
        void *solver_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SolverData *sd = (SolverData*) solver_data;
  ModelData *md = &(sd->model_data);
  ModelDeviceData *mdd = ( ModelDeviceData* ) ( md->model_dev_data );
  realtype time_step;
  PMC_SOLVER_C_FLOAT *J_data;

  // Update the state array with the current dependent variable values
  for (int i_spec=0, i_dep_var=0; i_spec<md->n_state_var; i_spec++) {
    if (md->var_type[i_spec]==CHEM_SPEC_VARIABLE) {
      if (NV_DATA_S(y)[i_dep_var] < 0.0) return 1;
      mdd->host_state[i_spec] = md->state[i_spec] = (PMC_C_FLOAT) (NV_DATA_S(y)[i_dep_var++]);
    } else {
      mdd->host_state[i_spec] = md->state[i_spec];
    }
  }

  // Get a pointer to the working Jacobian data array
  J_data = mdd->host_jac;
  
  // TODO Figure out how to keep the Jacobian from being redimensioned
  // Reset the Jacobian dimensions
  if (SM_NNZ_S(J)<SM_NNZ_S(md->J_init)) {
    SM_INDEXVALS_S(J) = (sunindextype*) realloc(SM_INDEXVALS_S(J),
              SM_NNZ_S(md->J_init)*sizeof(sunindextype));
    if (SM_INDEXVALS_S(J)==NULL) {
      printf("\n\nERROR allocating space for sparse matrix index values\n\n");
      exit(1);
    }
    SM_DATA_S(J) = (realtype*) realloc(SM_DATA_S(J),
              SM_NNZ_S(md->J_init)*sizeof(realtype));
    if (SM_DATA_S(J)==NULL) {
      printf("\n\nERROR allocating space for sparse matrix data\n\n");
      exit(1);
    }
  }
  SM_NNZ_S(J) = SM_NNZ_S(md->J_init);
  for (int i=0; i<SM_NNZ_S(J); i++) {
    J_data[i] = (PMC_SOLVER_C_FLOAT) 0.0;
    (SM_INDEXVALS_S(J))[i] = (SM_INDEXVALS_S(md->J_init))[i];
  }
  for (int i=0; i<=SM_NP_S(J); i++) {
    (SM_INDEXPTRS_S(J))[i] = (SM_INDEXPTRS_S(md->J_init))[i];
  } 

  // Update the aerosol representations
  aero_rep_update_state(md);

  // Run the sub models
  sub_model_calculate(md);

  // Run pre-Jacobian calculations
  rxn_pre_calc(md);

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);
  
  // Calculate the Jacobian for GPU rxns
  dim3 dimGrid( mdd->num_blocks );
  dim3 dimBlock( mdd->num_threads );
  rxn_gpu_calc_jac<<< dimGrid, dimBlock >>>(*mdd, (PMC_C_FLOAT) time_step);
  cudaDeviceSynchronize();

  // Calculate the Jacobian for the remaining rxns
  rxn_calc_jac(md, J_data, time_step);

  // Copy the working Jacobian back into the solver Jacobian
  for (int i=0; i<SM_NNZ_S(J); i++)
    SM_DATA_S(J)[i] = (realtype) (J_data[i]);

  return (0);

}

/** \brief Free GPU solver memory
  */
extern "C"
void phlex_gpu_solver_free(ModelDeviceData * model_device_data)
{
  rxn_gpu_solver_free( model_device_data->host_rxn_dev_data );
  HANDLE_ERROR( cudaFreeHost( model_device_data->host_state ) );
  HANDLE_ERROR( cudaFreeHost( model_device_data->host_deriv ) );
  HANDLE_ERROR( cudaFreeHost( model_device_data->host_jac ) );
  free( model_device_data );
}
