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

// State variable types (Must match parameters defined in pmc_chem_spec_data module)
#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

/** Create a new set of GPU solver data and add it to a SolverData object
 *
 * \param solver_data Solver data
 */
extern "C"
void phlex_gpu_solver_new( SolverData *solver_data )
{

  // Determine the current GPU memory usage
  size_t free_mem, total_mem;
  cudaMemGetInfo( &free_mem, &total_mem );
  printf("\nINFO (PartMC-8351420968): GPU total memory %zu, GPU free memory %zu\n",
    total_mem, free_mem);

  // Create a new SolverDeviceData object
  solver_data->solver_device_data = ( void * )
                                   malloc( sizeof( SolverDeviceData ) );
  if( solver_data->solver_device_data == NULL ) {
    printf("\n\nERROR allocating space for SolverDeviceData\n\n");
    exit( 1 );
  }

  // Get a pointer to the SolverDeviceData object
  SolverDeviceData *sdd = ( SolverDeviceData * ) 
                          solver_data->solver_device_data;
  
  // Save the number of states to solve
  sdd->n_states = solver_data->n_states;

  // Allocate the ModelDeviceData array
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &( sdd->host_model_dev_data ),
                               sdd->n_states * sizeof( ModelDeviceData ),
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer(
                               (void**) &( sdd->dev_model_dev_data ),
                               (void*) sdd->host_model_dev_data,
                               0
                             ) );

  // Start solving the first state on the first block
  sdd->deriv_threads     = 0;
  sdd->jac_threads       = 0;
  sdd->env_threads       = 0;
  int curr_deriv_block   = 0;
  int curr_jac_block     = 0;
  int curr_env_block     = 0;
  int curr_deriv_threads = 0;
  int curr_jac_threads   = 0;
  int curr_env_threads   = 0;
  int curr_deriv_id      = 0;
  int curr_jac_id        = 0;

  // Set up a ModelDeviceData object for each ModelData object
  for( int i_state = 0; i_state < sdd->n_states; i_state++ ) {

    // Get a pointer to the ModelData object
    ModelData * md = &( solver_data->model_data[ i_state ] );

    // Get a pointer to the device data object
    ModelDeviceData * mdd = &( sdd->host_model_dev_data[ i_state ] );

    // Set up the working state array
    HANDLE_ERROR( cudaHostAlloc( ( void** ) &( mdd->host_state ), 
                                 md->n_state_var * sizeof( PMC_C_FLOAT ),
                                 cudaHostAllocWriteCombined |
                                    cudaHostAllocMapped
                               ) );
    HANDLE_ERROR( cudaHostGetDevicePointer( 
                                 (void**) &(mdd->dev_state), 
                                 (void*) mdd->host_state,
                                 0
                               ) );

    // Set up the working environmental array
    HANDLE_ERROR( cudaHostAlloc( (void**) &(mdd->host_env), 
                                 md->n_env_var * sizeof(PMC_C_FLOAT),
                                 cudaHostAllocWriteCombined |
                                    cudaHostAllocMapped
                               ) );
    HANDLE_ERROR( cudaHostGetDevicePointer( 
                                 (void**) &(mdd->dev_env), 
                                 (void*) mdd->host_env,
                                 0
                               ) );

    // Set up the working derivative array
    HANDLE_ERROR( cudaHostAlloc( (void**) &(mdd->host_deriv),
                                 md->deriv_size * sizeof(PMC_SOLVER_C_FLOAT),
                                 cudaHostAllocMapped
                               ) );
    HANDLE_ERROR( cudaHostGetDevicePointer( 
                                 (void**) &(mdd->dev_deriv),
                                 (void*) mdd->host_deriv,
                                 0
                               ) );
    mdd->deriv_size = md->deriv_size;

    // Set up the working Jacobian data array
    HANDLE_ERROR( cudaHostAlloc( (void**) &(mdd->host_jac),
                                 ( md->jac_size > 0 ? md->jac_size : 1 ) 
                                      * sizeof(PMC_SOLVER_C_FLOAT),
                                 cudaHostAllocMapped
                               ) );
    HANDLE_ERROR( cudaHostGetDevicePointer( 
                                 (void**) &(mdd->dev_jac),
                                 (void*) mdd->host_jac,
                                 0
                               ) );
    mdd->jac_size = md->jac_size;

    // Initialize the reaction data
    rxn_gpu_solver_new( mdd, md->rxn_data );

    ////////////////////////////
    // Set ids for this state //
    ////////////////////////////
    
    int n_rxn = ( (RxnDeviceData*)(mdd->host_rxn_dev_data) )->n_rxn;
    
    // Advance to the next block once the deriv or jac array is too large 
    if( curr_deriv_id + mdd->deriv_size > MAX_SHARED_ARRAY_SIZE_ ||
        curr_deriv_threads + n_rxn > CUDA_MAX_THREADS ) {
      curr_deriv_block++;
      if( curr_deriv_threads > sdd->deriv_threads ) 
        sdd->deriv_threads = curr_deriv_threads;
      curr_deriv_threads = curr_deriv_id = 0;
    }
    if( curr_jac_id + mdd->jac_size > MAX_SHARED_ARRAY_SIZE_ ||
        curr_jac_threads + n_rxn > CUDA_MAX_THREADS ) {
      curr_jac_block++;
      if( curr_jac_threads > sdd->jac_threads ) 
        sdd->jac_threads = curr_jac_threads;
      curr_jac_threads = curr_jac_id = 0;
    }
    if( curr_env_threads + n_rxn > CUDA_MAX_THREADS ) {
      curr_env_block++;
      if( curr_env_threads > sdd->env_threads )
        sdd->env_threads = curr_env_threads;
      curr_env_threads = 0;
    }

    // Set the current state ids
    mdd->deriv_block    = curr_deriv_block;
    mdd->jac_block      = curr_jac_block;
    mdd->env_block      = curr_env_block;
    mdd->deriv_start_id = curr_deriv_id;
    mdd->jac_start_id   = curr_jac_id;

    // Advance the number of threads and array elements
    curr_deriv_threads += n_rxn;
    curr_jac_threads   += n_rxn;
    curr_env_threads   += n_rxn;
    curr_deriv_id      += mdd->deriv_size;
    curr_jac_id        += mdd->jac_size;
  }

  if( curr_deriv_threads > sdd->deriv_threads )
    sdd->deriv_threads = curr_deriv_threads;
  if( curr_jac_threads > sdd->jac_threads )
    sdd->jac_threads   = curr_jac_threads;
  if( curr_env_threads > sdd->env_threads )
    sdd->env_threads   = curr_env_threads;
  sdd->deriv_blocks    = curr_deriv_block + 1;
  sdd->jac_blocks      = curr_jac_block + 1;
  sdd->env_blocks      = curr_env_block + 1;
}

/** \brief Update the environmental state
 *
 * \param SolverData Solver data
 */
extern "C"
void phlex_gpu_solver_update_env_state( SolverData *solver_data )
{
  SolverDeviceData * sdd = ( SolverDeviceData * )
                           solver_data->solver_device_data;

  // Update the environmental state for GPU functions
  for( int i_state = 0; i_state < solver_data->n_states; i_state++ ) {
    ModelData * md = &( solver_data->model_data[ i_state ] );
    ModelDeviceData * mdd = &( sdd->host_model_dev_data[ i_state ] );
    for( int i_var = 0; i_var < md->n_env_var; i_var++ ) {
      mdd->host_env[ i_var ] = md->env[ i_var ];
    }
  }

  // Update the environmental state for reactions with GPU functions
  dim3 dimGrid( sdd->env_blocks );
  dim3 dimBlock( sdd->env_threads );
  rxn_gpu_update_env_state<<< dimGrid, dimBlock >>>( *sdd );
  cudaDeviceSynchronize();

  // Update the remaining reactions
  for( int i_state = 0; i_state < solver_data->n_states; i_state++ ) {
    ModelData * md = &( solver_data->model_data[ i_state ] );
    rxn_update_env_state( *md );
  }

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
int phlex_gpu_solver_f( realtype t, N_Vector y, N_Vector deriv,
          void *solver_data )
{
  SolverData *sd = (SolverData*) solver_data;
  SolverDeviceData *sdd = ( SolverDeviceData* ) ( sd->solver_device_data );
  realtype time_step;

  // Loop through the states to solve
  for( int i_dep_var = 0, i_state = 0; i_state < sd->n_states; i_state++ ) {
    
    // Get pointers to the ModelData and ModelDeviceData for this state
    ModelData *md = &( sd->model_data[ i_state ] );
    ModelDeviceData *mdd = &( sdd->host_model_dev_data[ i_state ] );

    // Update the state array with the current dependent variable values
    // Signal a recoverable error (positive return value) for negative 
    // concentrations.
    for( int j_spec = 0, j_dep_var = 0; j_spec < md->n_state_var; j_spec++ ) {
      if( md->var_type[ j_spec ] == CHEM_SPEC_VARIABLE ) {
        if( NV_DATA_S( y )[ i_dep_var ] < 0.0 ) return 1;
        md->state[ j_spec ] = 
                ( PMC_C_FLOAT ) ( NV_DATA_S( y )[ i_dep_var++ ] );
        mdd->host_deriv[ j_dep_var++ ] = 0.0;     
      }
    }

    // Update the aerosol representations
    aero_rep_update_state( *md );

    // Run the sub models
    sub_model_calculate( *md );

    // Run pre-derivative calculations
    rxn_pre_calc( *md );

    // Update the state array for use on the GPUs
    for( int j_spec = 0; j_spec < md->n_state_var; j_spec++ )
        mdd->host_state[ j_spec ] = md->state[ j_spec ];

  }

  // Get the current integrator time step (s)
  CVodeGetCurrentStep( sd->cvode_mem, &time_step );
  
  // Calculate the time derivative f(t,y) for GPU rxns
  dim3 dimGrid( sdd->deriv_blocks );
  dim3 dimBlock( sdd->deriv_threads );
  rxn_gpu_calc_deriv<<< dimGrid, dimBlock >>>( *sdd, (PMC_C_FLOAT) time_step );
  cudaDeviceSynchronize();

  // Loop through the states to solve
  for( int i_dep_var = 0, i_state = 0; i_state < sd->n_states; i_state++ ) {
  
    // Get pointers to the ModelData and ModelDeviceData for this state
    ModelData *md = &( sd->model_data[ i_state ] );
    ModelDeviceData *mdd = &( sdd->host_model_dev_data[ i_state ] );

    // Calculate the remaining time derivatives f(t,y)
    rxn_calc_deriv( *md, mdd->host_deriv, (PMC_C_FLOAT) time_step );
  
    // Copy working derivative array to solver derivative
    for( int i_spec = 0; i_spec < md->deriv_size; i_spec++ )
          NV_DATA_S( deriv )[ i_dep_var++ ] = ( realtype ) mdd->host_deriv[ i_spec ];
  
  }

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
int phlex_gpu_solver_Jac( realtype t, N_Vector y, N_Vector deriv, SUNMatrix J,
        void *solver_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  SolverData *sd = (SolverData*) solver_data;
  SolverDeviceData *sdd = ( SolverDeviceData* ) ( sd->solver_device_data );
  realtype time_step;

  // TODO Figure out how to keep the Jacobian from being redimensioned
  // Reset the Jacobian dimensions
  if ( SM_NNZ_S( J ) < SM_NNZ_S( sd->J_init ) ) {
    SM_INDEXVALS_S( J ) = ( sunindextype* ) realloc( SM_INDEXVALS_S( J ),
              SM_NNZ_S( sd->J_init ) * sizeof( sunindextype ) );
    if( SM_INDEXVALS_S( J ) == NULL ) {
      printf( "\n\nERROR allocating space for sparse matrix index values\n\n" );
      exit( 1 );
    }
    SM_DATA_S( J ) = ( realtype* ) realloc( SM_DATA_S( J ),
              SM_NNZ_S( sd->J_init ) * sizeof( realtype ) );
    if ( SM_DATA_S( J ) == NULL ) {
      printf( "\n\nERROR allocating space for sparse matrix data\n\n" );
      exit( 1 );
    }
  }
  SM_NNZ_S( J ) = SM_NNZ_S( sd->J_init );
  for( int i = 0; i < SM_NNZ_S( J ); i++ ) {
    ( SM_INDEXVALS_S( J ) )[ i ] = ( SM_INDEXVALS_S( sd->J_init ) )[ i ];
  }
  for( int i = 0; i <= SM_NP_S( J ); i++ ) {
    ( SM_INDEXPTRS_S( J ) )[ i ] = ( SM_INDEXPTRS_S( sd->J_init ) )[ i ];
  } 

  // Loop through the states to solve
  for( int i_dep_var = 0, i_state = 0; i_state < sd->n_states; i_state++ ) {
    
    // Get pointers to the ModelData and ModelDeviceData for this state
    ModelData *md = &( sd->model_data[ i_state ] );
    ModelDeviceData *mdd = &( sdd->host_model_dev_data[ i_state ] );

    // Reset the Jacobian data for this state
    for( int i_elem = 0; i_elem < md->jac_size; i_elem++ )
      mdd->host_jac[ i_elem ] = ( PMC_SOLVER_C_FLOAT ) ZERO;
  
    // Update the state array with the current dependent variable values
    // Signal a recoverable error (positive return value) for negative 
    // concentrations.
    for( int j_spec = 0; j_spec < md->n_state_var; j_spec++ ) {
      if( md->var_type[ j_spec ] == CHEM_SPEC_VARIABLE ) {
        if( NV_DATA_S( y )[ i_dep_var ] < 0.0 ) return 1;
        md->state[ j_spec ] = 
                ( PMC_C_FLOAT ) ( NV_DATA_S( y )[ i_dep_var++ ] );
      }
    }

    // Update the aerosol representations
    aero_rep_update_state( *md );

    // Run the sub models
    sub_model_calculate( *md );

    // Run pre-derivative calculations
    rxn_pre_calc( *md );

    // Update the state array for use on the GPUs
    for( int j_spec = 0; j_spec < md->n_state_var; j_spec++ )
        mdd->host_state[ j_spec ] = md->state[ j_spec ];

  }

  // Get the current integrator time step (s)
  CVodeGetCurrentStep(sd->cvode_mem, &time_step);
  
  // Calculate the Jacobian for GPU rxns
  dim3 dimGrid( sdd->jac_blocks );
  dim3 dimBlock( sdd->jac_threads );
  rxn_gpu_calc_jac<<< dimGrid, dimBlock >>>( *sdd, (PMC_C_FLOAT) time_step );
  cudaDeviceSynchronize();

  // Loop through the states to solve
  for( int i_jac_elem = 0, i_state = 0; i_state < sd->n_states; i_state++ ) {
  
    // Get pointers to the ModelData and ModelDeviceData for this state
    ModelData *md = &( sd->model_data[ i_state ] );
    ModelDeviceData *mdd = &( sdd->host_model_dev_data[ i_state ] );

    // Calculate the Jacobian for the remaining rxns
    rxn_calc_jac( *md, mdd->host_jac, time_step );

    // Copy the working Jacobian back into the solver Jacobian
    for( int i = 0; i < md->jac_size; i++ )
      SM_DATA_S( J )[ i_jac_elem++ ] = ( realtype )( mdd->host_jac[ i ] );

  }

  return (0);

}

/** \brief Free GPU solver device data memory from a void pointer
  */
extern "C"
void phlex_gpu_solver_solver_device_data_free_vp( void * solver_device_data )
{
  SolverDeviceData *sd = ( SolverDeviceData * ) solver_device_data;
  phlex_gpu_solver_solver_device_data_free( *sd );
  free( sd );
}

/** \brief Free GPU solver device data memory
  */
extern "C"
void phlex_gpu_solver_solver_device_data_free( 
          SolverDeviceData solver_device_data )
{
  for( int i_state = 0; i_state < solver_device_data.n_states; i_state++ ) {
    phlex_gpu_solver_model_device_data_free( 
              solver_device_data.host_model_dev_data[ i_state ] );
  } 
  HANDLE_ERROR( cudaFreeHost( solver_device_data.host_model_dev_data ) );
}

/** \brief Free GPU model device data memory
  */
extern "C"
void phlex_gpu_solver_model_device_data_free(
          ModelDeviceData model_device_data )
{
  rxn_gpu_solver_free( model_device_data.host_rxn_dev_data );
  HANDLE_ERROR( cudaFreeHost( model_device_data.host_state ) );
  HANDLE_ERROR( cudaFreeHost( model_device_data.host_deriv ) );
  HANDLE_ERROR( cudaFreeHost( model_device_data.host_jac ) );
}
