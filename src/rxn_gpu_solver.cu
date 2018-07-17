/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Reaction-specific functions for use by the solver
 *
*/
/** \file
 * \brief Reaction solver functions
*/
#include "cuda_util.h"
#include "phlex_solver.h"
#include "rxn_gpu_solver.h"
#include "rxn_solver.h"
#include "rxns/rxn_gpu_arrhenius.h"

// TODO Figure out how to use allocatable shared derivative array
#define MAX_DERIV_SIZE_ 500

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5
#define RXN_HL_PHASE_TRANSFER 6
#define RXN_AQUEOUS_EQUILIBRIUM 7
#define RXN_ZSR_AEROSOL_WATER 8
#define RXN_PDFITE_ACTIVITY 9
#define RXN_SIMPOL_PHASE_TRANSFER 10
#define RXN_CONDENSED_PHASE_ARRHENIUS 11

/** \brief Assemble a set of indices for each reaction to solve with GPUs
 *
 * \param rxn_dev_data Pointer to the device reaction data to set
 * \param host_rxn_data Pointer to the host reaction data
 */
void rxn_gpu_solver_new( RxnDeviceData * rxn_dev_data, void * host_rxn_data )
{
  // Get the number of reactions
  int *rxn_data = (int*) (host_rxn_data);
  int n_rxn = *(rxn_data++);

  // Count the number and size of reactions with GPU solver functions
  int n_gpu_rxn = 0;
  size_t size_gpu_rxn = sizeof(int);
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++ ) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);
    
    // Save the starting position of the current reaction's data
    char * first_datum = (char*) rxn_data;

      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          n_gpu_rxn++;
          rxn_data = (int*) rxn_arrhenius_skip(
                    (void*) rxn_data);
          size_gpu_rxn += sizeof(char) * ((char*)rxn_data - first_datum);
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                    (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }
  }

  // Allocate space for the reaction parameters
  HANDLE_ERROR( cudaHostAlloc( (void**) &((*rxn_dev_data).host_rxn_data),
                               size_gpu_rxn,
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( &((*rxn_dev_data).dev_rxn_data),
                                          (*rxn_dev_data).host_rxn_data,
                                          0 
                                        ) );

  // Allocate space for the starting index of each reaction
  HANDLE_ERROR( cudaHostAlloc( (void**) &((*rxn_dev_data).host_rxn_data_start),
                               n_gpu_rxn * sizeof(int),
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( &((*rxn_dev_data).dev_rxn_data_start),
                                          (*rxn_dev_data).host_rxn_data_start,
                                          0
                                        ) );
  
  int * host_rxn_data = (int*) (*rxn_dev_data).host_rxn_data;
  char * host_rxn_data_init = (char*) host_rxn_data;
  int * host_rxn_data_start = (*rxn_dev_data).host_rxn_data_start;
 
  // Set the number of reactions
  *(host_rxn_data++) = n_gpu_rxn;

  // Loop through the reactions, copying rxn data to the rxn data block
  for (int i_rxn = 0, i_gpu_rxn = 0; i_rxn < n_rxn; i_rxn++ ) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    char * host_datum;
    char * rxn_datum;

      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          host_rxn_data_start[ i_gpu_rxn++ ] = 
                  ((char*) host_rxn_data) - host_rxn_data_init;
          *(host_rxn_data++) = rxn_type;
          host_datum = (char*) host_rxn_data;
          rxn_datum = (char*) rxn_data;
          rxn_data = (int*) rxn_arrhenius_skip(
                    (void*) rxn_data);
          for ( ; datum < (char*) rxn_data; ) *(host_datum++) = *(rxn_datum++);
          host_rxn_data = (int*) host_datum;
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                    (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }
  }
}

/** \brief Calculate the time derivative for reactions with GPU solver functions
 * 
 * \param device_data Device data needed for solving
 * \param time_step Current solver time step (s)
 */
__global__ void rxn_gpu_calc_deriv( ModelDeviceData mdd, PMC_C_FLOAT time_step )
{

  // Get the reaction data
  RxnDeviceData * rd = ( RxnDeviceData* ) ( mdd.rxn_dev_data );

  // Get a unique device index
  int dev_id = blockIdx.x * blockDim.x + threadIdx.x;
  int dev_total = gridDim.x * blockDim.x;

  // Set up a shared derivative array
  __shared__ float shr_dev_deriv[ MAX_DERIV_SIZE_ ];

  // Initialize the derivative array
  for( int i_spec = threadIdx.x; i_spec < mdd.deriv_size; i_spec += blockDim.x )
    shr_dev_deriv[ i_spec ] = 0.0;
  __syncthreads();

  // Get the number of reactions
  int * rxn_data = (int*) (*rd).dev_rxn_data;
  int n_rxn = *(rxn_data);

  // Figure out what reactions to solve on this thread
  int rxn_start = dev_id * n_rxn / dev_total  + 
                  ( dev_id > n_rxn % dev_total ? n_rxn % dev_total : dev_id );
  int rxns_to_solve = n_rxn / dev_total +
                  ( dev_id < n_rxn % dev_total ? 1 : 0 );

  // Advance the rxn data pointer to the first reaction's data
  char first_rxn = ( (char*) rxn_data ) + 
                   (*rd).dev_rxn_data_start[ rxn_start ];
  rxn_data = (int*) first_rxn;
  for( int i_rxn = 0; i_rxn < rxns_to_solve; i_rxn++ ) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Add derivative contribution from appropriate reaction type
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_arrhenius_calc_deriv_contrib( 
                            rxn_data, mdd, shr_dev_deriv );
        break;
      default :
        printf("\nPartMC Internal Error: invalid rxn type in GPU solver.\n");
        exit 1;
    }
  }
  __syncthreads();

  // Add derivative contributions from this block to the primary deriv array
  for( int i_spec = threadIdx.x; i_spec < mdd.deriv_size; i_spec += blockDim.x )
    atomicAdd( &( mdd.dev_deriv[ i_spec ] ), shr_dev_deriv[ i_spec ] );

}

