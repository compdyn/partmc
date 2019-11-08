/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Interface Host-Device (CPU-GPU) to compute reaction-specific functions on GPU
 *
 */
extern "C" {
#include "camp_gpu_solver.h"
#include "rxns_gpu.h"

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5
#define RXN_HL_PHASE_TRANSFER 6
#define RXN_AQUEOUS_EQUILIBRIUM 7
#define RXN_SIMPOL_PHASE_TRANSFER 10
#define RXN_CONDENSED_PHASE_ARRHENIUS 11
#define RXN_FIRST_ORDER_LOSS 12
#define RXN_EMISSION 13
#define RXN_WET_DEPOSITION 14

#define STREAM_RXN_ENV_GPU 0
#define STREAM_ENV_GPU 1
#define STREAM_DERIV_GPU 2

//GPU async stream related variables to ensure robustness
int n_solver_objects=0; //Number of solver_new_gpu calls
//cudaStream_t *stream_gpu; //GPU streams to async computation/data movement
int n_streams = 16;
//todo: this is a supposition, should be fixed with some way of know how many solver
// //objects we have and how many streams we need allocate in stream array
//(problem is cudaStream_t not found in .h like modelData )
int n_streams_limit = 128;

//Gpu hardware info
int max_n_gpu_thread;
int max_n_gpu_blocks;

//Debug info
int counterDeriv;       // Total calls to f()
int counterJac;         // Total calls to Jac()
clock_t timeDeriv;      // Compute time for calls to f()
clock_t timeJac;        // Compute time for calls to Jac()
clock_t timeDerivKernel; // Compute time for calls to f() kernel
clock_t timeDerivSend;
clock_t timeDerivReceive;
clock_t timeDerivCPU;
clock_t t1;             //Auxiliar time counter
clock_t t3;


static void HandleError(cudaError_t err,
                        const char *file,
                        int line) {
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err),
           file, line);
    exit(EXIT_FAILURE);
  }
}

static void HandleError2(const char *file,
                         int line) {
  cudaError_t err;
  err=cudaGetLastError();
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err),
           file, line);    exit(EXIT_FAILURE);
  }
}

/** \brief Allocate GPU solver variables
 *
 * \param n_dep_var number of solver variables per grid cell
 * \param n_state_var Number of variables on the state array per grid cell
 * \param n_rxn Number of reactions to include
 * \param n_rxn_int_param Total number of integer reaction parameters
 * \param n_rxn_float_param Total number of floating-point reaction parameters
 * \param n_cells Number of grid cells to solve simultaneously
 */
void solver_new_gpu_cu(ModelData *model_data, int n_dep_var,
                       int n_state_var, int n_rxn,
                       int n_rxn_int_param, int n_rxn_float_param, int n_rxn_env_param,
                       int n_cells) {
  //TODO: Select what % of data we want to compute on GPU simultaneously with CPU remaining %

  //Lengths
  model_data->state_size = n_state_var * n_cells * sizeof(double);
  model_data->deriv_size = n_dep_var * n_cells * sizeof(double);
  model_data->env_size = PMC_NUM_ENV_PARAM_ * n_cells * sizeof(double); //Temp and pressure
  model_data->rxn_env_data_size = n_rxn_env_param * n_cells * sizeof(double);
  model_data->rxn_env_data_idx_size = (n_rxn+1) * sizeof(int);
  model_data->small_data = 0;
  model_data->implemented_all = true;
  counterDeriv=0;
  counterJac=0;
  timeDeriv=0;
  timeDerivKernel=0;
  timeDerivSend=0;
  timeDerivReceive=0;
  timeJac=0;

  //Allocate streams array and update variables related to streams
  model_data->model_data_id = n_solver_objects;
  //if(n_solver_objects==0){
    //stream_gpu = (cudaStream_t *)malloc(n_streams_limit * sizeof(cudaStream_t));
    model_data->stream_gpu = (cudaStream_t *)malloc(n_streams * sizeof(cudaStream_t));
  //}
  n_solver_objects++;

  //Detect if we are working with few data values
  if (n_dep_var*n_cells < DATA_SIZE_LIMIT_OPT){
    model_data->small_data = 1;
  }

  //Set working GPU: we have 4 gpu available on power9. as default, it should be assign to gpu 0
  int device=0;
  cudaSetDevice(device);

  //Set GPU properties
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);

  //Set max threads without triggering too many resources error
  max_n_gpu_thread = prop.maxThreadsPerBlock/2;
  max_n_gpu_blocks = prop.maxGridSize[1];
  int n_blocks = (n_rxn + max_n_gpu_thread - 1) / max_n_gpu_thread;

  //GPU allocation
  cudaMalloc((void **) &model_data->deriv_gpu_data, model_data->deriv_size);
  cudaMalloc((void **) &model_data->state_gpu, model_data->state_size);
  cudaMalloc((void **) &model_data->env_gpu, model_data->env_size);
  cudaMalloc((void **) &model_data->rxn_env_data_gpu, model_data->rxn_env_data_size);
  cudaMalloc((void **) &model_data->rxn_env_data_idx_gpu, model_data->rxn_env_data_idx_size);

  //GPU allocation few data on pinned memory
  if(model_data->small_data){
    //Notice auxiliar variables are created because we
    // can't pin directly variables initialized before
    cudaMallocHost((void**)&model_data->deriv_aux, model_data->deriv_size);
  }
  else{
    model_data->deriv_aux = (realtype *)malloc(model_data->deriv_size);
  }

  printf("small_data:%d\n", model_data->small_data);

  //GPU create streams
  for (int i = 0; i < n_streams; ++i)
    HANDLE_ERROR( cudaStreamCreate(&model_data->stream_gpu[i]) );

  // Warning if exceeding GPU limits
  if( n_blocks > max_n_gpu_blocks){
    printf("\nWarning: More blocks assigned: %d than maximum block numbers: %d",
           n_blocks, max_n_gpu_blocks);
  }

#ifdef PMC_DEBUG_PRINT
  print_gpu_specs();
#endif

}

/** \brief Set reaction data on GPU prepared structure. RXN data is divided
 * into two different matrix, per double and int data respectively. Matrix are
 * reversed to improve memory access on GPU.
 *
 * \param md Pointer to the model data
 */

void solver_set_rxn_data_gpu(ModelData *model_data) {

  int n_rxn = model_data->n_rxn;
  int n_cells = model_data->n_cells;
  unsigned int int_max_length = 0;
  unsigned int double_max_length = 0;

  //RXN lengths
  unsigned int int_lengths[n_rxn];
  unsigned int double_lengths[n_rxn];

  //Number of extra values added to square matrix(zeros and -1's)
  unsigned int n_zeros[n_rxn];

  //Position on the matrix for each row
  unsigned int rxn_position[n_rxn];

  //Get lengths for int and double arrays
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

    // Set a WARNING if the reaction is not implemented yet on GPU
    bool implemented = false;
    int rxn_type = model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]];

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        implemented = false;
        break;
      case RXN_ARRHENIUS :
        implemented = true;
        break;
      case RXN_CMAQ_H2O2 :
        implemented = true;
        break;
      case RXN_CMAQ_OH_HNO3 :
        implemented = true;
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        implemented = false;
        break;
      case RXN_EMISSION :
        implemented = true;
        break;
      case RXN_FIRST_ORDER_LOSS :
        implemented = true;
        break;
      case RXN_HL_PHASE_TRANSFER :
        implemented = false;
        break;
      case RXN_PHOTOLYSIS :
        implemented = true;
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        implemented = false;
        break;
      case RXN_TROE :
        implemented = true;
        break;
      case RXN_WET_DEPOSITION :
        implemented = true;
        break;
    }
    if(!implemented){
#ifdef FAILURE_DETAIL
      printf("WARNING: Reaction type %d is not fully implemented on GPU. Computing on CPU...\n", rxn_type);
#endif
      model_data->implemented_all=false;
    }

    //Get RXN lengths
    int_lengths[i_rxn] = model_data->rxn_int_indices[i_rxn+1] - model_data->rxn_int_indices[i_rxn];
    double_lengths[i_rxn] = model_data->rxn_float_indices[i_rxn+1] - model_data->rxn_float_indices[i_rxn];

    //Update max size
    if(int_lengths[i_rxn]>int_max_length) int_max_length=int_lengths[i_rxn];
    if(double_lengths[i_rxn]>double_max_length) double_max_length=double_lengths[i_rxn];

    //Set initial position
    rxn_position[i_rxn] = i_rxn;

  }

  //Add a for to search the biggest distance int_max_length (ptrs[i] - ptrs[i-1]

  //Total lengths of rxn structure
  unsigned int rxn_int_length=n_rxn*int_max_length;
  unsigned int rxn_double_length=n_rxn*double_max_length;

  //Allocate int and double rxn data separately
  //Add -1 to avoid access and have a square matrix
  int *int_pointer = (int *) malloc(rxn_int_length * sizeof(int));
  memset(int_pointer, -1, rxn_int_length * sizeof(int));

  //Add 0 to avoid access and have a square matrix
  double *double_pointer = (double*)calloc(rxn_double_length, sizeof(double));

  //GPU allocation
  cudaMalloc((void **) &model_data->int_pointer_gpu, rxn_int_length * sizeof(int));
  cudaMalloc((void **) &model_data->double_pointer_gpu, rxn_double_length * sizeof(double));

  //Update number of zeros added on each reaction
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++)
    n_zeros[i_rxn] = (int_max_length - int_lengths[i_rxn]) +
                     (double_max_length - double_lengths[i_rxn]);

  //NOTE: no improvement on doing the sorting or not for gpu seems.
  //Sort by lengths
  //BubbleSort RXN by ascendant number of zeros for performance reasons
  //Fix reordered rxn give wrong values
  //bubble_sort_gpu(n_zeros, rxn_position, n_rxn);

  //Copy into gpu rxn data
  //Follows the rxn_position order
  //Rxn matrix is reversed to improve memory access on GPU
  //Matrix order is [int_length][n_rxn]

  int rxn_env_data_idx_aux[n_rxn];

  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    int i_pos=rxn_position[i_rxn];//i_rxn;//rxn_position[i_rxn];//for bubblesort
    for (int j = 0; j < int_lengths[i_pos]; j++){
      int *rxn_int_data = &(model_data->rxn_int_data[model_data->rxn_int_indices[i_pos]]);
      int_pointer[n_rxn*j + i_rxn] = rxn_int_data[j];
    }
    for (int j = 0; j < double_lengths[i_pos]; j++) {
      double *rxn_float_data = &(model_data->rxn_float_data[model_data->rxn_float_indices[i_pos]]);
      double_pointer[n_rxn*j + i_rxn] = rxn_float_data[j];
    }
    //Reorder the rate indices
    //Todo update on main code the rxn_env_data to read consecutively in cpu
    rxn_env_data_idx_aux[i_rxn] = model_data->rxn_env_idx[i_pos];
  }

  //Save data to GPU
  HANDLE_ERROR(cudaMemcpy(model_data->int_pointer_gpu, int_pointer, rxn_int_length*sizeof(int), cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(model_data->double_pointer_gpu, double_pointer, rxn_double_length*sizeof(double), cudaMemcpyHostToDevice));

  //Set rxn_env_data-idx
  HANDLE_ERROR(cudaMemcpy(model_data->rxn_env_data_idx_gpu, rxn_env_data_idx_aux, model_data->rxn_env_data_idx_size, cudaMemcpyHostToDevice));

  //Allocate jacobian
  model_data->jac_size = model_data->n_per_cell_solver_jac_elem * n_cells * sizeof(double);
  cudaMalloc((void **) &model_data->jac_gpu_data, model_data->jac_size);

  if(model_data->small_data){
    cudaMallocHost((void**)&model_data->jac_aux, model_data->jac_size);
  }

  free(int_pointer);
  free(double_pointer);

}


void rxn_update_env_state_gpu(ModelData *model_data){

  // Get a pointer to the derivative data
  int n_cells = model_data->n_cells;
  int n_rxn = model_data->n_rxn;
  int n_rxn_threads = n_rxn*n_cells; //Reaction group per number of repetitions/cells
  double *state = model_data->total_state;
  double *rxn_env_data = model_data->rxn_env_data;
  double *env = model_data->total_env;
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);

  //Faster, use for few values
  if (model_data->small_data){
    //This method of passing them as a function parameter has a theoric maximum of 4kb of data
    model_data->rxn_env_data_gpu= rxn_env_data;
    model_data->env_gpu= env;
  }
    //Slower, use for large values
  else{
//TODO: seems no improvement, do some profiling

/*HANDLE_ERROR(cudaMemcpyAsync(model_data->rxn_env_data_gpu, rxn_env_data,
            model_data->rxn_env_data_size, cudaMemcpyHostToDevice,
                 model_data->stream_gpu[STREAM_RXN_ENV_GPU]));
    HANDLE_ERROR(cudaMemcpyAsync(model_data->env_gpu, env, model_data->env_size,
            cudaMemcpyHostToDevice,
                 model_data->stream_gpu[STREAM_ENV_GPU]));*/

  }

}

/** \brief GPU function: Solve derivative
 *
 * \param state_init Pointer to first value of state array
 * \param deriv_init Pointer to first value of derivative array
 * \param time_step Current time step being computed (s)
 * \param deriv_length_cell Derivative length for one cell
 * \param model_data->state_size_cell Derivative length for one cell
 * \param n_rxn Number of reactions to include
 * \param n_cells_gpu Number of cells to compute
 * \param model_data->int_pointer Pointer to integer reaction data
 * \param model_data->double_pointer Pointer to double reaction data
 * \param rxn_env_data_init Pointer to first value of reaction rates
 */
__global__ void solveDerivative(double *state_init, double *deriv_init,
                                double time_step, int deriv_length_cell, int state_size_cell,
                                int rxn_env_data_size_cell, int n_rxn, int n_cells,
                                int *int_pointer, double *double_pointer,
                                double *rxn_env_data_init, int *rxn_env_data_idx,
                                double *env_init) //Interface CPU/GPU
{
  //Get thread id
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  //Maximum number of threads to compute all reactions
  if(index < n_rxn*n_cells){

    //Thread index for deriv and state,
    // till we don't finish all reactions of a cell, we stay on same index
    int i_cell=index/n_rxn;
    int i_rxn=index%n_rxn;

    //Another option: compute first the cells and then the reactions (seems working fine but no speedup)
    //TODO: reorder rxn_env_data (first n_cells) to make this efficient (atm is worst for large n_cells)
    //int i_cell=index%n_cells;
    //int i_rxn=index/n_cells;

    //Get indices of each reaction
    double *rxn_float_data = (double *) &(((double *) double_pointer)[i_rxn]);
    int *int_data = (int *) &(((int *) int_pointer)[i_rxn]); //Same indices for each cell
    int rxn_type = int_data[0];
    int *rxn_int_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *deriv_data = &( deriv_init[deriv_length_cell*i_cell]);
    double *state = &( state_init[state_size_cell*i_cell]);

    ModelData model_data;
    model_data.grid_cell_state = &( state_init[state_size_cell*i_cell]);
    model_data.grid_cell_env = &( env_init[PMC_NUM_ENV_PARAM_*i_cell]);
    model_data.n_rxn = n_rxn;

    //Get indices for rates
    double *rxn_env_data = &(rxn_env_data_init
            [rxn_env_data_size_cell*i_cell+rxn_env_data_idx[i_rxn]]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                                       rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
        //                                             rxn_float_data, rxn_env_data,time_stepn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(md, rxn_env_data,
        //        state, deriv_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                        rxn_float_data, rxn_env_data,time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(&model_data, deriv_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data,time_step);
        break;
    }
    __syncthreads();
  }

}

static void print_derivative_3(N_Vector deriv) {
  // printf(" deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
    printf(" deriv: % -le", NV_DATA_S(deriv)[i]);
    printf(" index: %d \n", i);
  }
}

/** \brief Calculate the time derivative \f$f(t,y)\f$ on GPU
 *
 * \param model_data Pointer to the model data
 * \param deriv NVector to hold the calculated vector
 * \param time_step Current model time step (s)
 */
void rxn_calc_deriv_gpu(ModelData *model_data, N_Vector deriv, realtype time_step) {

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int n_cells = model_data->n_cells;
  int n_rxn = model_data->n_rxn;
  int n_rxn_threads = n_rxn*n_cells; //Reaction group per number of repetitions/cells
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);
  double *state = model_data->total_state;
  double *rxn_env_data = model_data->rxn_env_data;
  double *env = model_data->total_env;

  //TODO: cudaMemcpyAsync : Since cudamemcpy is taking like 6 times of gpu compute time,
  //use async copy of parts of the code to go copying other parts meanwhile computing
  //(a carefully tuned CUDA program that uses streams and cudaMemcpyAsync to efficiently overlap execution with data
  // transfers may very well perform)

  t1 = clock();


  //TODO: overlap state data transfer with kernel exec dividing into multiple async streams
  //Faster, use for few values
  if (model_data->small_data){
    //This method of passing them as a function parameter has a theoric maximum of 4kb of data
    model_data->state_gpu= state;
  }
    //Slower, use for large values
  else{
    HANDLE_ERROR(cudaMemcpy(model_data->state_gpu, state, model_data->state_size, cudaMemcpyHostToDevice));
    //todo rxn_env_data only change each time_step, not each deriv iteration

    HANDLE_ERROR(cudaMemcpy(model_data->env_gpu, env, model_data->env_size, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(model_data->rxn_env_data_gpu, rxn_env_data, model_data->rxn_env_data_size, cudaMemcpyHostToDevice));

  }

  //todo this async later
  //HANDLE_ERROR(cudaMemset(model_data->deriv_gpu_data, 0.0, model_data->deriv_size));
  timeDerivSend += (clock() - t1);
  clock_t t2 = clock();

  //TODO: execute this asyncrhonous and let CPU thinks be computed in the meanwhile
  //TODO: compute this and jac in different streams (maybe we can parallelism them)
  cudaDeviceSynchronize();
  solveDerivative << < (n_blocks), max_n_gpu_thread >> >
   (model_data->state_gpu, model_data->deriv_gpu_data, time_step, model_data->n_per_cell_dep_var,
   model_data->n_per_cell_state_var, model_data->n_rxn_env_data,
   n_rxn, n_cells, model_data->int_pointer_gpu, model_data->double_pointer_gpu,
   model_data->rxn_env_data_gpu, model_data->rxn_env_data_idx_gpu, model_data->env_gpu);

  cudaDeviceSynchronize();
  timeDerivKernel += (clock() - t2);
  t3 = clock();

  //TODO: make this memcpyasync and add contributions from this to cpu contrib (Create aux deriv_data to store gpu memcpy and later add to deriv_data)
  //Use pinned memory for few values
  if (model_data->small_data){
    HANDLE_ERROR(cudaMemcpy(model_data->deriv_aux, model_data->deriv_gpu_data, model_data->deriv_size, cudaMemcpyDeviceToHost));
    memcpy(deriv_data, model_data->deriv_aux, model_data->deriv_size);
  }
  else {
    //HANDLE_ERROR(cudaMemcpyAsync(model_data->deriv_aux, model_data->deriv_gpu_data,
    //        model_data->deriv_size, cudaMemcpyDeviceToHost,
    //        model_data->stream_gpu[STREAM_DERIV_GPU]));
HANDLE_ERROR(cudaMemcpy(deriv_data, model_data->deriv_gpu_data, model_data->deriv_size, cudaMemcpyDeviceToHost));
  }

  timeDerivReceive += (clock() - t3);
  timeDeriv += (clock() - t1);
  t3 = clock();

  //((double)t)/CLOCKS_PER_SEC;

}

void rxn_fusion_deriv_gpu(ModelData *model_data, N_Vector deriv) {


  timeDerivCPU += (clock() - t1);

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);

  cudaDeviceSynchronize();

  HANDLE_ERROR(cudaMemsetAsync(model_data->deriv_gpu_data, 0.0,
          model_data->deriv_size, model_data->stream_gpu[STREAM_DERIV_GPU]));




  if (model_data->small_data){
  }
  else {
    for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
      //Add to deriv the auxiliar contributions from gpu
      //deriv_data[i] += model_data->deriv_aux[i];
    }
  }

}


#ifdef PMC_USE_SUNDIALS
void rxn_calc_deriv_cpu(ModelData *model_data, double *deriv_data,
                    realtype time_step) {

  //clock_t t = clock();

  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);
    double *rxn_env_data =
        &(model_data->grid_cell_rxn_env_data[model_data->rxn_env_idx[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(model_data, deriv_data,
                                                   rxn_int_data, rxn_float_data,
                                                   rxn_env_data, time_step);
        break;
      case RXN_ARRHENIUS:
        rxn_gpu_arrhenius_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                         rxn_float_data, rxn_env_data,
                                         time_step);
        break;
      case RXN_CMAQ_H2O2:
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                         rxn_float_data, rxn_env_data,
                                         time_step);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(model_data, deriv_data,
                                            rxn_int_data, rxn_float_data,
                                            rxn_env_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
            model_data, deriv_data, rxn_int_data, rxn_float_data, rxn_env_data,
            time_step);
        break;
      case RXN_EMISSION:
        rxn_gpu_emission_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                        rxn_float_data, rxn_env_data,
                                        time_step);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_gpu_first_order_loss_calc_deriv_contrib(model_data, deriv_data,
                                                rxn_int_data, rxn_float_data,
                                                rxn_env_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER:
        //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(model_data, deriv_data,
        //                                         rxn_int_data, rxn_float_data,
        //                                         rxn_env_data, time_step);
        break;
      case RXN_PHOTOLYSIS:
        rxn_gpu_photolysis_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                          rxn_float_data, rxn_env_data,
                                          time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
        //   model_data, deriv_data, rxn_int_data, rxn_float_data, rxn_env_data,
        //    time_step);
        break;
      case RXN_TROE:
        rxn_gpu_troe_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                    rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_WET_DEPOSITION:
        rxn_gpu_wet_deposition_calc_deriv_contrib(model_data, deriv_data,
                                              rxn_int_data, rxn_float_data,
                                              rxn_env_data, time_step);
        break;
    }
  }

  //timeDeriv += (clock()- t);

}
#endif

/** \brief GPU function: Solve jacobian
 *
 * \param state_init Pointer to first value of state array
 * \param jac_init Pointer to first value of jacobian array
 * \param time_step Current time step being computed (s)
 * \param jac_length_cell jacobian length for one cell
 * \param model_data->state_size_cell jacobian length for one cell
 * \param n_rxn Number of reactions to include
 * \param n_cells_gpu Number of cells to compute
 * \param model_data->int_pointer Pointer to integer reaction data
 * \param model_data->double_pointer Pointer to double reaction data
 * \param rxn_env_data_init Pointer to first value of reaction rates
 */
__global__ void solveJacobian(double *state_init, double *jac_init,
                              double time_step, int jac_length_cell, int state_size_cell,
                              int rxn_env_data_size_cell, int n_rxn,
                              int n_cells, int *int_pointer, double *double_pointer,
                              double *rxn_env_data_init, int *rxn_env_data_idx) //Interface CPU/GPU
{
  //Get thread id
  /*int index = blockIdx.x * blockDim.x + threadIdx.x;

  //Maximum number of threads to compute all reactions
  if(index < n_rxn*n_cells){

    //Thread index for jac and state,
    // till we don't finish all reactions of a cell, we stay on same index
    int i_cell=index/n_rxn;
    int i_rxn=index%n_rxn;

    //Get indices of each reaction
    int *int_data = (int *) &(((int *) int_pointer)[i_rxn]); //Same indices for each cell
    double *float_data = (double *) &(((double *) double_pointer)[i_rxn]);
    int rxn_type = int_data[0];
    int *rxn_int_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *jac_data = &( jac_init[jac_length_cell*i_cell]);
    double *state = &( state_init[state_size_cell*i_cell]);

    //Get indices for rates
    double *rxn_env_data = &(rxn_env_data_init
    [rxn_env_data_size_cell*i_cell+rxn_env_data_idx[i_rxn]]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_jac_contrib(rxn_env_data,
        //        state, jac_data, rxn_int_data, rxn_float_data, time_step, n_rxn);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_jac_contrib(rxn_env_data,
                                           state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_jac_contrib(rxn_env_data,
                                           state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_jac_contrib(rxn_env_data,
                                              state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        //rxn_gpu_condensed_phase_arrhenius_calc_jac_contrib(rxn_env_data,
        //        state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_jac_contrib(rxn_env_data,
                                          state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_jac_contrib(rxn_env_data,
                                                  state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_jac_contrib(rxn_env_data,
        //        state, jac_data, rxn_int_data, rxn_float_data, time_step, n_rxn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_jac_contrib(rxn_env_data,
                                            state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_jac_contrib(rxn_env_data,
        //        state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_jac_contrib(rxn_env_data,
                                      state, jac_data, rxn_int_data, rxn_float_data, time_step, n_rxn);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_jac_contrib(rxn_env_data,
                                                state, jac_data, rxn_int_data, rxn_float_data, time_step,n_rxn);
        break;
    }
    __syncthreads();
  }
   */

}

/** \brief Calculate the Jacobian on GPU
 *
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current model time step (s)
 */

void rxn_calc_jac_gpu(ModelData *model_data, SUNMatrix jac, realtype time_step) {

  //TODO: Fix jacobian with jac_ids...

  // Get a pointer to the jacobian data
  double *jac_data = SM_DATA_S(jac);
  int n_cells = model_data->n_cells;
  int n_rxn = model_data->n_rxn;
  int n_rxn_threads = n_rxn*n_cells; //Reaction group per number of repetitions/cells
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);
  double *state = model_data->total_state;
  double *rxn_env_data = model_data->rxn_env_data;

  //TODO: in jacobian is not need to pass the state? should be the same than deriv state?

  //Faster, use for few values
  if (model_data->small_data){
    //This method of passing them as a function parameter has a theoric maximum of 4kb of data
    model_data->state_gpu= state;
    model_data->rxn_env_data_gpu= rxn_env_data;
  }
    //Slower, use for large values
  else{
    HANDLE_ERROR(cudaMemcpy(model_data->state_gpu, state, model_data->state_size, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(model_data->rxn_env_data_gpu, rxn_env_data, model_data->rxn_env_data_size, cudaMemcpyHostToDevice));
  }

  HANDLE_ERROR(cudaMemset(model_data->jac_gpu_data, 0, model_data->jac_size));

  solveJacobian << < (n_blocks), max_n_gpu_thread >> >
    (model_data->state_gpu, model_data->jac_gpu_data, time_step, model_data->n_per_cell_rxn_jac_elem,
    model_data->n_per_cell_state_var, model_data->n_rxn_env_data,
    n_rxn, n_cells, model_data->int_pointer_gpu, model_data->double_pointer_gpu, model_data->rxn_env_data_gpu, model_data->rxn_env_data_idx_gpu);

  cudaDeviceSynchronize();// Secure cuda synchronization

  //Use pinned memory for few values
  if (model_data->small_data){
    HANDLE_ERROR(cudaMemcpy(model_data->jac_aux, model_data->jac_gpu_data, model_data->jac_size, cudaMemcpyDeviceToHost));
    memcpy(jac_data, model_data->jac_aux, model_data->jac_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(jac_data, model_data->jac_gpu_data, model_data->jac_size, cudaMemcpyDeviceToHost));
  }

}


/** \brief Free GPU data structures
 */
void free_gpu_cu(ModelData *model_data) {

  printf("timeDeriv %lf\n", (((double)timeDeriv) * 1000) / CLOCKS_PER_SEC);
  printf("timeDerivSend %lf\n", (((double)timeDerivSend) * 1000) / CLOCKS_PER_SEC);
  printf("timeDerivKernel %lf\n", (((double)timeDerivKernel) * 1000) / CLOCKS_PER_SEC);
  printf("timeDerivReceive %lf\n", (((double)timeDerivReceive) * 1000) / CLOCKS_PER_SEC);
  printf("timeDerivCPU %lf\n", (((double)timeDerivCPU) * 1000) / CLOCKS_PER_SEC);

  for (int i = 0; i < n_streams; ++i)
    HANDLE_ERROR( cudaStreamDestroy(model_data->stream_gpu[i]) );
/*

  */
  //free(model_data->jac_aux);
  HANDLE_ERROR(cudaFree(model_data->int_pointer_gpu));
  HANDLE_ERROR(cudaFree(model_data->double_pointer_gpu));
  HANDLE_ERROR(cudaFree(model_data->deriv_gpu_data));
  //HANDLE_ERROR(cudaFree(jac_gpu_data));

  if(model_data->small_data){
  }
  else{
    free(model_data->deriv_aux);
    HANDLE_ERROR(cudaFree(model_data->state_gpu));
    HANDLE_ERROR(cudaFree(model_data->env_gpu));
    HANDLE_ERROR(cudaFree(model_data->rxn_env_data_gpu));
    HANDLE_ERROR(cudaFree(model_data->rxn_env_data_idx_gpu));

  }

/*
  HANDLE_ERROR(cudaFree(int_pointer_gpu));
  HANDLE_ERROR(cudaFree(double_pointer_gpu));
  HANDLE_ERROR(cudaFree(deriv_gpu_data));
  HANDLE_ERROR(cudaFree(jac_gpu_data));

  if(small_data){
  }
  else{
    HANDLE_ERROR(cudaFree(state_gpu));
    HANDLE_ERROR(cudaFree(rxn_env_data_gpu));
    HANDLE_ERROR(cudaFree(rxn_env_data_idx_gpu));
  }
*/
}

/* Auxiliar functions */

void bubble_sort_gpu(unsigned int *n_zeros, unsigned int *rxn_position, int n_rxn){

  int tmp,s=1,i_rxn=n_rxn;

  while(s){
    s=0;
    for (int i = 1; i < i_rxn; i++) {
      //Few zeros go first
      if (n_zeros[i] < n_zeros[i - 1]) {
        //Swap positions
        tmp = rxn_position[i];
        rxn_position[i] = rxn_position[i - 1];
        rxn_position[i - 1] = tmp;

        tmp = n_zeros[i];
        n_zeros[i] = n_zeros[i - 1];
        n_zeros[i - 1] = tmp;
        s=1;
      }
    }
    i_rxn--;
  }

}

/* Prints */

void print_gpu_specs() {

  printf("GPU specifications \n");

  int nDevices;
  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Peak Memory Bandwidth (GB/s): %f\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printf("  maxGridSize: %d\n", prop.maxGridSize[1]);
    printf("  maxThreadsPerBlock: %d\n", prop.maxThreadsPerBlock);
    printf("  maxThreadsDim: %d\n", prop.maxThreadsDim[1]);
    printf("  totalGlobalMem: %zu\n", prop.totalGlobalMem);
    printf("  sharedMemPerBlock: %zu\n", prop.sharedMemPerBlock); //bytes
    printf("  multiProcessorCount: %d\n", prop.multiProcessorCount);
  }

}

}
