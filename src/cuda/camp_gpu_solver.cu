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

//Gpu hardware info
int max_n_gpu_thread;
int max_n_gpu_blocks;

//General data
double *deriv_gpu_data;
double *deriv_cpu;
double *jac_gpu_data;
double *jac_cpu;
size_t deriv_size;
size_t jac_size;
size_t state_size;
size_t env_size;
size_t rate_constants_size;
size_t rate_constants_idx_size;

bool few_data = 0;
bool implemented_all = true;
int *int_pointer;
int *int_pointer_gpu;
double *double_pointer;
double *double_pointer_gpu;
double *state_gpu;
double *state_cpu;
double *env_gpu;
double *rate_constants_gpu;
double *rate_constants_cpu;
int *rate_constants_idx_gpu;

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
void solver_new_gpu_cu(int n_dep_var,
                       int n_state_var, int n_rxn,
                       int n_rxn_int_param, int n_rxn_float_param, int n_rxn_env_param,
                       int n_cells) {

  //Lengths
  state_size = n_state_var*n_cells * sizeof(double);
  deriv_size = n_dep_var*n_cells * sizeof(double);
  env_size = PMC_NUM_ENV_PARAM_*n_cells * sizeof(double); //Temp and pressure
  rate_constants_size = n_rxn_env_param * n_cells * sizeof(double);
  rate_constants_idx_size = (n_rxn+1) * sizeof(int);

  //Detect if we are working with few data values
  if (n_dep_var*n_cells < DATA_SIZE_LIMIT_OPT){
    few_data = 1;
  }
 //few_data = 1;

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
  cudaMalloc((void **) &deriv_gpu_data, deriv_size);
  cudaMalloc((void **) &state_gpu, state_size);
  cudaMalloc((void **) &env_gpu, env_size);
  cudaMalloc((void **) &rate_constants_gpu, rate_constants_size);
  cudaMalloc((void **) &rate_constants_idx_gpu, rate_constants_idx_size);

  //GPU allocation few data on pinned memory
  if(few_data){
    //Notice auxiliar variables are created because we
    // can't pin directly variables initialized before
    cudaMallocHost((void**)&rate_constants_cpu, rate_constants_size);
    cudaMallocHost((void**)&deriv_cpu, deriv_size);
  }

  // Warning if exceeding GPU limits
  if( n_blocks > max_n_gpu_blocks){
    printf("\nWarning: More blocks assigned: %d than maximum block numbers: %d",
           n_blocks, max_n_gpu_blocks);
  }

#ifdef PMC_DEBUG_PRINT
  print_gpu_specs();
#endif

}

/** \brief Allocate Jacobian on GPU
*
* \param n_jac_elem Number of data elements on the jacobian
* \param n_cells Number of cells to compute
*/
void allocate_jac_gpu(int n_jac_elem, int n_cells){

  jac_size = n_jac_elem * n_cells * sizeof(double);
  cudaMalloc((void **) &jac_gpu_data, jac_size);

  if(few_data){
    cudaMallocHost((void**)&jac_cpu, jac_size);
  }

}

/** \brief Set reaction data on GPU prepared structure. RXN data is divided
 * into two different matrix, per double and int data respectively. Matrix are
 * reversed to improve memory access on GPU.
 *
 * \param md Pointer to the model data
 */

void solver_set_rxn_data_gpu(ModelData *model_data) {

  int n_rxn = model_data->n_rxn;
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
      implemented_all=false;
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

  //todo checkout chem_mod after 99 merge, save the changed files of this branch (like gpu) copy into checkout chem_mod

  //Add a for to search the biggest distance int_max_length (ptrs[i] - ptrs[i-1]

  //Total lengths of rxn structure
  unsigned int rxn_int_length=n_rxn*int_max_length;
  unsigned int rxn_double_length=n_rxn*double_max_length;

  //Allocate int and double rxn data separately
  //Add -1 to avoid access and have a square matrix
  int_pointer = (int *) malloc(rxn_int_length * sizeof(int));
  memset(int_pointer, -1, rxn_int_length * sizeof(int));

  //Add 0 to avoid access and have a square matrix
  double_pointer = (double*)calloc(rxn_double_length, sizeof(double));

  //GPU allocation
  cudaMalloc((void **) &int_pointer_gpu, rxn_int_length * sizeof(int));
  cudaMalloc((void **) &double_pointer_gpu, rxn_double_length * sizeof(double));

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

  int rate_constants_idx_aux[n_rxn];

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
    //Todo update on main code the rate_constants to read consecutively in cpu
    rate_constants_idx_aux[i_rxn] = model_data->rxn_env_idx[i_pos];
  }

  //Save data to GPU
  HANDLE_ERROR(cudaMemcpy(int_pointer_gpu, int_pointer, rxn_int_length*sizeof(int), cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(double_pointer_gpu, double_pointer, rxn_double_length*sizeof(double), cudaMemcpyHostToDevice));

  //Set rate_constants-idx
  HANDLE_ERROR(cudaMemcpy(rate_constants_idx_gpu, rate_constants_idx_aux, rate_constants_idx_size, cudaMemcpyHostToDevice));

}

/** \brief GPU function: Solve derivative
 *
 * \param state_init Pointer to first value of state array
 * \param deriv_init Pointer to first value of derivative array
 * \param time_step Current time step being computed (s)
 * \param deriv_length_cell Derivative length for one cell
 * \param state_size_cell Derivative length for one cell
 * \param n_rxn Number of reactions to include
 * \param n_cells_gpu Number of cells to compute
 * \param int_pointer Pointer to integer reaction data
 * \param double_pointer Pointer to double reaction data
 * \param rate_constants_init Pointer to first value of reaction rates
 */
__global__ void solveDerivative(double *state_init, double *deriv_init,
                                double time_step, int deriv_length_cell, int state_size_cell,
                                int rate_constants_size_cell,
                                int n_rxn, int n_cells, int *int_pointer, double *double_pointer,
                                double *rate_constants_init, int *rate_constants_idx) //Interface CPU/GPU
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
    //TODO: reorder rate_constants (first n_cells) to make this efficient (atm is worst for large n_cells)
    //int i_cell=index%n_cells;
    //int i_rxn=index/n_cells;

    //Get indices of each reaction
    double *float_data = (double *) &(((double *) double_pointer)[i_rxn]);
    int *int_data = (int *) &(((int *) int_pointer)[i_rxn]); //Same indices for each cell
    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *deriv_data = &( deriv_init[deriv_length_cell*i_cell]);
    double *state = &( state_init[state_size_cell*i_cell]);

    //Get indices for rates
    double *rate_constants = &(rate_constants_init
            [rate_constants_size_cell*i_cell+rate_constants_idx[i_rxn]]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(rate_constants,
                state, deriv_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(rate_constants,
                                             state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(rate_constants,
                                             state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(rate_constants,
                                                state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(rate_constants,
                state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(rate_constants,
                                            state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(rate_constants,
                                                    state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(rate_constants,
        //        state, deriv_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(rate_constants,
                                              state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(md, rate_constants,
        //        state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(rate_constants,
                                        state, deriv_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(rate_constants,
                                                  state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
    }
    __syncthreads();
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
  double *rate_constants = model_data->rxn_env_data;

//TODO: not working with new arrhenius test FOR SOME REASON SEG FAULT
//Pa ARREGLARLO poner deriv size y todo eso en un struct gpu y asi el multi cell y one cell core tienen su struct y sus cosas
//(parece que al ser global variable pues no lo pilla bien y no sobreescribe x algun motivo)

/*

  //Faster, use for few values
  if (few_data){
    //This method of passing them as a function parameter has a theoric maximum of 4kb of data
    state_gpu= state;
    rate_constants_gpu= rate_constants;
  }
    //Slower, use for large values
  else{
    HANDLE_ERROR(cudaMemcpy(state_gpu, state, state_size, cudaMemcpyHostToDevice));
    //todo rate_constants only change each time_step, not each deriv iteration
    HANDLE_ERROR(cudaMemcpy(rate_constants_gpu, rate_constants, rate_constants_size, cudaMemcpyHostToDevice));
  }

  HANDLE_ERROR(cudaMemset(deriv_gpu_data, 0.0, deriv_size));


  cudaDeviceSynchronize();// Secure cuda synchronization

  //Use pinned memory for few values
  if (!few_data){
    HANDLE_ERROR(cudaMemcpy(deriv_cpu, deriv_gpu_data, deriv_size, cudaMemcpyDeviceToHost));
    memcpy(deriv_data, deriv_cpu, deriv_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(deriv_data, deriv_gpu_data, deriv_size, cudaMemcpyDeviceToHost));
  }

  printf("holaa %d, %d\n", deriv_size/sizeof(double), n_cells);

  //TODO Calculate on CPU non gpu-translated type reactions (HL & SIMPOL phase transfer) (or wait till v2.0 with C++)

 */

}

/** \brief GPU function: Solve jacobian
 *
 * \param state_init Pointer to first value of state array
 * \param jac_init Pointer to first value of jacobian array
 * \param time_step Current time step being computed (s)
 * \param jac_length_cell jacobian length for one cell
 * \param state_size_cell jacobian length for one cell
 * \param n_rxn Number of reactions to include
 * \param n_cells_gpu Number of cells to compute
 * \param int_pointer Pointer to integer reaction data
 * \param double_pointer Pointer to double reaction data
 * \param rate_constants_init Pointer to first value of reaction rates
 */
//TODO: fix jacGPU once matt let me modify rxn_solver and reduce jacobians (in v2.0 or before if needed)
__global__ void solveJacobian(double *state_init, double *jac_init,
                              double time_step, int jac_length_cell, int state_size_cell,
                              int rate_constants_size_cell, int n_rxn,
                              int n_cells, int *int_pointer, double *double_pointer,
                              double *rate_constants_init, int *rate_constants_idx) //Interface CPU/GPU
{
  //Get thread id
  int index = blockIdx.x * blockDim.x + threadIdx.x;

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
    int *rxn_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *jac_data = &( jac_init[jac_length_cell*i_cell]);
    double *state = &( state_init[state_size_cell*i_cell]);

    //Get indices for rates
    double *rate_constants = &(rate_constants_init
    [rate_constants_size_cell*i_cell+rate_constants_idx[i_rxn]]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_jac_contrib(rate_constants,
        //        state, jac_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_ARRHENIUS :
        //rxn_gpu_arrhenius_calc_jac_contrib(rate_constants,
        //                                   state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_jac_contrib(rate_constants,
                                           state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_jac_contrib(rate_constants,
                                              state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        //rxn_gpu_condensed_phase_arrhenius_calc_jac_contrib(rate_constants,
        //        state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_jac_contrib(rate_constants,
                                          state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_jac_contrib(rate_constants,
                                                  state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_jac_contrib(rate_constants,
        //        state, jac_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_jac_contrib(rate_constants,
                                            state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_jac_contrib(rate_constants,
        //        state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_jac_contrib(rate_constants,
                                      state, jac_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_jac_contrib(rate_constants,
                                                state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
    }
    __syncthreads();
  }
}

/** \brief Calculate the Jacobian on GPU
 *
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current model time step (s)
 */

void rxn_calc_jac_gpu(ModelData *model_data, SUNMatrix jac, realtype time_step) {

  // Get a pointer to the jacobian data
  double *jac_data = SM_DATA_S(jac);
  int n_cells = model_data->n_cells;
  int n_rxn = model_data->n_rxn;
  int n_rxn_threads = n_rxn*n_cells; //Reaction group per number of repetitions/cells
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);
  double *state = model_data->total_state;
  double *rate_constants = model_data->rxn_env_data;

  /*

  //Faster, use for few values
  if (few_data){
    //This method of passing them as a function parameter has a theoric maximum of 4kb of data
    state_gpu= state;
    rate_constants_gpu= rate_constants;
  }
    //Slower, use for large values
  else{
    HANDLE_ERROR(cudaMemcpy(state_gpu, state, state_size, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(rate_constants_gpu, rate_constants, rate_constants_size, cudaMemcpyHostToDevice));
  }

  HANDLE_ERROR(cudaMemset(jac_gpu_data, 0, jac_size));

  solveJacobian << < (n_blocks), max_n_gpu_thread >> >
    (state_gpu, jac_gpu_data, time_step, model_data->n_per_cell_rxn_jac_elem,
    model_data->n_per_cell_state_var, model_data->n_rxn_env_data,
    n_rxn, n_cells, int_pointer_gpu, double_pointer_gpu, rate_constants_gpu, rate_constants_idx_gpu);

  cudaDeviceSynchronize();// Secure cuda synchronization

  //Use pinned memory for few values
  if (few_data){
    HANDLE_ERROR(cudaMemcpy(jac_cpu, jac_gpu_data, jac_size, cudaMemcpyDeviceToHost));
    memcpy(jac_data, jac_cpu, jac_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(jac_data, jac_gpu_data, jac_size, cudaMemcpyDeviceToHost));
  }

   */

}


/** \brief Free GPU data structures
 */
void free_gpu_cu() {

printf("free");

  HANDLE_ERROR(cudaFree(int_pointer_gpu));
  HANDLE_ERROR(cudaFree(double_pointer_gpu));
  HANDLE_ERROR(cudaFree(deriv_gpu_data));
  //HANDLE_ERROR(cudaFree(jac_gpu_data));

  if(few_data){
  }
  else{
    HANDLE_ERROR(cudaFree(state_gpu));
    HANDLE_ERROR(cudaFree(rate_constants_gpu));
    HANDLE_ERROR(cudaFree(rate_constants_idx_gpu));
  }
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
