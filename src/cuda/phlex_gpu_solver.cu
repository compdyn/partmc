//Test
extern "C" {
#include "phlex_gpu_solver.h"
#include "rxns_gpu.h"

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

bool few_data = 0;
unsigned int *start_rxn_param;
unsigned int *dev_start_rxn_param;
unsigned int int_max_size = 0;
unsigned int double_max_size = 0;
int *int_pointer;
int *int_pointer_gpu;
double *double_pointer;
double *double_pointer_gpu;
double *state_gpu;
double *state_cpu;
double *env_gpu;
double *rate_constants_gpu;
double *rate_constants_cpu;

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
                       int n_rxn_int_param, int n_rxn_float_param,
                       int n_cells) {

  //Sizes
  size_t start_size = (n_rxn+1) * sizeof(unsigned int);
  state_size = n_state_var*n_cells * sizeof(double);
  deriv_size = n_dep_var*n_cells * sizeof(double);
  env_size = PMC_NUM_ENV_PARAM_*n_cells * sizeof(double); //Temp and pressure
  rate_constants_size = n_rxn * n_cells * sizeof(double);
  start_rxn_param = (unsigned int *) malloc(start_size); //Started indexes of rxn arrays

  //Detect if we are working with few data values
  if (n_dep_var*n_cells < DATA_SIZE_LIMIT_OPT){
    few_data = 1;
  }

  //Set working GPU: we have 4 gpu available on power9. as default, it should be assign to gpu 0
  int device=0;
  cudaSetDevice(device);

  //Set GPU properties
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  max_n_gpu_thread = prop.maxThreadsPerBlock/2; //TODO: Fix bug insufficient resources without the /2 (1024 threads)
  max_n_gpu_blocks = prop.maxGridSize[1];
  int n_blocks = (n_rxn + max_n_gpu_thread - 1) / max_n_gpu_thread;

  //GPU allocation
  cudaMalloc((void **) &dev_start_rxn_param, start_size);
  cudaMalloc((void **) &deriv_gpu_data, deriv_size);
  cudaMalloc((void **) &state_gpu, state_size);
  cudaMalloc((void **) &env_gpu, env_size);
  cudaMalloc((void **) &rate_constants_gpu, rate_constants_size);

  //GPU allocation few data
  if(few_data){
    //Pinned memory

    //Can't pin directly variables initialized before, auxiliar variables are created
    cudaMallocHost((void**)&rate_constants_cpu, rate_constants_size);
    cudaMallocHost((void**)&deriv_cpu, deriv_size);
  }

  // Warning if exceeding GPU limits

  if( n_blocks > max_n_gpu_blocks){
    printf("\nWarning: More blocks assigned: %d than maximum block numbers: %d",
           n_blocks, max_n_gpu_blocks);
  }

#ifdef PMC_DEBUG
  print_gpu_data();
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

  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);
  void *rxn_param = (void *) rxn_data;
  int *float_data = (int *) rxn_data;
  unsigned int int_size = 0;
  unsigned int int_total_size = 0;
  unsigned int double_size = 0;
  unsigned int double_total_size = 0;

  size_t start_size = (n_rxn+1) * sizeof(unsigned int); //Added 1 to simplify things
  unsigned int int_sizes[start_size];
  unsigned int double_sizes[start_size];

  //Get sizes for int and double arrays
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

    //Reaction distances between pointers rows
    start_rxn_param[i_rxn] = (unsigned int) ((int *) rxn_data - (int *) rxn_param);

    int *rxn_start = rxn_data;

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        float_data = (int *) rxn_gpu_aqueous_equilibrium_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_aqueous_equilibrium_skip((void *) rxn_data);
        break;
      case RXN_ARRHENIUS :
        float_data = (int *) rxn_gpu_arrhenius_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_arrhenius_skip((void *) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        float_data = (int *) rxn_gpu_CMAQ_H2O2_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_CMAQ_H2O2_skip((void *) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        float_data = (int *) rxn_gpu_CMAQ_OH_HNO3_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_CMAQ_OH_HNO3_skip((void *) rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        float_data = (int *) rxn_gpu_condensed_phase_arrhenius_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_condensed_phase_arrhenius_skip((void *) rxn_data);
        break;
      case RXN_EMISSION :
        float_data = (int *) rxn_gpu_emission_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_emission_skip((void *) rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        float_data = (int *) rxn_gpu_first_order_loss_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_first_order_loss_skip((void *) rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        float_data = (int*) rxn_gpu_HL_phase_transfer_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_HL_phase_transfer_skip((void *) rxn_data);
        break;
      case RXN_PDFITE_ACTIVITY :
        float_data = (int *) rxn_gpu_PDFiTE_activity_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_PDFiTE_activity_skip((void *) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        float_data = (int *) rxn_gpu_photolysis_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_photolysis_skip((void *) rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        float_data = (int *) rxn_gpu_SIMPOL_phase_transfer_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_SIMPOL_phase_transfer_skip((void *) rxn_data);
        break;
      case RXN_TROE :
        float_data = (int *) rxn_gpu_troe_int_size(
                (void *) rxn_data);
        rxn_data =(int*)rxn_gpu_troe_skip((void *) rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        float_data = (int *) rxn_gpu_wet_deposition_int_size(
                (void *) rxn_data);
        rxn_data =(int*)rxn_gpu_wet_deposition_skip((void *) rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        float_data = (int *) rxn_gpu_ZSR_aerosol_water_int_size(
                (void *) rxn_data);
        rxn_data = (int*) rxn_gpu_ZSR_aerosol_water_skip((void *) rxn_data);
        break;
    }

    int_size = (unsigned int) ((int *) float_data - (int *) rxn_start);
    int_total_size += int_size;
    int_sizes[i_rxn+1] = int_total_size;
    if(int_size>int_max_size) int_max_size=int_size;

    double_size = (unsigned int) ((double *) rxn_data - (double *) float_data);
    double_total_size += double_size;
    double_sizes[i_rxn+1] = double_total_size;
    if(double_size>double_max_size) double_max_size=double_size;

  }

  int_sizes[0]=0;
  double_sizes[0]=0;
  unsigned int rxn_int_size=n_rxn*int_max_size;
  unsigned int rxn_double_size=n_rxn*double_max_size;

  //Allocate int and double rxn data separately
  int_pointer = (int *) malloc(rxn_int_size * sizeof(int));
  memset(int_pointer, -1, rxn_int_size * sizeof(int));
  double_pointer = (double*)calloc(rxn_double_size, sizeof(double));
  cudaMalloc((void **) &int_pointer_gpu, rxn_int_size * sizeof(int));
  cudaMalloc((void **) &double_pointer_gpu, rxn_double_size * sizeof(double));

  //Copy into gpu rxn data
  //Rxn matrix is reversed to improve memory access on gpu
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

    int_size = int_sizes[i_rxn+1] - int_sizes[i_rxn];
    double_size = double_sizes[i_rxn+1] - double_sizes[i_rxn];

    for (int j = 0; j < int_size; j++)
      int_pointer[n_rxn*j+i_rxn] = ((int *) rxn_param)[start_rxn_param[i_rxn] + j]; //Matrix order [int_size][n_rxn]

    for (int j = 0; j < double_size; j++) {
      double *float_data = (double *) &(((int *) rxn_param)[start_rxn_param[i_rxn] + int_size]);
      double_pointer[n_rxn*j+i_rxn] = float_data[j];//Matrix order [int_size][n_rxn]
    }

  }

  //TODO: Avoid zeros or improve zeros access:
  // Quick sort to reorganize rows starting with low number of zeros to a lot of zeros in row

  #ifdef PMC_DEBUG
    printf(" Zeros_added_int_%: %f ", ((double) (n_rxn*int_max_size))/(n_rxn*int_max_size-int_total_size));
    printf(" Zeros_added_double_%: %f\n ", ((double) (n_rxn*double_max_size))/(n_rxn*double_max_size-double_total_size));
  #endif

  //Save data to GPU
  HANDLE_ERROR(cudaMemcpy(int_pointer_gpu, int_pointer, rxn_int_size*sizeof(int), cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(double_pointer_gpu, double_pointer, rxn_double_size*sizeof(double), cudaMemcpyHostToDevice));

  //Set flag to zero to avoid repetitions, since this data is fixed all time
  solver_set_gpu_sizes = 0;

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
          double time_step, int deriv_length_cell, int state_size_cell, int n_rxn,
          int n_cells, int *int_pointer, double *double_pointer,
          double *rate_constants_init) //Interface CPU/GPU
{
  //Get thread id
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  //Maximum number of threads to compute all reactions
  if(index < n_rxn*n_cells){

    //Thread index for deriv and state,
    // till we don't finish all reactions of a cell, we stay on same index
    int cell=index/n_rxn;

    //Get indices of each reaction
    int *int_data = (int *) &(((int *) int_pointer)[index%n_rxn]); //Same indices for each cell
    double *float_data = (double *) &(((double *) double_pointer)[index%n_rxn]);
    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *deriv_data = &( deriv_init[deriv_length_cell*cell]);
    double *state = &( state_init[state_size_cell*cell]);

    //Get indices for rates
    double *rate_constants = &( rate_constants_init[index]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(rate_constants,
        //        state, deriv_data, (void *) rxn_data, float_data, time_step, n_rxn);
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
        //rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(rate_constants,
        //        state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
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
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_deriv_contrib(rate_constants,
                state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(rate_constants,
                state, deriv_data, (void *) rxn_data, float_data, time_step,n_rxn);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(rate_constants,
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
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(rate_constants,
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
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];
  int n_rxn_threads = rxn_data[0]*n_cells; //Reaction group per number of repetitions/cells
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);
  double *state = model_data->state;
  double *rate_constants = model_data->rate_constants;

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

  HANDLE_ERROR(cudaMemset(deriv_gpu_data, 0, deriv_size));

  solveDerivative << < (n_blocks), max_n_gpu_thread >> >
    (state_gpu, deriv_gpu_data, time_step, model_data->n_dep_var, model_data->n_state_var,
    n_rxn, n_cells, int_pointer_gpu, double_pointer_gpu, rate_constants_gpu);

  cudaDeviceSynchronize();// Secure cuda synchronization

  //Use pinned memory for few values
  if (few_data){
    HANDLE_ERROR(cudaMemcpy(deriv_cpu, deriv_gpu_data, deriv_size, cudaMemcpyDeviceToHost));
    memcpy(deriv_data, deriv_cpu, deriv_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(deriv_data, deriv_gpu_data, deriv_size, cudaMemcpyDeviceToHost));
  }

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

__global__ void solveJacobian(double *state_init, double *jac_init,
        double time_step, int jac_length_cell, int state_size_cell, int n_rxn,
        int n_cells, int *int_pointer, double *double_pointer,
        double *rate_constants_init) //Interface CPU/GPU
{
  //Get thread id
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  //Maximum number of threads to compute all reactions
  if(index < n_rxn*n_cells){

    //Thread index for jac and state,
    // till we don't finish all reactions of a cell, we stay on same index
    int cell=index/n_rxn;

    //Get indices of each reaction
    int *int_data = (int *) &(((int *) int_pointer)[index%n_rxn]); //Same indices for each cell
    double *float_data = (double *) &(((double *) double_pointer)[index%n_rxn]);
    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1*n_rxn]);

    //Get indices for concentrations
    double *jac_data = &( jac_init[jac_length_cell*cell]);
    double *state = &( state_init[state_size_cell*cell]);

    //Get indices for rates
    double *rate_constants = &( rate_constants_init[index]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_jac_contrib(rate_constants,
        //        state, jac_data, (void *) rxn_data, float_data, time_step, n_rxn);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_jac_contrib(rate_constants,
                  state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
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
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_jac_contrib(rate_constants,
                state, jac_data, (void *) rxn_data, float_data, time_step,n_rxn);
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
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_jac_contrib(rate_constants,
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
  int n_cells = model_data->n_cells;
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];
  int n_rxn_threads = rxn_data[0]*n_cells; //Reaction group per number of repetitions/cells
  int n_blocks = ((n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread);
  double *jac_data = SM_DATA_S(jac);
  double *state = model_data->state;
  double *rate_constants = model_data->rate_constants;

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
    (state_gpu, jac_gpu_data, time_step, model_data->n_jac_elem, model_data->n_state_var,
    n_rxn, n_cells, int_pointer_gpu, double_pointer_gpu, rate_constants_gpu);

  cudaDeviceSynchronize();// Secure cuda synchronization

  //Use pinned memory for few values
  if (few_data){
    HANDLE_ERROR(cudaMemcpy(jac_cpu, jac_gpu_data, jac_size, cudaMemcpyDeviceToHost));
    memcpy(jac_data, jac_cpu, jac_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(jac_data, jac_gpu_data, jac_size, cudaMemcpyDeviceToHost));
  }

}



/** \brief Free GPU data structures
 */
void free_gpu_cu() {

  free(start_rxn_param);

  HANDLE_ERROR(cudaFree(int_pointer_gpu));
  HANDLE_ERROR(cudaFree(double_pointer_gpu));
  HANDLE_ERROR(cudaFree(deriv_gpu_data));
  HANDLE_ERROR(cudaFree(jac_gpu_data));
  HANDLE_ERROR(cudaFree(dev_start_rxn_param));

  if(few_data){
  }
  else{
    HANDLE_ERROR(cudaFree(state_gpu));
    HANDLE_ERROR(cudaFree(rate_constants_gpu));
  }

}

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

/* Functions on development */

/** \brief GPU function: Update reaction data for new environmental state. Currently in development:
 **  Not time improvement for now.
 *
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
__global__ void updateEnvRxnBlock(double *rate_constants_init, int n_rxn, int n_rxn_threads,
                                  int n_cells_gpu, int *int_pointer, double *env_init, double *double_pointer){

  int index = blockIdx.x * blockDim.x + threadIdx.x;

  if (index < n_rxn_threads) {

    int env_size_cell = PMC_NUM_ENV_PARAM_; //Temperature and pressure
    int cell = index / n_rxn;

    int *int_data = (int *) &(((int *) int_pointer)[index % n_rxn]);
    double *float_data = (double *) &(((double *) double_pointer)[index % n_rxn]);

    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1 * n_rxn]);

    double *env = &( env_init[env_size_cell * cell]);
    double *rate_constants = &( rate_constants_init[index]);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_TROE :
        rxn_gpu_troe_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_update_env_state(rate_constants, n_rxn, float_data, env, rxn_data);
        break;
    }
    __syncthreads();
  }
}

/** \brief Update reaction data on GPU for new environmental state. Currently in development:
 **  Not time improvement for now.
 *
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
void rxn_update_env_state_gpu(ModelData *model_data, double *env){

  int n_cells = model_data->n_cells;
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];
  int n_rxn_threads = rxn_data[0]*n_cells;

  //Note: Pinned memory or function parameter don't achieve an improvement

  HANDLE_ERROR(cudaMemcpy(env_gpu, env, env_size, cudaMemcpyHostToDevice));

  updateEnvRxnBlock << < (n_rxn_threads + max_n_gpu_thread - 1) / max_n_gpu_thread, max_n_gpu_thread >> >
                                                                                    (rate_constants_gpu, n_rxn, n_rxn_threads, n_cells, int_pointer_gpu, env_gpu, double_pointer_gpu);
  cudaDeviceSynchronize();

  HANDLE_ERROR(cudaMemcpy(model_data->rate_constants, rate_constants_gpu, rate_constants_size, cudaMemcpyDeviceToHost));

}


}
