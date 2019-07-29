//Test
extern "C" {
#include "phlex_gpu_solver.h"
#include "rxns_gpu.h"
//#include "phlex_solver.h"
//}

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

#define CHEM_SPEC_VARIABLE 1

//Gpu hardware info
int max_shared_memory_block_double;
int max_n_gpu_thread;
int max_n_gpu_blocks;

//General data
ModelDatagpu *mdgpu;
double *derivgpu_data;
double *deriv_cpu;
size_t deriv_size;
unsigned int countergpu = 0;

size_t state_size;
size_t env_size;
size_t rate_constants_size;
bool solver_set_gpu_sizes = 1;
int *int_pointer;
double *double_pointer;
int *int_pointer_gpu;
double *double_pointer_gpu;
unsigned int *start_rxn_param;
unsigned int *dev_start_rxn_param;
unsigned int int_max_size = 0;
unsigned int double_max_size = 0;
double *state_gpu;
double *state_cpu;
double *env_gpu;
double *rate_constants_gpu;
double *rate_constants_cpu;
bool few_data = 0;

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

void solver_new_gpu_cu(int n_dep_var,
                       int n_state_var, int *var_type, int n_rxn,
                       int n_rxn_int_param, int n_rxn_float_param,
                       int n_cells) { //Ttodo: not necessary pass this parameters, there are on md

  //Sizes
  size_t start_size = (n_rxn+1) * sizeof(unsigned int);
  state_size = n_state_var*n_cells * sizeof(double); 
  deriv_size = n_dep_var*n_cells * sizeof(double);
  env_size = 2*n_cells * sizeof(double); //Temp and pressure
  rate_constants_size = n_rxn * n_cells * sizeof(double);

  //Create started indexes of arrays
  start_rxn_param = (unsigned int *) malloc(start_size);

#ifndef PMC_DEBUG
  printf("n_rxn: %d " , n_rxn);
  printf("n_state_var: %d" ,n_state_var*n_cells);
  printf("n_dep_var: %d" ,n_dep_var*n_cells);
#endif

  //Detect if we are working with few data values
  if (n_dep_var*n_cells < DATA_SIZE_LIMIT_OPT)
    few_data = 1;

  //GPU allocation
  cudaMalloc((void **) &dev_start_rxn_param, start_size);

  //realtype *deriv_data = N_VGetArrayPointer(sd->y);
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterPortable));//pinned

  cudaMalloc((void **) &derivgpu_data, deriv_size);
  cudaMalloc((void **) &state_gpu, state_size);
  cudaMalloc((void **) &env_gpu, env_size);
  cudaMalloc((void **) &rate_constants_gpu, rate_constants_size);

  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterMapped));//pinned, not work properly
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterDefault));
  //HANDLE_ERROR(cudaHostGetDevicePointer((void**) &(derivgpu_data), (void*)deriv_data, 0));

  if(few_data){
    cudaMallocHost((void**)&rate_constants_cpu, rate_constants_size);
    cudaMallocHost((void**)&deriv_cpu, deriv_size);//pinned
  }

  //deriv_cpu = (double *) malloc(deriv_size);
  //cudaMallocHost((void**)&state_cpu, state_size);//pinned

  //cudaMalloc((void **) &mdgpu, sizeof(*mdgpu));

  //HANDLE_ERROR(cudaHostGetDevicePointer((void**) &(derivgpu_data), (void*)deriv_data, 0));

  int nDevices;

  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
#ifdef PMC_DEBUG
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
    printf("  sharedMemPerBlock: %zu\n", prop.sharedMemPerBlock);
    printf("  multiProcessorCount: %d\n", prop.multiProcessorCount);
#endif
  }

  //Set Gpu to work: we have 4 gpu available on power9. as default, it should be assign to gpu 0
  int device = 0;
  //cudaSetDevice( device ); //TODO: Selecting different gpu raise illegal memory access
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);

  //max_shared_memory_block_double = prop.sharedMemPerBlock/sizeof(double); //TODO: Fix static declaration of shared memory
  max_n_gpu_thread = prop.maxThreadsPerBlock/2;//TODO: fix bug too many resources without the /2
  max_n_gpu_blocks = prop.maxGridSize[1]; //test if is 3d or 1d the maximum

  int n_blocks = (n_rxn + max_n_gpu_thread - 1) / max_n_gpu_thread;

// Control exceeding limits

  if( n_blocks > max_n_gpu_blocks){
    printf("\nWarning: More blocks assigned: %d than maximum block numbers: %d",
           n_blocks, max_n_gpu_blocks);
  }

  if (max_n_gpu_thread > MAX_SHARED_MEMORY_BLOCK_DOUBLE)
    printf("\nWarning: More threads assigned: %d than maximum shared memory: %d",
           max_n_gpu_thread, MAX_SHARED_MEMORY_BLOCK_DOUBLE);

}

void solver_update_state_gpu(ModelDatagpu *md) {//HANDLE_ERROR(cudaMemcpy(mdgpu->state, md->state, state_size*sizeof(int), cudaMemcpyHostToDevice));
}

void solver_set_data_gpu(ModelDatagpu *model_data) {
  //Get rxn sizes
  if (solver_set_gpu_sizes) {

    //HANDLE_ERROR(cudaHostRegister(model_data->state, state_size, cudaHostRegisterPortable)); pinned

    // Get the number of reactions
    int *rxn_data = (int *) (model_data->rxn_data);
    int n_rxn = *(rxn_data++);
    void *rxn_param = (void *) rxn_data;
    int *float_data = (int *) rxn_data;
    unsigned int int_size = 0;
    unsigned int int_total_size = 0;
    unsigned int double_size = 0;
    unsigned int double_total_size = 0;

    size_t start_size = (n_rxn+1) * sizeof(unsigned int);
    unsigned int int_sizes[start_size];
    unsigned int double_sizes[start_size];

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

    //Ttodo: best option is put sizes array of rxn inside rxn matrix taking advantage of read memory, or avoid zeros
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
    //Rxn matrix is rotated to improve memory access on gpu
    for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

      int_size = int_sizes[i_rxn+1] - int_sizes[i_rxn];
      double_size = double_sizes[i_rxn+1] - double_sizes[i_rxn];

      for (int j = 0; j < int_size; j++)
        int_pointer[n_rxn*j+i_rxn] = ((int *) rxn_param)[start_rxn_param[i_rxn] + j]; //[int_size][n_rxn]

      for (int j = 0; j < double_size; j++) {
        double *float_data = (double *) &(((int *) rxn_param)[start_rxn_param[i_rxn] + int_size]);
        double_pointer[n_rxn*j+i_rxn] = float_data[j];//[int_size][n_rxn]
      }

    }
    //TODO: Quick sort to reorganize rows starting with low number of zeros to a lot of zeros in row

#ifndef PMC_DEBUG
    printf(" Zeros_added_int_%: %f ", ((double) (n_rxn*int_max_size))/(n_rxn*int_max_size-int_total_size));
    printf(" Zeros_added_double_%: %f\n ", ((double) (n_rxn*double_max_size))/(n_rxn*double_max_size-double_total_size));
#endif

    HANDLE_ERROR(cudaMemcpy(int_pointer_gpu, int_pointer, rxn_int_size*sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(double_pointer_gpu, double_pointer, rxn_double_size*sizeof(double), cudaMemcpyHostToDevice));

    solver_set_gpu_sizes = 0;
  }
}

__global__ void updateEnvRxnBlock(double *rate_constants_init, int n_rxn, int n_rxn_total_threads,
       int n_cells_gpu, int *int_pointer, double *env_init, double *double_pointer){

  int index = blockIdx.x * blockDim.x + threadIdx.x;

  //I think shared is not useful since we are accessing rate_constant only 1 time thread
  //__shared__ double rate_constants[MAX_SHARED_MEMORY_BLOCK_DOUBLE];

  if (index < n_rxn_total_threads) {

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

void rxn_update_env_state_gpu(ModelDatagpu *model_data, double *env){

  int n_cells = model_data->n_cells;
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];
  int n_rxn_total_threads = rxn_data[0]*n_cells;

  if (few_data){
    //env_gpu=env;//Not clear improvement for the data size tested
    HANDLE_ERROR(cudaMemcpy(env_gpu, env, env_size, cudaMemcpyHostToDevice));
  }
  else{
    HANDLE_ERROR(cudaMemcpy(env_gpu, env, env_size, cudaMemcpyHostToDevice));
  }

  updateEnvRxnBlock << < (n_rxn_total_threads + max_n_gpu_thread - 1) / max_n_gpu_thread, max_n_gpu_thread >> >
  (rate_constants_gpu, n_rxn, n_rxn_total_threads, n_cells, int_pointer_gpu, env_gpu, double_pointer_gpu);
  cudaDeviceSynchronize();

  //This takes so much time, let alive for the moment only for testing purposes!!
  /*
  if (few_data){
    //HANDLE_ERROR(cudaMemcpy(rate_constants_cpu, rate_constants_gpu, rate_constants_size, cudaMemcpyDeviceToHost));
    //memcpy(model_data->rate_constants, rate_constants_cpu, rate_constants_size);
    //Seems pinned memory is not giving any improve
    HANDLE_ERROR(cudaMemcpy(model_data->rate_constants, rate_constants_gpu, rate_constants_size, cudaMemcpyDeviceToHost));
  }
  else{
    HANDLE_ERROR(cudaMemcpy(model_data->rate_constants, rate_constants_gpu, rate_constants_size, cudaMemcpyDeviceToHost));
  }*/

}

//TODO: FIX THIS BUG (GUILLERMO OR SOMEONE) (this works, but if move this into another file with maxthreads=1024, it crash
__device__ void rxn_gpu_tmp_arrhenius2(
      ModelDatagpu *model_data,
      double *deriv, int *rxn_data, double * double_pointer_gpu,
      double time_step, int n_rxn)
{

  double *state = model_data->state;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the reaction rate
  double rate = float_data[6*n_rxn];
  for (int i_spec=0; i_spec<int_data[0]; i_spec++) rate *= state[int_data[(2 + i_spec)*n_rxn]-1];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<int_data[0]; i_spec++, i_dep_var++) {
      if (int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn] < 0) continue;
      //deriv[DERIV_ID_(i_dep_var)] -= rate;
      atomicAdd(&(deriv[int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn]]),-rate);
    }
    for (int i_spec=0; i_spec<int_data[1*n_rxn]; i_spec++, i_dep_var++) {
      if (int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn] < 0) continue;

      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-rate*float_data[(7 + i_spec)*n_rxn]*time_step <=
          state[int_data[(2 + int_data[0] + i_spec)*n_rxn]-1]) {
        //deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
        atomicAdd(&(deriv[int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn]]),
                  rate*float_data[(7 + i_spec)*n_rxn]);
      }
    }
  }
}

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *state_init, double *deriv,
          double time_step, int deriv_length_cell, int state_size_cell, int n_rxn,
          int n_cells_gpu, int *int_pointer, double *double_pointer,
          double *rate_constants_init, int n_blocks) //Interface CPU/GPU
{
   //rxn_gpu_tmp_arrhenius( //To test the bug with large threads number
          //model_data, deriv, int_pointer, double_pointer, time_step, n_rxn
     //     );

  //Get thread id
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ double deriv_init[MAX_SHARED_MEMORY_BLOCK_DOUBLE];

  int deriv_length = deriv_length_cell * n_cells_gpu;
  //Repeat shared memory usage if is less than deriv_length
  int repeat_shared = (deriv_length/MAX_SHARED_MEMORY_BLOCK_DOUBLE)+1;
  int n_rxn_total_threads = n_rxn * n_cells_gpu / repeat_shared;
  //extern __shared__ double deriv_shr[];
  //double *deriv_init = &deriv_shr[0];

  if (tid < n_rxn_total_threads) {

    //TODO: for big repeat_shared results changes a bit, investigate this
    for(int i=0;i<repeat_shared;i++){

      //Update index to access data
      int index = tid + i*n_rxn_total_threads;

      //Cell corresponding acces deriv and state data
      int cell=index/n_rxn;

      //Index to access is less than shared deriv size
      if (deriv_length_cell*cell < MAX_SHARED_MEMORY_BLOCK_DOUBLE){

        if (tid < MAX_SHARED_MEMORY_BLOCK_DOUBLE)
          deriv_init[tid] = 0.0;

        __syncthreads();

        int *int_data = (int *) &(((int *) int_pointer)[index%n_rxn]);
        double *float_data = (double *) &(((double *) double_pointer)[index%n_rxn]);
        int rxn_type = int_data[0];
        int *rxn_data = (int *) &(int_data[1*n_rxn]);

        double *deriv_data = &( deriv_init[deriv_length_cell*cell]);
        double *state = &( state_init[state_size_cell*cell]);
        double *rate_constants = &( rate_constants_init[index]);

        switch (rxn_type) {
          case RXN_AQUEOUS_EQUILIBRIUM :
            //rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(rate_constants,
            //        state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_ARRHENIUS :
            rxn_gpu_arrhenius_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_CMAQ_H2O2 :
            rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_CMAQ_OH_HNO3 :
            rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_CONDENSED_PHASE_ARRHENIUS :
            rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_EMISSION :
            rxn_gpu_emission_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_FIRST_ORDER_LOSS :
            rxn_gpu_first_order_loss_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_HL_PHASE_TRANSFER :
            //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(rate_constants,
            //        state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_PDFITE_ACTIVITY :
            rxn_gpu_PDFiTE_activity_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_PHOTOLYSIS :
            rxn_gpu_photolysis_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_SIMPOL_PHASE_TRANSFER :
            //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(rate_constants,
            //        state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_TROE :
            rxn_gpu_troe_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_WET_DEPOSITION :
            rxn_gpu_wet_deposition_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;
          case RXN_ZSR_AEROSOL_WATER :
            rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(rate_constants,
                    state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
            break;

            //case RXN_ARRHENIUS :
              //rxn_gpu_arrhenius_calc_deriv_contrib(rate_constants,
              //        state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);

            //rxn_gpu_tmp_arrhenius(
            //     model_data, deriv_data, rxn_data, float_data, time_step, n_rxn);

            //break;
        }
        __syncthreads();

        //Save in the big global memory the shared memory data calculated
          if (tid < MAX_SHARED_MEMORY_BLOCK_DOUBLE)
            deriv[index] = deriv_init[tid];
      }
    }
  }

  //if (threadIdx.x < deriv_length)
  //  deriv[index] = deriv_init[index];

}

void rxn_calc_deriv_gpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step) {

  // Get a pointer to the derivative data
  int n_cells = model_data->n_cells;
  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];
  int n_rxn_total_threads = rxn_data[0]*n_cells;
  double *state = model_data->state;
  double *rate_constants = model_data->rate_constants;

  /*if(countergpu==29) {
    for (int i_rxn=1; i_rxn<n_rxn; i_rxn++) {
      if(rxn_type==RXN_AQUEOUS_EQUILIBRIUM){
        printf(" type1: %d  ", int_data);
      }
    }
  }
  countergpu++;*/

  if (few_data){
    state_gpu= state;//Faster, use for few values
  }
  else{
    HANDLE_ERROR(cudaMemcpy(state_gpu, state, state_size, cudaMemcpyHostToDevice));//Slower, use for large values
  }

  rate_constants_gpu= rate_constants;

//(n_rxn*n_cells + max_n_gpu_thread - 1) / max_n_gpu_thread
  int n_blocks = ((n_rxn_total_threads + max_n_gpu_thread - 1) / max_n_gpu_thread)+1;
  solveRxnBlock << < (n_blocks-1), max_n_gpu_thread >> >
    (mdgpu, state_gpu, derivgpu_data, time_step, model_data->n_dep_var, model_data->n_state_var,
    n_rxn, n_cells, int_pointer_gpu, double_pointer_gpu, rate_constants_gpu, n_blocks);

  cudaDeviceSynchronize();//retrieve errors (But don't retrieve anything for me)

  if (few_data){
    HANDLE_ERROR(cudaMemcpy(deriv_cpu, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));
    memcpy(deriv_data, deriv_cpu, deriv_size);
  }
  else {
    HANDLE_ERROR(cudaMemcpy(deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));
  }

}

void free_gpu_cu() {

  //HANDLE_ERROR( cudaFreeHost( derivgpu_data ) );
  //HANDLE_ERROR( cudaFreeHost( mdgpu ) );

  HANDLE_ERROR(cudaFree( int_pointer_gpu ));
  HANDLE_ERROR(cudaFree( double_pointer_gpu ));
  HANDLE_ERROR(cudaFree(derivgpu_data));
  HANDLE_ERROR(cudaFree(dev_start_rxn_param));

  //HANDLE_ERROR(cudaFree(state_gpu)); //Invalid device pointer
  //HANDLE_ERROR(cudaFree(env_gpu));
  //HANDLE_ERROR(cudaFree(rate_constants_gpu));
  free(start_rxn_param);

  if(few_data){
  //  free(deriv_cpu); //TODO: fix segmentation fault
  }

  //free(int_pointer_gpu); //DO not, segmentation fault
  //free(double_pointer_gpu);//DO not, segmentation fault
  //free(deriv_cpu); DO not, segmentation fault
  //cudaFree( derivgpu_data ); //DO not, segmentation fault
  //cudaFree( mdgpu );
  //HANDLE_ERROR( cudaFreeHost( deriv ) );

}


//Todo: multiples gpus (be careful with setDevice):
/*
cudaSetDevice( 0 );
kernel<<<...>>>(...);
cudaMemcpyAsync(...);
cudaSetDevice( 1 );
kernel<<<...>>>(...);
 */

}
