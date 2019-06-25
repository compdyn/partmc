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

//TtODO: check if a double is necessary, or can we use a float, since cvode
// check a tolerance error to converge and maybe this error only reach float
// (and also the test comparison between ebi and pmc uses a tolerance, maybe float)

ModelDatagpu *mdgpu;
ModelDatagpu *mdcpu;
double *derivgpu_data;
double *deriv_cpu;
size_t deriv_size; //unsigned int
unsigned int countergpu = 0;

size_t state_size;
bool solver_set_gpu_sizes = 1;
int *short_pointer;
double *double_pointer;
int *short_pointer_gpu;
double *double_pointer_gpu;
unsigned int *start_rxn_param;
unsigned int *dev_start_rxn_param;
unsigned int *int_sizes;
unsigned int *int_sizes_gpu;
unsigned int *double_sizes;
unsigned int *double_sizes_gpu;
unsigned int int_max_size = 0;
unsigned int double_max_size = 0;
double *state_gpu;

static void HandleError(cudaError_t err,
                        const char *file,
                        int line) {
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err),
           file, line);
    exit(EXIT_FAILURE);
  }
}

void solver_new_gpu_cu(SolverDatagpu *sd, int n_dep_var,
                       int n_state_var, int *var_type, int n_rxn,
                       int n_rxn_int_param, int n_rxn_float_param) {


  /*HANDLE_ERROR(cudaHostAlloc((void **) &md,//no work well
                             sizeof(ModelDatagpu),
                             cudaHostAllocWriteCombined |
                             cudaHostAllocMapped
  ));
  HANDLE_ERROR(cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void **) &(mdgpu),
          (void *) md,
          0));

  HANDLE_ERROR(cudaHostAlloc((void **) &md->state,
                             state_size,
                             cudaHostAllocWriteCombined |
                             cudaHostAllocMapped
  ));
  HANDLE_ERROR(cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void **) &(mdgpu->state),
          (void *) md->state,
          0));*/

  ModelDatagpu *md = &sd->model_data;

  //Sizes
  size_t start_size = (n_rxn+1) * sizeof(unsigned int);
  state_size = n_state_var * sizeof(double);
  deriv_size = n_dep_var * sizeof(sd->y);
  printf("n_rxn: %d" , n_rxn);
  //printf("n_rxn_float_param: %d ", n_rxn_float_param);
  //printf("n_rxn_int_param: %d ", n_rxn_int_param);
  //printf("n_state_var: %d" ,n_state_var);

  //Create started indexes of arrays
  start_rxn_param = (unsigned int *) malloc(start_size);
  int_sizes = (unsigned int *) malloc(start_size);
  double_sizes = (unsigned int *) malloc(start_size);

  //GPU allocation
  cudaMalloc((void **) &dev_start_rxn_param, start_size);
  cudaMalloc((void **) &int_sizes_gpu, start_size);
  cudaMalloc((void **) &double_sizes_gpu, start_size);



  //realtype *deriv_data = N_VGetArrayPointer(sd->y);
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterPortable));//pinned


  cudaMalloc((void **) &derivgpu_data, deriv_size);
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterMapped));//pinned, not work properly
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterDefault));
  //HANDLE_ERROR(cudaHostGetDevicePointer((void**) &(derivgpu_data), (void*)deriv_data, 0));

  cudaMallocHost((void**)&deriv_cpu, deriv_size);//pinned

  //cudaMalloc((void **) &state_gpu, state_size);
  //cudaMalloc((void **) &mdgpu, sizeof(*mdgpu));

  //HANDLE_ERROR(cudaHostGetDevicePointer((void**) &(derivgpu_data), (void*)deriv_data, 0));

  if (deriv_size > MAX_SHARED_MEMORY_BLOCK_DOUBLE)
#ifndef FAILURE_DETAIL
    printf("\nMore solver variables(deriv): %d than maximum shared memory: %d",
           deriv_size, MAX_SHARED_MEMORY_BLOCK_DOUBLE);
#endif

}

void solver_update_state_gpu(ModelDatagpu *md) {//HANDLE_ERROR(cudaMemcpy(mdgpu->state, md->state, state_size*sizeof(int), cudaMemcpyHostToDevice));
}

void get_rxn_pointers (int *short_pointer0, double *double_pointer0,
      unsigned int *int_sizes0, unsigned int *double_sizes0){

  short_pointer0 = short_pointer;
  double_pointer0 = double_pointer;
  int_sizes0 = int_sizes;
  double_sizes0 = double_sizes;
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

    for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

      //Reaction distances between pointers rows
      start_rxn_param[i_rxn] = (unsigned int) ((int *) rxn_data - (int *) rxn_param);

      int *rxn_start = rxn_data;

      // Get the reaction type
      int rxn_type = *(rxn_data++);

      //TTODO: Join functions to retrieve float and int sizes
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

    //Allocate int and double rxn data separately
    short_pointer = (int *) malloc(n_rxn*int_max_size * sizeof(int));
    double_pointer = (double *) malloc(n_rxn*double_max_size * sizeof(double));
    cudaMalloc((void **) &short_pointer_gpu, n_rxn*int_max_size  * sizeof(int));
    cudaMalloc((void **) &double_pointer_gpu, n_rxn*double_max_size  * sizeof(double));

    //TODO: Transporse rxn matrix to better memory acces
    for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
      //TTODO: use short_pointer in reactions

      int_size = int_sizes[i_rxn+1] - int_sizes[i_rxn];
      double_size = double_sizes[i_rxn+1] - double_sizes[i_rxn];

      for (int j = 0; j < int_size; j++)
        short_pointer[int_sizes[i_rxn] + j] = ((int *) rxn_param)[start_rxn_param[i_rxn] + j];
      for (int j = int_size; j < (int_max_size-int_size); j++)
        short_pointer[int_sizes[i_rxn] + j] = -1;

      for (int j = 0; j < double_size; j++) {
        double *float_data = (double *) &(((int *) rxn_param)[start_rxn_param[i_rxn] + int_size]);
        double_pointer[double_sizes[i_rxn] + j] = float_data[j];
      }
      for (int j = double_size; j < (double_max_size-double_size); j++)
        double_pointer[double_sizes[i_rxn] + j] = 0;

      printf(" int_offset: %d ", (int_max_size-int_size));
      printf(" double_offset: %d  \n", (double_max_size-double_size));

      //Ttodo: Quick sort to reorganize rows starting with low number of zeros to a lot of zeros in row
    }

    HANDLE_ERROR(cudaMemcpy(short_pointer_gpu, short_pointer, int_total_size*sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(double_pointer_gpu, double_pointer, double_total_size*sizeof(double), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(int_sizes_gpu, int_sizes, n_rxn*sizeof(unsigned int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(double_sizes_gpu, double_sizes, n_rxn*sizeof(unsigned int), cudaMemcpyHostToDevice));

    //TTODO: Create global size arrays (int and float) with each size for reaction instead the defines on reaction files
    //TODO: Some rxn values changes on monarch each iteration(temperature/pressure-update function) by update functions
    solver_set_gpu_sizes = 0;
  }
}

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *deriv,
                              double time_step, int deriv_length, int n_rxn,
                              int *short_pointer, double *double_pointer,
                              unsigned int *int_sizes, unsigned int *double_sizes) //Interface CPU/GPU
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ double deriv_data[MAX_SHARED_MEMORY_BLOCK_DOUBLE];

  if (index < deriv_length){ //This produces seg.fault for some large values seems
    deriv_data[index] = 0.0;
    deriv[index] = 0.0;
  }
  //TODO: Comment guillermo this initialization

  //for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x)
  //{
  //  deriv_data[index] = 0.0;
  //  deriv[index] = 0.0;
  //}

  __syncthreads();

  if (index < n_rxn) {

    int *int_data = (int *) &(((int *) short_pointer)[int_sizes[index]]);
    double *float_data = (double *) &(((double *) double_pointer)[double_sizes[index]]);

    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1]);

    //Notice till this, without executing calc_deriv switch, gpu version takes 2 secs more
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(
        //        model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        //atomicAdd(&(deriv[0]), 1);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(
                  model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(
                  model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(
        //        model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
        //        model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(
                model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(
               model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(
               model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
        break;

      //case RXN_ARRHENIUS :
      //  rxn_gpu_arrhenius_calc_deriv_contrib(
      //          model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      //break;

//atomicAdd(&(deriv[index]), float_data[0]);
    }
  }
  __syncthreads();

  if (index < deriv_length){//seems better but i think dont work for large values
    deriv[index] = deriv_data[index];
  }
  //for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x)
    //atomicAdd(&(deriv[i_spec]), deriv_data[i_spec]);
}

void rxn_calc_deriv_gpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step) {

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];

  //TODO: Only 1 time call reaction_int_size,get int sizes, assign values to rxn_datagpu with only doubles (no array of ints and doubles)

  /*if(countergpu==29) {
    for (int i_rxn=1; i_rxn<n_rxn; i_rxn++) {
      int * int_data = (int*) &(((int*)  rxn_data)[start_rxn_param[ i_rxn ]] );
      int rxn_type = int_data[0];
      int * rxn_data = (int*)&(int_data[1]);
     // int int_data = *((int *) &(((int *) short_pointer)[int_sizes[i_rxn - 1]]));

      if(rxn_type==RXN_AQUEOUS_EQUILIBRIUM){
        //double *float_data =  (double*) &(((int*)rxn_data)[
        //        2+(((int*)rxn_data)[0]+2)*(((int*)rxn_data)[0]+((int*)rxn_data)[1]) ]);
        //double * float_data2 = (double *) &(((double *) double_pointer)[double_sizes[i_rxn-1]]);
        //printf(" type1: %f  ", float_data[6]);
        //printf(" type2: %f\n", float_data2[6]);
        int int_data = *((int *) &(((int *) short_pointer)[int_sizes[i_rxn - 1]]));
        printf(" type1: %d  ", int_data);
      }
    }
  }
  countergpu++;*/

  mdgpu = model_data; //This synchronyze state with state_gpu

  //HANDLE_ERROR(cudaMemcpy(state_gpu, model_data->state, state_size, cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy(&(mdgpu->state), &state_gpu, sizeof(mdgpu), cudaMemcpyHostToDevice));
  //TODO:speedup with this memcpy is 0.776 infront mdgpu=model_data

  //TTODO:Seems for big sizes only works with if thread==2 deriv=0, despite it gives different results

solveRxnBlock << < (n_rxn + MAX_N_GPU_THREAD - 1) / MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >> >
    (mdgpu, derivgpu_data, time_step, NV_LENGTH_S(deriv), //deriv_size/sizeof(double)
            n_rxn, short_pointer_gpu, double_pointer_gpu, int_sizes_gpu, double_sizes_gpu);
  cudaDeviceSynchronize();//retrieve errors (But don't retrieve anything for me)
  //HANDLE_ERROR(cudaMemcpy(deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));//0.5secs
  HANDLE_ERROR(cudaMemcpy(deriv_cpu, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));//0.29secs

  memcpy(deriv_data, deriv_cpu, deriv_size); //This is so fast(0.01secs)

}

void solveRxncpu(ModelDatagpu *model_data, double *deriv_data,
                   double time_step, int *int_data, double *float_data, int deriv_length)
{

  int rxn_type = int_data[0];
  int *rxn_data = (int *) &(int_data[1]);

  switch (rxn_type) {
    case RXN_AQUEOUS_EQUILIBRIUM :
      rxn_cpu_aqueous_equilibrium_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_ARRHENIUS :
      rxn_cpu_arrhenius_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_CMAQ_H2O2 :
      rxn_cpu_CMAQ_H2O2_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_CMAQ_OH_HNO3 :
      rxn_cpu_CMAQ_OH_HNO3_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_CONDENSED_PHASE_ARRHENIUS :
      rxn_cpu_condensed_phase_arrhenius_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_EMISSION :
      rxn_cpu_emission_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_FIRST_ORDER_LOSS :
      rxn_cpu_first_order_loss_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_HL_PHASE_TRANSFER :
//rxn_cpu_HL_phase_transfer_calc_deriv_contrib(
//        model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_PDFITE_ACTIVITY :
      rxn_cpu_PDFiTE_activity_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_PHOTOLYSIS :
      rxn_cpu_photolysis_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_SIMPOL_PHASE_TRANSFER :
//rxn_cpu_SIMPOL_phase_transfer_calc_deriv_contrib(
//        model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_TROE :
      rxn_cpu_troe_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_WET_DEPOSITION :
      rxn_cpu_wet_deposition_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;
    case RXN_ZSR_AEROSOL_WATER :
      rxn_cpu_ZSR_aerosol_water_calc_deriv_contrib(
              model_data, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length);
      break;

  }

}

void rxn_calc_deriv_cpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step)
{

  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];

  //Case i_rxn=0
  solveRxncpu(model_data, deriv_data, time_step, short_pointer, double_pointer, NV_LENGTH_S(deriv));

  for (int i_rxn = 1; i_rxn < n_rxn; i_rxn++) {

    int *int_data = (int *) &(((int *) short_pointer)[int_sizes[i_rxn - 1]]);
    double *float_data = (double *) &(((double *) double_pointer)[double_sizes[i_rxn - 1]]);

    solveRxncpu(model_data, deriv_data, time_step, int_data, float_data, NV_LENGTH_S(deriv));
  }
}

void free_gpu_cu() {

  //HANDLE_ERROR( cudaFreeHost( derivgpu_data ) );
  //HANDLE_ERROR( cudaFreeHost( mdgpu ) );

  cudaFree( short_pointer_gpu );
  cudaFree( double_pointer_gpu );
  cudaFree( int_sizes_gpu );
  cudaFree( double_sizes_gpu );

  //cudaFree( derivgpu_data ); //TODO: ESTO PETA 
  //cudaFree( mdgpu );
  //HANDLE_ERROR( cudaFreeHost( deriv ) );

}

}
