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
size_t deriv_size;
unsigned int countergpu = 0;

size_t state_size;
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

//int num_cells = 1;

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

void solver_new_gpu_cu(SolverDatagpu *sd, int n_dep_var,
                       int n_state_var, int *var_type, int n_rxn,
                       int n_rxn_int_param, int n_rxn_float_param,
                       int num_cells_aux) { //Ttodo: not necessary pass this parameters, there are on md


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
  state_size = n_state_var * sizeof(double); //TODO: Adapt state to has the same size as deriv?
  deriv_size = n_dep_var * sizeof(sd->y);
  printf("n_rxn: %d " , n_rxn);
  //printf("n_rxn_float_param: %d ", n_rxn_float_param);
  //printf("n_rxn_int_param: %d ", n_rxn_int_param);
  printf("n_state_var: %d" ,n_state_var); //TODO: Fix n_state_var received is so big (47 when two cells)

  //Create started indexes of arrays
  start_rxn_param = (unsigned int *) malloc(start_size);

  //GPU allocation
  cudaMalloc((void **) &dev_start_rxn_param, start_size);

  //realtype *deriv_data = N_VGetArrayPointer(sd->y);
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterPortable));//pinned


  cudaMalloc((void **) &derivgpu_data, deriv_size);
  cudaMalloc((void **) &state_gpu, state_size);
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterMapped));//pinned, not work properly
  //HANDLE_ERROR(cudaHostRegister(deriv_data, deriv_size, cudaHostRegisterDefault));
  //HANDLE_ERROR(cudaHostGetDevicePointer((void**) &(derivgpu_data), (void*)deriv_data, 0));

  cudaMallocHost((void**)&deriv_cpu, deriv_size);//pinned
  //cudaMallocHost((void**)&state_cpu, state_size);//pinned

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

    //TODO: Ask mario if is better save sizes of rxn or fill with -1 and put an if inside for to avoid it
    //Ttodo: best option is put sizes array of rxn insdie rxn matrix taking advantage of read memory, or avoid zeros
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

    for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {

      int_size = int_sizes[i_rxn+1] - int_sizes[i_rxn];
      double_size = double_sizes[i_rxn+1] - double_sizes[i_rxn];

      for (int j = 0; j < int_size; j++)
        int_pointer[n_rxn*j+i_rxn] = ((int *) rxn_param)[start_rxn_param[i_rxn] + j]; //[int_size][n_rxn]

      for (int j = 0; j < double_size; j++) {
        double *float_data = (double *) &(((int *) rxn_param)[start_rxn_param[i_rxn] + int_size]);
        double_pointer[n_rxn*j+i_rxn] = float_data[j];//[int_size][n_rxn]
      }

      //printf(" int_offset: %d ", (int_max_size-int_size));
      //printf(" double_offset: %d  \n", (double_max_size-double_size));

      //Ttodo: Quick sort to reorganize rows starting with low number of zeros to a lot of zeros in row
    }

    printf(" Zeros_added_int_%: %f ", ((double) (n_rxn*int_max_size))/(n_rxn*int_max_size-int_total_size));
    printf(" Zeros_added_double_%: %f\n ", ((double) (n_rxn*double_max_size))/(n_rxn*double_max_size-double_total_size));

    HANDLE_ERROR(cudaMemcpy(int_pointer_gpu, int_pointer, rxn_int_size*sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(double_pointer_gpu, double_pointer, rxn_double_size*sizeof(double), cudaMemcpyHostToDevice));

    //TODO: Update Some rxn values, changes on monarch each iteration(temperature/pressure-update function) by update functions
    solver_set_gpu_sizes = 0;
  }
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

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *state, double *deriv,
          double time_step, int deriv_length, int state_size, int n_rxn_total_cells,
          int num_cells_gpu, int *int_pointer, double *double_pointer, unsigned int int_max_size,
          unsigned int double_max_size) //Interface CPU/GPU
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

   //rxn_gpu_tmp_arrhenius(
          //model_data, deriv, int_pointer, double_pointer, time_step, n_rxn
     //     );

  __shared__ double deriv_data2[MAX_SHARED_MEMORY_BLOCK_DOUBLE];

  if (threadIdx.x < deriv_length){ //This produces seg.fault for some large values seems
    deriv_data2[index] = 0.0;
  }

  //if(index==2)//!!DON'T DELETE, it gives bad result but allows take measurements for large input values
  //for (int i=0; i<deriv_length; i++)
  //deriv_data[i]=0.0;

  //if (threadIdx.x < deriv_length)//dont work?
    //for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x) {
      //deriv_data[index] = 0.0; }
  __syncthreads();

  if (index < n_rxn_total_cells) {

    int state_size_cell = state_size/num_cells_gpu;

    int deriv_length_cell = deriv_length/num_cells_gpu;
    int n_rxn = n_rxn_total_cells/num_cells_gpu;
    int cell=index/n_rxn;

    int *int_data = (int *) &(((int *) int_pointer)[index%n_rxn]);
    double *float_data = (double *) &(((double *) double_pointer)[index%n_rxn]);

    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[n_rxn]);

    double *deriv_data = &( deriv_data2[deriv_length_cell*cell]);
    //state= state+(deriv_length_cell*cell); //TODO: Fix this different size with deriv maybe?
    state= state+(state_size_cell*cell);

    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        //rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(
        //        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_HL_PHASE_TRANSFER :
        //rxn_gpu_HL_phase_transfer_calc_deriv_contrib(
        //        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        //rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
        //        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
        break;

        //case RXN_ARRHENIUS :
          //rxn_gpu_arrhenius_calc_deriv_contrib(
          //        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);

//        rxn_gpu_tmp_arrhenius(
//                model_data, deriv_data, rxn_data, float_data, time_step, n_rxn);

        //break;
    }
  }
  __syncthreads();

  if (threadIdx.x < deriv_length)
  deriv[index] = deriv_data2[index];

  //if (threadIdx.x < deriv_length) {
  //  for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x)//IT WORK for large rxn input data
  //    atomicAdd(&(deriv[i_spec]), deriv_data[i_spec]);
  //}
}

void rxn_calc_deriv_gpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step) {

  // Get a pointer to the derivative data
  int num_cells = model_data->num_cells;
  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn_total_cells = rxn_data[0]*num_cells;
  double *state = model_data->state;


  /*if(countergpu==29) {
    for (int i_rxn=1; i_rxn<n_rxn; i_rxn++) {
      if(rxn_type==RXN_AQUEOUS_EQUILIBRIUM){
        printf(" type1: %d  ", int_data);
      }
    }
  }
  countergpu++;*/

  //memcpy(state_cpu, state, state_size);
  //HANDLE_ERROR(cudaMemcpy(state_gpu, state, state_size, cudaMemcpyHostToDevice));//Slower, use for large values

  //mdgpu = model_data; //Faster, use for few values
  state_gpu= state; //Faster, use for few values

  //Test to solve a bug with operations
  //rxn_gpu_tmp_arrhenius << < (n_rxn + MAX_N_GPU_THREAD - 1) / MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >> >
   //(
   //mdgpu, derivgpu_data,int_pointer_gpu, double_pointer_gpu, time_step, n_rxn

   //     mdgpu, state, derivgpu_data, time_step, NV_LENGTH_S(deriv),
   //     n_rxn, int_pointer_gpu, double_pointer_gpu, int_max_size, double_max_size
    //);

  solveRxnBlock << < (n_rxn_total_cells + MAX_N_GPU_THREAD - 1) / MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >> >
    (mdgpu, state_gpu, derivgpu_data, time_step, NV_LENGTH_S(deriv), model_data->n_state_var,
    n_rxn_total_cells, num_cells, int_pointer_gpu, double_pointer_gpu, int_max_size, double_max_size);

  cudaDeviceSynchronize();//retrieve errors (But don't retrieve anything for me)
  HANDLE_ERROR2();//0.5secs

  //HANDLE_ERROR(cudaMemcpy(deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));//0.5secs

  HANDLE_ERROR(cudaMemcpy(deriv_cpu, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost));//0.29secs
  memcpy(deriv_data, deriv_cpu, deriv_size); //This is so fast(0.01secs)
}

void free_gpu_cu() {

  //HANDLE_ERROR( cudaFreeHost( derivgpu_data ) );
  //HANDLE_ERROR( cudaFreeHost( mdgpu ) );

  HANDLE_ERROR(cudaFree( int_pointer_gpu ));
  HANDLE_ERROR(cudaFree( double_pointer_gpu ));
  HANDLE_ERROR(cudaFree(derivgpu_data));
  HANDLE_ERROR(cudaFree(dev_start_rxn_param));


  free(start_rxn_param);

  //free(int_pointer_gpu); //DO not, segmentation fault
  //free(double_pointer_gpu);//DO not, segmentation fault
  //free(deriv_cpu); DO not, segmentation fault
  //cudaFree( derivgpu_data ); //TODO: ESTO PETA
  //cudaFree( mdgpu );
  //HANDLE_ERROR( cudaFreeHost( deriv ) );

}

/*
void solveRxncpu(ModelDatagpu *model_data, double *deriv_data,
                   double time_step, int *int_data, double *float_data, int deriv_length, int n_rxn)
{

  int rxn_type = int_data[0];
  int *rxn_data = (int *) &(int_data[1]);

  switch (rxn_type) {
    case RXN_AQUEOUS_EQUILIBRIUM :
      rxn_cpu_aqueous_equilibrium_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_ARRHENIUS :
      rxn_cpu_arrhenius_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_CMAQ_H2O2 :
      rxn_cpu_CMAQ_H2O2_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_CMAQ_OH_HNO3 :
      rxn_cpu_CMAQ_OH_HNO3_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_CONDENSED_PHASE_ARRHENIUS :
      rxn_cpu_condensed_phase_arrhenius_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_EMISSION :
      rxn_cpu_emission_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_FIRST_ORDER_LOSS :
      rxn_cpu_first_order_loss_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_HL_PHASE_TRANSFER :
//rxn_cpu_HL_phase_transfer_calc_deriv_contrib(
//        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_PDFITE_ACTIVITY :
      rxn_cpu_PDFiTE_activity_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_PHOTOLYSIS :
      rxn_cpu_photolysis_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_SIMPOL_PHASE_TRANSFER :
//rxn_cpu_SIMPOL_phase_transfer_calc_deriv_contrib(
//        model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_TROE :
      rxn_cpu_troe_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_WET_DEPOSITION :
      rxn_cpu_wet_deposition_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;
    case RXN_ZSR_AEROSOL_WATER :
      rxn_cpu_ZSR_aerosol_water_calc_deriv_contrib(
              model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);
      break;

  }

}

void rxn_calc_deriv_cpu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step)
{

  realtype *deriv_data = N_VGetArrayPointer(deriv);
  int *rxn_data = (int *) (model_data->rxn_data);
  int n_rxn = rxn_data[0];

  //Case i_rxn=0
  solveRxncpu(model_data, deriv_data, time_step, int_pointer, double_pointer, NV_LENGTH_S(deriv), n_rxn);

  for (int i_rxn = 1; i_rxn < n_rxn; i_rxn++) {

    int *int_data = (int *) &(((int *) int_pointer)[int_sizes[i_rxn - 1]]);
    double *float_data = (double *) &(((double *) double_pointer)[double_sizes[i_rxn - 1]]);

    solveRxncpu(model_data, deriv_data, time_step, int_data, float_data, NV_LENGTH_S(deriv), n_rxn);
  }
}
*/


}
