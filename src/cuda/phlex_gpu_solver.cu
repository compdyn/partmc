//Test
extern "C"{
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

/*#define NUM_REACT_ int_data[0]
#define NUM_PROD_ int_data[1]
#define A_ float_data[0]
#define B_ float_data[1]
#define C_ float_data[2]
#define D_ float_data[3]
#define E_ float_data[4]
#define CONV_ float_data[5]
#define RATE_CONSTANT_ float_data[6]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 7
#define REACT_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x]-1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+(NUM_REACT_+2)*(NUM_REACT_+NUM_PROD_))
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+NUM_PROD_)
*/
//#define

const int N = 16; 
const int blocksize = 16;

//TODO: check if a double is necessary, or can we use a float, since cvode
// check a tolerance error to converge and maybe this error only reach float
// (and also the test comparison between ebi and pmc uses a tolerance, maybe float)

SolverDatagpu *sdgpu;
ModelDatagpu *mdgpu;
ModelDatagpu *mdcpu;
double * derivgpu_data; //Realtype is a double(8 bytes)
size_t deriv_size; //unsigned int
unsigned int countergpu=0;
unsigned int c1=0;
unsigned int c2=0;
unsigned int c3=0;
unsigned int c4=0;
unsigned int c5=0;
unsigned int c6=0;
unsigned int c7=0;
unsigned int c8=0;
unsigned int c9=0;
unsigned int c10=0;
unsigned int c11=0;
unsigned int c12=0;
unsigned int c13=0;
unsigned int c14=0;
size_t rxn_data_size;
size_t state_size;
void * rxn_data_gpu;

bool solver_set_gpu_sizes = 1;
int * short_pointer;
double * double_pointer;
int * short_pointer_gpu;
double * double_pointer_gpu;
unsigned int * start_rxn_param;
unsigned int * dev_start_rxn_param;
unsigned int * int_sizes;
unsigned int * int_sizes_gpu;
unsigned int * double_sizes;
unsigned int * double_sizes_gpu;


__global__ void hello(char *a, int *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
            file, line );
    exit( EXIT_FAILURE );
  }
}

void solver_new_gpu_cu(SolverDatagpu *sd, int n_dep_var,
     int n_state_var, int *var_type, int n_rxn,
     int n_rxn_int_param, int n_rxn_float_param, int n_aero_phase,
     int n_aero_phase_int_param, int n_aero_phase_float_param,
     int n_aero_rep, int n_aero_rep_int_param, int n_aero_rep_float_param,
     int n_sub_model, int n_sub_model_int_param, int n_sub_model_float_param){

  ModelDatagpu * md= &sd->model_data;
  rxn_data_size=(n_rxn_int_param + 1 + n_rxn) * sizeof(int)
                + n_rxn_float_param * sizeof(double);
  state_size=n_state_var * sizeof(double);

  // Allocate the ModelDeviceData array this works, let only this cudahostalloc of md/mdgpu activated
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &md,
                               sizeof(ModelDatagpu) ,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu ),
          (void*) md,
          0));

  //This aproax wont work for data missalignment, need cudamalloc and cudamemcpy
  /*HANDLE_ERROR( cudaHostAlloc( ( void** ) &md->rxn_data,
                               rxn_data_size,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu->rxn_data ),
          (void*) md->rxn_data,
          0));*/

  HANDLE_ERROR( cudaHostAlloc( ( void** ) &md->state,
                               state_size,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu->state ),
          (void*) md->state,
          0));

  //Pinned memory works without soa in rxn memcopying
  //HANDLE_ERROR(cudaMallocHost((void**)&(mdcpu), sizeof(ModelDatagpu)));
  //HANDLE_ERROR(cudaMallocHost((void**)&(mdcpu->rxn_data), rxn_data_size));
  //HANDLE_ERROR(cudaMallocHost((void**)&(mdgpu->state), state_size));
  //HANDLE_ERROR(cudaMallocHost((void**)&(mdcpu->rxn_data), rxn_data_size));
  //HANDLE_ERROR(cudaMallocHost((void**)&(mdcpu->state), state_size));
  cudaMalloc((void**)&(rxn_data_gpu), rxn_data_size);
  //cudaMalloc((void**)&mdgpu, sizeof(*mdgpu));

  //cudaMalloc((void**)&mdgpu, sizeof(*mdgpu));
  //cudaMalloc((void**)&(mdgpu->state), state_size);
  //cudaMalloc((void**)&(mdgpu->rxn_data), rxn_data_size);

  deriv_size = n_dep_var*sizeof(sd->y);//This size should be ok (same of chem_mod)
  cudaMalloc((void**) &derivgpu_data, deriv_size);

  if(deriv_size>MAX_SHARED_MEMORY_BLOCK_DOUBLE)
#ifndef FAILURE_DETAIL
  printf("\nMore solver variables(deriv): %d than maximum shared memory: %d",
         deriv_size,MAX_SHARED_MEMORY_BLOCK_DOUBLE);
#endif

}

void solver_update_gpu(ModelDatagpu *md){

  //This is for state array more or less

  //mdgpu->state = md->state;
  //mdgpu->rxn_data = md->rxn_data;
  //mdgpu->env = md->env;

  //memcpy( mdgpu, md, sizeof(ModelDatagpu) );
  //memcpy( mdcpu->rxn_data, md->rxn_data, rxn_data_size );
  //memcpy( mdgpu->state, md->state, state_size );

  //cudaMemcpy(mdgpu, md, sizeof(ModelDatagpu), cudaMemcpyHostToDevice );
  //cudaMemcpy( rxn_data_gpu, md->rxn_data, rxn_data_size, cudaMemcpyHostToDevice );
  //cudaMemcpy( mdgpu->state, md->state, state_size, cudaMemcpyHostToDevice );

}

void solver_set_data_gpu(ModelDatagpu *model_data) {

  //Get rxn sizes
  if(solver_set_gpu_sizes){

    // Get the number of reactions
    int *rxn_data = (int*) (model_data->rxn_data);
    int n_rxn = *(rxn_data++);
    void * rxn_param = (void*)rxn_data;
    int * float_data = (int*)rxn_data;
    unsigned int int_size = 0;
    unsigned int int_total_size = 0;
    unsigned int double_size=0;
    unsigned int double_total_size = 0;
    size_t start_size = n_rxn * sizeof(unsigned int);

    //Create started indexes of arrays
    //TTODO: Check cudahostalloc instead of alloc and memcpy
    start_rxn_param=(unsigned int*)malloc(start_size);
    int_sizes=(unsigned int*)malloc(start_size);
    double_sizes=(unsigned int*)malloc(start_size);
    cudaMalloc((void**) &dev_start_rxn_param, start_size);
    cudaMalloc((void**) &int_sizes_gpu, start_size);
    cudaMalloc((void**) &double_sizes_gpu, start_size);

    for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

      //Reaction distances between pointers rows
      start_rxn_param[i_rxn]=(unsigned int) ((int*) rxn_data - (int*) rxn_param);

      int *rxn_start = rxn_data;

      // Get the reaction type
      int rxn_type = *(rxn_data++);

      //TTODO: Join functions to retrieve float and int sizes
      // Call the appropriate function
      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          float_data = (int*) rxn_gpu_aqueous_equilibrium_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_aqueous_equilibrium_skip(
                  (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          float_data = (int*) rxn_gpu_arrhenius_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_arrhenius_skip(
                  (void*) rxn_data);
          break;
        case RXN_CMAQ_H2O2 :
          float_data = (int*) rxn_gpu_CMAQ_H2O2_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_CMAQ_H2O2_skip(
                  (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          float_data = (int*) rxn_gpu_CMAQ_OH_HNO3_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_CMAQ_OH_HNO3_skip(
                  (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          float_data = (int*) rxn_gpu_condensed_phase_arrhenius_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_condensed_phase_arrhenius_skip(
                  (void*) rxn_data);
          break;
        case RXN_EMISSION :
          float_data = (int*) rxn_gpu_emission_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_emission_skip(
                  (void*) rxn_data);
          break;
        case RXN_FIRST_ORDER_LOSS :
          float_data = (int*) rxn_gpu_first_order_loss_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_first_order_loss_skip(
                  (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          float_data = (int*) rxn_gpu_HL_phase_transfer_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_HL_phase_transfer_skip(
                  (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          float_data = (int*) rxn_gpu_PDFiTE_activity_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_PDFiTE_activity_skip(
                  (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          float_data = (int*) rxn_gpu_photolysis_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_photolysis_skip(
                  (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          float_data = (int*) rxn_gpu_SIMPOL_phase_transfer_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_SIMPOL_phase_transfer_skip(
                  (void*) rxn_data);
          break;
        case RXN_TROE :
          float_data = (int*) rxn_gpu_troe_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_troe_skip(
                  (void*) rxn_data);
          break;
        case RXN_WET_DEPOSITION :
          float_data = (int*) rxn_gpu_wet_deposition_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_wet_deposition_skip(
                  (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          float_data = (int*) rxn_gpu_ZSR_aerosol_water_int_size(
                  (void*) rxn_data);
          rxn_data = (int*) rxn_gpu_ZSR_aerosol_water_skip(
                  (void*) rxn_data);
          break;
      }

      int_size = (unsigned int)((int*) float_data - (int*) rxn_start );
      int_total_size += int_size;
      int_sizes[i_rxn]=int_total_size;

      double_size=(unsigned int)((double*) rxn_data - (double*) float_data);
      double_total_size += double_size;
      double_sizes[i_rxn]=double_total_size;

    }

    //Allocate int and double rxn data separately
    short_pointer = (int*) malloc(int_total_size * sizeof(int));
    double_pointer = (double*) malloc(double_total_size * sizeof(double));
    cudaMalloc((void**) &short_pointer_gpu, int_total_size * sizeof(int));
    cudaMalloc((void**) &double_pointer_gpu, double_total_size * sizeof(double));

    for (int j = 0; j < int_sizes[0]; j++)
      short_pointer[j] = ((int*) rxn_param)[j];

    for (int j = 0; j < double_sizes[0]; j++){
      double *float_data = (double*) &(((int*)rxn_param)[
              int_sizes[0]]);
      double_pointer[j] = float_data[j];
    }

    for (int i_rxn=1; i_rxn<n_rxn; i_rxn++) {
      //Get rxn data
      //TODO: use short_pointer in reactions

      int_size= int_sizes[i_rxn]-int_sizes[i_rxn-1];
      double_size= double_sizes[i_rxn]-double_sizes[i_rxn-1];

      for (int j = 0; j < int_size; j++)
        short_pointer[int_sizes[i_rxn-1] + j] = ((int *) rxn_param)[start_rxn_param[i_rxn]+j];

      for (int j = 0; j < double_size; j++){
        double *float_data = (double*) &(((int*)rxn_param)[start_rxn_param[i_rxn]+int_size]);
        double_pointer[double_sizes[i_rxn-1] + j] = float_data[j];
      }
    }

    HANDLE_ERROR(cudaMemcpy( short_pointer_gpu, short_pointer, int_total_size * sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy( double_pointer_gpu, double_pointer, double_total_size * sizeof(double), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy( int_sizes_gpu, int_sizes, n_rxn*sizeof(unsigned int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy( double_sizes_gpu, double_sizes, n_rxn*sizeof(unsigned int), cudaMemcpyHostToDevice));


    //TTODO: Create global size arrays (int and float) with each size for reaction instead the defines
    //TODO: Secure 100% size not changes on every iteration on monarch, it shouldn't since in new solver we allocate rxn one time
    solver_set_gpu_sizes=0;
  }

}

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *deriv,//double *deriv
        double time_step,
        int deriv_length, void * rxn_data_gpu2, int n_rxn,
        int * short_pointer, double * double_pointer,
        unsigned int * int_sizes, unsigned int * double_sizes) //Interface CPU/GPU
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int * int_data;
  double *float_data;

  __shared__ double deriv_data[ MAX_SHARED_MEMORY_BLOCK_DOUBLE ];

  //if(index==2)
  //for (int i=0; i<deriv_length; i++)
    //deriv_data[i]=model_data->state[0];

  if(index<deriv_length)
    deriv_data[index]=deriv[index];

  __syncthreads();

  if(index==0){
    int_data  = (int*) short_pointer;
    float_data = (double*) double_pointer;

  }else {
    if (index < n_rxn) {
      int_data = (int *) &(((int *) short_pointer)[int_sizes[index - 1]]);
      float_data = (double *) &(((double *) double_pointer)[double_sizes[index - 1]]);
    }
  }

  if(index < n_rxn){

    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[1]);

    //int * rxn_data2 = (int*) &(((int*) rxn_data_gpu2)[dev_start_rxn_param[ 2 ]] );

    //void * rxn_data3 = (void*)&(((int*) &(((int*) rxn_data_gpu2)[dev_start_rxn_param[ index ]] ))[1]);
    //void * rxn_data3 = (void*)&(((int*) &(((int*) rxn_data_gpu2)[dev_start_rxn_param[ index ]] ))[0]);

    //int * rxn_data2 = ((int*)dev_rxn_param)+dev_start_rxn_param[index];
    //int * rxn_data= (rxn_data2++);
    //int rxn_type = *((int*) &(((int*) dev_rxn_param)[dev_start_rxn_param[ index ]] ));

    //int * rxn_data = (int*)&(rxn_data2[1]);
    //int rxn_type = *(((int*)dev_rxn_param)+dev_start_rxn_param[index]);

    //double *float_data =  (double*) &(((int*)rxn_data)[
    //      2+(((int*)rxn_data)[0]+2)*(((int*)rxn_data)[0]+((int*)rxn_data)[1]) ]);

    //atomicAdd(&(deriv_data[0]),float_data[6]);
    //}

    //Notice till this, without executing calc_deriv switch, gpu version takes 2 secs more
    // Call the appropriate function
    switch (rxn_type) {
      /*case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_gpu_HL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, double_pointer_gpu, time_step);
        break;*/
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(
              model_data, deriv_data, (void*) rxn_data, float_data, time_step);


        //double *float_data =  (double*) &(rxn_data[
        //        2+(rxn_data[0]+2)*(rxn_data[0]+rxn_data[1]) ]);//fine

        //if (rate)
        //if(index<deriv_length)
        //for ( int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x )
        //atomicAdd(&(deriv[index]), float_data[6]);
        break;
    }
  }
  __syncthreads();
  for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x)// add also the
    atomicAdd(&( deriv[ i_spec ] ), deriv_data[ i_spec ] );
}

void rxn_calc_deriv_gpu_cu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step){

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);
  void * rxn_param = (void*)rxn_data;

  //TODO: Only 1 time call reaction_int_size,get int sizes, assign values to rxn_datagpu with only doubles (no array of ints and doubles)
  //y refactor de int *int_data = (int*) rxn_data; a cast de double para k cuadre en gpu

  /*if(countergpu==29) {
    for (int i_rxn=1; i_rxn<n_rxn; i_rxn++) {

      int * int_data = (int*) &(((int*) rxn_param)[start_rxn_param[ i_rxn ]] );
      int rxn_type = int_data[0];
      int * rxn_data = (int*)&(int_data[1]);
     // int int_data = *((int *) &(((int *) short_pointer)[int_sizes[i_rxn - 1]]));

      if(rxn_type==RXN_ARRHENIUS){
        double *float_data =  (double*) &(((int*)rxn_data)[
                2+(((int*)rxn_data)[0]+2)*(((int*)rxn_data)[0]+((int*)rxn_data)[1]) ]);
        double * float_data2 = (double *) &(((double *) double_pointer)[double_sizes[i_rxn-1]]);
        printf(" type1: %f  ", float_data[6]);
        printf(" type2: %f\n", float_data2[6]);
      }
    }
  }*/
  countergpu++;

//TODO: Init deriv to 0 inside gpu instead of memcpy
  //cudaMemcpy( dev_start_rxn_param, start_rxn_param, start_size, cudaMemcpyHostToDevice );

  //HANDLE_ERROR(cudaMemcpy( mdgpu->rxn_data, model_data->rxn_data, rxn_data_size, cudaMemcpyHostToDevice ));
  //HANDLE_ERROR(cudaMemcpy( rxn_data_gpu2, (void*)rxn_data2, rxn_data_size- sizeof(int), cudaMemcpyHostToDevice ));

  //void * rxn_data_gpu;
  //cudaMalloc((void**)&(rxn_data_gpu), rxn_data_size);
  //HANDLE_ERROR(cudaMemcpy( rxn_data_gpu, (void*)rxn_param, rxn_data_size, cudaMemcpyHostToDevice ));//fix this size

  //cudaMemcpy( mdgpu->rxn_data, model_data->rxn_data, rxn_data_size, cudaMemcpyHostToDevice );
  //HANDLE_ERROR(cudaMemcpy( mdgpu->state, mdcpu->state, state_size, cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy( &mdgpu, &model_data, rxn_data_size, cudaMemcpyHostToDevice));
  //cudaMemcpy( &(mdgpu->rxn_data), &(model_data->rxn_data), sizeof(model_data->rxn_data), cudaMemcpyHostToDevice );
  //cudaMemcpy( &(mdgpu->state), &(model_data->state), sizeof(model_data->state), cudaMemcpyHostToDevice );

  mdgpu=model_data;
  //mdgpu->state = model_data->state;
  //mdgpu->rxn_data = model_data->rxn_data;

  //int *rxn_data2 = (int*) (mdcpu->rxn_data);
  //int n_rxn2 = *(rxn_data2++);
  //void *dev_rxn_param= (void*)rxn_data2;

  cudaMemcpy( derivgpu_data, deriv_data, deriv_size, cudaMemcpyHostToDevice );
  solveRxnBlock<<< (n_rxn+MAX_N_GPU_THREAD-1)/MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >>>
      ( mdgpu, derivgpu_data, time_step, NV_LENGTH_S(deriv), //deriv_size/sizeof(double)//rxn_data_gpu
      (void*)rxn_data_gpu, n_rxn, short_pointer_gpu, double_pointer_gpu,
      int_sizes_gpu, double_sizes_gpu);
  //cudaDeviceSynchronize();//retrieve errors (But don't retrieve anything for me)
  cudaMemcpy( deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost );

  //solveRxnBlockCpu( mdgpu, deriv_data, time_step, start_rxn_param,
  //            NV_LENGTH_S(deriv), (void*)rxn_data2, n_rxn, short_pointer, double_pointer,
  //            int_sizes, double_sizes); //deriv_size/sizeof(double)

//    for (int i=0; i<NV_LENGTH_S(deriv); i++) {
     
   //NV_DATA_S(deriv)[i]=0.01;
   // deriv_data[i]=0.01;
	//printf(" deriv gpu: % -le\n", NV_DATA_S(derivgpu)[i]);
  //  }

  //free(start_rxn_param);
  //cudaFree( dev_start_rxn_param );
}

void free_gpu_cu() {

  //HANDLE_ERROR( cudaFreeHost( derivgpu_data ) );
  //HANDLE_ERROR( cudaFreeHost( mdgpu ) );

  //cudaFree( derivgpu_data ); //TODO: ESTO PETA 
  //cudaFree( mdgpu );
  //HANDLE_ERROR( cudaFreeHost( deriv ) );

}

void printfCUDA(int aggg)
{
	char a[N] = "Hello \0\0\0\0\0\0";
	int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	char *ad;
	int *bd;
	const int csize = N*sizeof(char);
	const int isize = N*sizeof(int);

	printf("%s", a);

	cudaMalloc( (void**)&ad, csize );
	cudaMalloc( (void**)&bd, isize );
	cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice );
	cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice );

	dim3 dimBlock( blocksize, 1 );
	dim3 dimGrid( 1, 1 );
	hello<<<dimGrid, dimBlock>>>(ad, bd);
	cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost );
	cudaFree( ad );
	cudaFree( bd );

	printf("%s nano\n", a);
}
} //extern C
//extern "C"
//void printfCUDA()
//{printf("hola\n");}