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

//#define

const int N = 16; 
const int blocksize = 16;

//TODO: check if a double is necessary, or can we use a float, since cvode
// check a tolerance error to converge and maybe this error only reach float
// (and also the test comparison between ebi and pmc uses a tolerance, maybe float)

//TTODO: Global variables? ew
SolverDatagpu *sdgpu;
ModelDatagpu *mdgpu;
double * derivgpu_data; //Realtype is a double(8 bytes)
unsigned short deriv_size;
unsigned short countergpu=0;
unsigned short c1=0;
unsigned short c2=0;
unsigned short c3=0;
unsigned short c4=0;
unsigned short c5=0;
unsigned short c6=0;
unsigned short c7=0;
unsigned short c8=0;
unsigned short c9=0;
unsigned short c10=0;
unsigned short c11=0;
unsigned short c12=0;
unsigned short c13=0;
unsigned short c14=0;
unsigned int rxn_data_size;
unsigned int state_size;

__global__ void hello(char *a, int *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}

//extern "C"

//Planning: Add pointer to gpu data allocated in sd pointer
//or call a new gpu_new function on fortran, which call to a .c function and this function to cu (I prefer this in fact instead of be annoying with sd pointers, more modular)

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
/*
  sdgpu = (SolverDatagpu*)malloc(sizeof(SolverDatagpu));
    //SolverData *sd = (SolverData*) malloc(sizeof(SolverData));
  if (sdgpu==NULL) {
    printf("\n\nERROR allocating space for SolverDatagpu\n\n");
    exit(1);
  }

  sdgpu = (SolverDatagpu*)sd;
  sdgpu->curr_J_guess = true;
  printf(" soy yo el test: %i\n", sdgpu->curr_J_guess);
*/
/*
  // Allocate the ModelDeviceData array
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &sd,
                               sizeof(SolverDatagpu) ,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer(
          (void**) &( sdgpu ),
          (void*) sd,
          0));
*/

  ModelDatagpu * md= &sd->model_data;
  rxn_data_size=(n_rxn_int_param + 1 + n_rxn) * sizeof(int)
                + n_rxn_float_param * sizeof(double);
  state_size=n_state_var * sizeof(double);

  // Allocate the ModelDeviceData array
  /*HANDLE_ERROR( cudaHostAlloc( ( void** ) &md,
                               sizeof(ModelDatagpu) ,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu ),
          (void*) md,
          0));

  // Allocate the ModelDeviceData array
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &md->rxn_data,
                               rxn_data_size,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu->rxn_data ),
          (void*) md->rxn_data,
          0));

  HANDLE_ERROR( cudaHostAlloc( ( void** ) &md->state,
                               state_size,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu->state ),
          (void*) md->state,
          0));*/

  //Pinned memory
  HANDLE_ERROR(cudaMallocHost((void**)&mdgpu, sizeof(ModelDatagpu)));
  //HANDLE_ERROR(cudaMallocHost((void**)&md->rxn_data, rxn_data_size));
  HANDLE_ERROR(cudaMallocHost((void**)&mdgpu->state, state_size));

  //mdgpu->state=md->state;
  //cudaMalloc((void**)&mdgpu, sizeof(*mdgpu));
  //cudaMalloc((void**)&(mdgpu->state), state_size);

  deriv_size = n_dep_var*sizeof(sd->y);//This size should be ok (same of chem_mod)
  cudaMalloc((void**) &derivgpu_data, deriv_size);

  //printf(" soy yo el test: %d\n", deriv_size);

  if(deriv_size>MAX_SHARED_MEMORY_BLOCK_DOUBLE)
//#ifndef FAILURE_DETAIL
  printf("\nMore solver variables(deriv): %d than maximum shared memory: %d",
         deriv_size,MAX_SHARED_MEMORY_BLOCK_DOUBLE);
//#endif
}

void solver_update_gpu(ModelDatagpu *md){

  //mdgpu->state = md->state;
  //mdgpu->env = md->env;

  //Pinned
  //memcpy( mdgpu->rxn_data, md->rxn_data, rxn_data_size );
  memcpy( mdgpu->state, md->state, state_size );

  //cudaMemcpy(mdgpu, md, sizeof(ModelDatagpu), cudaMemcpyHostToDevice );
  //cudaMemcpy( mdgpu->rxn_data, md->rxn_data, rxn_data_size, cudaMemcpyHostToDevice );
  //cudaMemcpy( mdgpu->state, md->state, state_size, cudaMemcpyHostToDevice );

  //cudaMemcpy( mdgpu->rxn_data, md->rxn_data, rxn_data_size, cudaMemcpyHostToDevice );
  //cudaMemcpy( &(mdgpu->state), &md->state, sizeof(mdgpu->state), cudaMemcpyHostToDevice );

}

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *deriv,//double *deriv
        double time_step,unsigned int * dev_start_rxn_param,
        int deriv_length) //Interface CPU/GPU
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  // Get the number of reactions
  //int *dev_rxn_param = (int*) (model_data->rxn_data);
  //int n_rxn = *(dev_rxn_param++);

  __shared__ double deriv_data[ MAX_SHARED_MEMORY_BLOCK_DOUBLE ];

  //for ( int i_spec = threadIdx.x; i_spec < MAX_SHARED_MEMORY_BLOCK_DOUBLE; i_spec ++ )
    //deriv_data[ i_spec ] = 0.01;
  //TTODO: Check if works fusioning this two fors

  if(threadIdx.x==1)
  for (int i=0; i<deriv_length; i++)
    //deriv_data[i]=0.000000001;
    deriv_data[i]=model_data->state[0];

  __syncthreads();
   
  for ( int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x )
    atomicAdd( &( deriv_data[ i_spec ] ), deriv[ i_spec ] );
 
  //if(index<NV_LENGTH_S(deriv)){
    //NV_DATA_S(deriv)[threadIdx.x]=0.01;
    //printf(" deriv cpu: %le \n", NV_DATA_S(deriv)[i]);
  //}

  __syncthreads();

  // Loop through the reactions advancing the rxn_data pointer each time
  /*if(index<n_rxn) { //Maybe not necessary

    // Get the reaction type
    //int rxn_type = *(rxn_data);
    //int rxn_type = (void*) &(((int*) dev_rxn_param)[dev_start_rxn_param[ index ]] );

    //This don't need to work since the data type is different
    //int * rxn_data = (int*) &(((int*) dev_rxn_param)[dev_start_rxn_param[ index ]] );
    int * rxn_data = ((int*)dev_rxn_param)+dev_start_rxn_param[index];

    //int * rxn_data = (int*) dev_rxn_param[dev_start_rxn_param[ index ]];
    int rxn_type = *(rxn_data++);

    //Notice till this, without executing calc_deriv switch, gpu version takes 2 secs more
    // Call the appropriate function
    switch (rxn_type) {*/
      /*case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_gpu_aqueous_equilibrium_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_gpu_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_gpu_CMAQ_H2O2_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_gpu_CMAQ_OH_HNO3_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_gpu_condensed_phase_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_EMISSION :
        rxn_gpu_emission_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_gpu_first_order_loss_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_gpu_HL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_gpu_PDFiTE_activity_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_gpu_photolysis_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_gpu_SIMPOL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_gpu_troe_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_gpu_wet_deposition_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_gpu_ZSR_aerosol_water_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;

	/*case RXN_ARRHENIUS :
        //rxn_gpu_arrhenius_calc_deriv_contrib(
          //      model_data, deriv_data, (void*) rxn_data, time_step);
        break;

     }
  }*/

  __syncthreads();

  //if(threadIdx.x<deriv_length)
  for (int i_spec = threadIdx.x; i_spec < deriv_length; i_spec += blockDim.x)// add also the
    atomicAdd(&( deriv[ i_spec ] ), deriv_data[ i_spec ] );

}



void rxn_calc_deriv_gpu_cu(ModelDatagpu *model_data, N_Vector deriv, realtype time_step){

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  //Create started indexes of arrays
  size_t size = n_rxn * sizeof(unsigned int);
  unsigned int *start_rxn_param;
  unsigned int *dev_start_rxn_param;
  HANDLE_ERROR( cudaHostAlloc( (void**) &start_rxn_param, size,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_start_rxn_param,
                                          start_rxn_param, 0 ) );

  void * rxn_param = rxn_data;
  //TTODO: Put this loop at new_phlex_solver since data matrix don't differ in iterations
  // Loop through the reactions advancing the rxn_data pointer each time
 for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    //Reaction distances between pointers rows
    start_rxn_param[i_rxn]=(unsigned int) ( (int*) rxn_data - (int*) rxn_param );
    // Get the reaction type
    int rxn_type = *(rxn_data++);


    /*if(countergpu==29){
    //printf(" type: %d\n", start_rxn_param[i_rxn]);//makes sense
      int * rxn_data2 = ((int*)rxn_param)+start_rxn_param[i_rxn];
      //int * rxn_data = (int*) dev_rxn_param[dev_start_rxn_param[ index ]];
      int rxn_type2 = *(rxn_data2++);
      printf(" type1: %d   ", rxn_data2);//it works
      printf(" type2: %d\n", rxn_data);
    }*/

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_gpu_aqueous_equilibrium_skip(
                (void*) rxn_data);c1++;
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_arrhenius_skip(
                (void*) rxn_data);c2++;
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_gpu_CMAQ_H2O2_skip(
                (void*) rxn_data);c3++;
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_gpu_CMAQ_OH_HNO3_skip(
                (void*) rxn_data);c4++;
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_condensed_phase_arrhenius_skip(
                (void*) rxn_data);c5++;
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_gpu_emission_skip(
                (void*) rxn_data);c6++;
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_gpu_first_order_loss_skip(
                (void*) rxn_data);c7++;
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_gpu_HL_phase_transfer_skip(
                (void*) rxn_data);c8++;
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_gpu_PDFiTE_activity_skip(
                (void*) rxn_data);c9++;
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_gpu_photolysis_skip(
                (void*) rxn_data);c10++;
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_gpu_SIMPOL_phase_transfer_skip(
                (void*) rxn_data);c11++;
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_gpu_troe_skip(
                (void*) rxn_data);c12++;
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_gpu_wet_deposition_skip(
                (void*) rxn_data);c13++;
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_gpu_ZSR_aerosol_water_skip(
                (void*) rxn_data);c14++;
        break;
    }
  }
  /*if(countergpu==10000){
    printf("types:");
    printf(", %hu:",&c1);
    printf(", %hu:",&c2);
    printf(", %hu:",&c3);
    printf(", %hu:",&c4);
    printf(", %hu:",&c5);
    printf(", %hu:",&c6);
    printf(", %hu:",&c7);
    printf(", %hu:",&c8);
    printf(", %hu:",&c9);
    printf(", %hu:",&c10);
    printf(", %hu:",&c11);
    printf(", %hu:",&c12);
    printf(", %hu:",&c13);
    printf(", %hu:",&c14);
    printf("\n");
  }*/
  if(countergpu==29) {
    //double *state = model_data->state;

    //for (int i=0; i<2; i++)
    printf(" state: %f\n", model_data->state[0]);

  }
  countergpu++;

//TODO: Init deriv to 0 inside gpu instead of memcpy
//deriv_size
  cudaMemcpy( derivgpu_data, deriv_data, deriv_size, cudaMemcpyHostToDevice );

  solveRxnBlock<<< (n_rxn+MAX_N_GPU_THREAD-1)/MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >>>
      ( mdgpu, derivgpu_data, time_step, dev_start_rxn_param,
      NV_LENGTH_S(deriv)); //deriv_size/sizeof(double)

  cudaMemcpy( deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost );

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
