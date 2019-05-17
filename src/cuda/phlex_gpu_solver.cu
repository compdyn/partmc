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
double* derivgpu_data; //Realtype is a double(8 bytes)
unsigned short deriv_size;
unsigned short counter=0;

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

void solver_new_gpu_cu(SolverDatagpu *sd, int n_dep_var){
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

  // Allocate the ModelDeviceData array
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &md,
                               sizeof(ModelDatagpu) ,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( mdgpu ),
          (void*) md,
          0));
/*
  N_Vector * deriv2;
  N_Vector aux = N_VNew_Serial(n_dep_var);
  deriv2 = &aux;

  // Allocate the ModelDeviceData array
  HANDLE_ERROR( cudaHostAlloc( ( void** ) &deriv2,
                               sizeof(N_Vector) ,
                               cudaHostAllocWriteCombined |
                               cudaHostAllocMapped
  ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( //Connected gpu and cpu data in principle
          (void**) &( derivgpu ),
          (void*) deriv2,
          0));
*/
  deriv_size = n_dep_var*sizeof(sd->y);//This size should be ok (same of chem_mod)
  cudaMalloc(&derivgpu_data, deriv_size);

  //printf(" soy yo el test: %d\n", deriv_size);

  if(deriv_size>MAX_SHARED_MEMORY_BLOCK_DOUBLE)
#ifndef FAILURE_DETAIL
  printf("\nMore solver variables(deriv): %d than maximum shared memory: %d",
         deriv_size,MAX_SHARED_MEMORY_BLOCK_DOUBLE);
#endif
}

void solver_initialize_gpu_cu(void *sd) {



}

void solver_run_gpu_cu(void *sd) {
}

__global__ void solveRxnBlock(ModelDatagpu *model_data, double *deriv,
        void *dev_rxn_param, double time_step,unsigned int * dev_start_rxn_param
        ,unsigned short deriv_size_gpu) //Interface CPU/GPU
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data);

  __shared__ double deriv_data[ MAX_SHARED_MEMORY_BLOCK_DOUBLE ];

  for ( int i_spec = threadIdx.x; i_spec < deriv_size_gpu; i_spec += blockDim.x )
    deriv_data[ i_spec ] = 0.0;
  //TTODO: Check if works fusioning this two fors
  __syncthreads();

  for ( int i_spec = threadIdx.x; i_spec < deriv_size_gpu; i_spec += blockDim.x )
    atomicAdd( &( deriv_data[ i_spec ] ), deriv[ i_spec ] );

  __syncthreads();

  // Loop through the reactions advancing the rxn_data pointer each time
  if(index<n_rxn) { //Maybe not necessary

    // Get the reaction type
    //int rxn_type = *(rxn_data);
    //int rxn_type = (void*) &(((int*) dev_rxn_param)[dev_start_rxn_param[ index ]] );
    rxn_data = (int*) &(((int*) dev_rxn_param)[dev_start_rxn_param[ index ]] );
    int rxn_type = *(rxn_data++);

    //Notice till this, without executing calc_deriv switch, gpu version takes 2 secs more
    /*
    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
    }
     */
  }

  __syncthreads();

  for (int i_spec = threadIdx.x; i_spec < deriv_size_gpu; i_spec += blockDim.x)// add also the
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

  int * rxn_param = rxn_data;
  //TTODO: Put this loop at new_phlex_solver since data matrix don't differ in iterations
  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    //Reaction distances between pointers rows
    start_rxn_param[i_rxn]=(unsigned int) ( (int*) rxn_data - (int*) rxn_param );
    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_gpu_aqueous_equilibrium_skip(
                (void*) rxn_data);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_arrhenius_skip(
                (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_gpu_CMAQ_H2O2_skip(
                (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_gpu_CMAQ_OH_HNO3_skip(
                (void*) rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_condensed_phase_arrhenius_skip(
                (void*) rxn_data);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_gpu_emission_skip(
                (void*) rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_gpu_first_order_loss_skip(
                (void*) rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_gpu_HL_phase_transfer_skip(
                (void*) rxn_data);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_gpu_PDFiTE_activity_skip(
                (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_gpu_photolysis_skip(
                (void*) rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_gpu_SIMPOL_phase_transfer_skip(
                (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_gpu_troe_skip(
                (void*) rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_gpu_wet_deposition_skip(
                (void*) rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_gpu_ZSR_aerosol_water_skip(
                (void*) rxn_data);
        break;
    }
  }
  ////cudaMemcpy( dev_deriv, deriv, deriv_size, cudaMemcpyHostToDevice ); //crec que no fa falta si ta link el punter de dades

  //N_Vector aux= *derivgpu;
  //realtype *derivgpu_data = N_VGetArrayPointer(aux);


  /*if(counter==1){

    //printf(" soy yo el test: %d\n", sizeof(N_Vector));
    //deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
    //printf(" soy yo el test: %f\n", deriv_data[0]);
    //derivgpu_data=(double*)deriv_data;
    //printf(" soy yo el test: %f\n", derivgpu_data[2]);

  }
  counter++;*/
//TODO: Init deriv to 0 inside gpu instead of memcpy
  cudaMemcpy( derivgpu_data, deriv_data, deriv_size, cudaMemcpyHostToDevice );

      solveRxnBlock<<< (n_rxn+MAX_N_GPU_THREAD-1)/MAX_N_GPU_THREAD, MAX_N_GPU_THREAD >>>
      ( model_data, derivgpu_data, (void*) rxn_param, time_step, dev_start_rxn_param,
      deriv_size);

  cudaMemcpy( deriv_data, derivgpu_data, deriv_size, cudaMemcpyDeviceToHost );

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
