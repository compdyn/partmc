
/** GPU Solver test
 *
 * Simulates phlex-chem time derivative function with GPU support
 */
#include <cuda.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/// Maximum number of GPU blocks to use in solver
#define MAX_N_GPU_BLOCK 10
/// Max Number of GPU threads to use in each block
#define MAX_N_GPU_THREAD 600
/// Number of repeats to perform for each configuration for time averaging
#define N_REPEAT 1000 //In fact, its not 1000 repeat, it's 34k repeats...
/// Small number for floating point comparisons
#define SMALL_NUMBER 1.0e-5
/// Number of state variables (species) / Deriv and state array size
#define N_SPEC 100
/// Number of reactions
#define N_RXN 5000
/// Minimum number of species involved in a reaction
#define MIN_RXN_SIZE 3
/// Maximum number of species involved in a reaction
#define MAX_RXN_SIZE 30

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

/** Print a row of derivatives
 *
 * @param n_gpus Number of blocks
 * @param n_threads Number of threads
 * @param comp_time Computational time (s)
 * @param deriv Derivative array
 */
void printDeriv( int n_blocks, int n_threads, float comp_time, float * deriv )
{
  printf( "%d\t%d\t%f", n_blocks, n_threads, comp_time );
//  for( int i_spec = 0; i_spec < N_SPEC; i_spec++ )
//    printf( "\t%le", deriv[ i_spec ] );
  printf( "\n" );
}

/** Assert that two calculated derivative arrays are essentially equal
 *
 * @param deriv1 First array derivative to compare
 * @param deriv2 Second derivative array to compare
 */
void assertEqual( float * deriv1, float * deriv2 )
{
  for ( int i_spec = 0; i_spec < N_SPEC; i_spec++ ) {
    if ( deriv1[ i_spec ] == 0.0 ) {
      if ( deriv2[ i_spec ] == 0.0 ) continue;
      printf( "\n\nMismatched derivative element %d: %le %le\n\n",
          i_spec, deriv1[ i_spec ], deriv2[ i_spec ] );
      exit( 1 );
    }
    if ( fabs( ( deriv1[ i_spec ] - deriv2[ i_spec ] ) / deriv1[ i_spec ] ) <
        SMALL_NUMBER ) continue;
    printf( "\n\nMismatched derivative element %d: %le %le\n\n",
        i_spec, deriv1[ i_spec ], deriv2[ i_spec ] );
    exit( 1 );
  }
}

/** Set up a set of random reactions
 *
 * @param rxn_param Pointer to hold the reaction parameters
 * @param dev_rxn_param Pointer to hold the reaction parameters on the GPUs
 */
void newReactionSet( void ** rxn_param, void ** dev_rxn_param )
{
  // Set the size of each reaction and sum the total size for all rxns
  int rxn_size[ N_RXN ];
  int rxn_size_int = 0;
  int rxn_size_float = 0;
  for( int i_rxn = 0; i_rxn < N_RXN; i_rxn++ ) {
    rxn_size[ i_rxn ] = ( rand() % ( MAX_RXN_SIZE + 1 - MIN_RXN_SIZE ) ) +
                        MIN_RXN_SIZE;
    // sum the total size including the space required to align the float values
    rxn_size_int += rxn_size[ i_rxn ] +
	( ( rxn_size[ i_rxn ] * sizeof(int) ) % sizeof(double) ) / sizeof(int);
    rxn_size_float += rxn_size[ i_rxn ];
  }

  // Allocate space for the reaction parameters
  // Parameter block holds size (i.e. number of species involved) in each rxn,
  // the index of the beginning of the floating point values in the rxn data,
  // the index in the state array of each species in the reaction, and a
  // floating-point value to multiply the species concentration by during
  // solving.
  size_t rxn_param_size = 2 * N_RXN * sizeof(int) +
                          rxn_size_int * sizeof(int) +
                          rxn_size_float * sizeof(double);
  HANDLE_ERROR( cudaHostAlloc( (void**) rxn_param,
				rxn_param_size,
				cudaHostAllocWriteCombined |
					cudaHostAllocMapped ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( dev_rxn_param, *rxn_param, 0 ) );

  // Set the values of the reaction parameters
  int * rxn_ptr = (int*) *rxn_param;
  for ( int i_rxn = 0; i_rxn < N_RXN; i_rxn++ ) {
    int * rxn_start = rxn_ptr;
    *(rxn_ptr++) = rxn_size[ i_rxn ];
    *(rxn_ptr) = rxn_size[ i_rxn ] + 2 +
	( ( rxn_size[ i_rxn ] * sizeof(int) ) % sizeof(double) ) / sizeof(int);
    double * rxn_flt_ptr = (double*) &( rxn_start[ *rxn_ptr++ ] );
    for ( int i_spec = 0; i_spec < rxn_size[ i_rxn ]; i_spec++ )
      *(rxn_ptr++) = rand() % ( N_SPEC + 1 );
    for ( int i_spec = 0; i_spec < rxn_size[ i_rxn ]; i_spec++ )
      *(rxn_flt_ptr++) = 1.0e-3 * ( rand() % 2000000 );
    rxn_ptr = (int*) rxn_flt_ptr;
  }
}

/** Solve for a single reaction
 *
 * @param rxn_ptr Pointer to the reaction parameters. Will be advanced to next
 *                reaction data
 * @param state Pointer to the state array
 * @param deriv Pointer to the derivative array
 */
void solveRxn( void ** rxn_ptr, float * state, float * deriv)
{
  int n_rxn_spec = ( (int*) *rxn_ptr )[0];
  int dbl_index = ( (int*) *rxn_ptr )[1];
  int * spec_ids = &( ( (int*) *rxn_ptr )[2] );
  double * rxn_params = (double*) &( ( (int*) *rxn_ptr )[ dbl_index ] );

  // The actual operations performed depend on the reaction type, but only
  // involve reading from this reaction's parameters in rxn_ptr, reading the
  // value of any number of elements in the state array and adding values to
  // any number of elements in the deriv array.
  for ( int i_spec = 0; i_spec < n_rxn_spec; i_spec++ ) {
    deriv[ spec_ids[ i_spec ] ] += state[ spec_ids[ i_spec ] ] *
                                   rxn_params[ i_spec ];
  }

  *rxn_ptr = (void*) &( rxn_params[ n_rxn_spec ] );
}

/** Solve the mechanism without GPUs
 *
 * @param state Pointer to the state array
 * @param rxn_param Pointer to the reaction parameter block
 * @param deriv Pointer to the derivative array to calculate
 * @return Computational time (s)
 */
float solveNoGpus( float * state, void * rxn_param, float * deriv )
{
  clock_t comp_begin = clock();

  for (int i_repeat = 0; i_repeat < N_REPEAT; i_repeat++ ) {
    for (int i_spec = 0; i_spec < N_SPEC; i_spec++ ) deriv[ i_spec ] = 0.0;
    void * rxn_ptr = rxn_param;
    for (int i_rxn = 0; i_rxn < N_RXN; i_rxn++ )
      solveRxn( &rxn_ptr, state, deriv );
  }

  clock_t comp_end = clock();
  return ( (float) ( comp_end - comp_begin ) ) / CLOCKS_PER_SEC;
}

/** Skip over a reaction advancing the reaction data pointer
 *
 * @param rxn_ptr Pointer to the reactions parameters
 * @return rxn_ptr advanced by the size of the reaction data
 */
void * skipRxn( void * rxn_ptr )
{
  int n_rxn_spec = ( (int*) rxn_ptr )[0];
  int dbl_index = ( (int*) rxn_ptr )[1];
  double * rxn_params = (double*) &( ( (int*) rxn_ptr )[ dbl_index ] );
  return (void*) &( rxn_params[ n_rxn_spec ] );
}

/** Solve for a single reaction on the device
 *
 * @param rxn_ptr Pointer to the reaction parameters. Will be advanced to next
 *                reaction data
 * @param state Pointer to the state array
 * @param deriv Pointer to the derivative array
 */
__device__ void devSolveRxn( void ** rxn_ptr, float * state, float * deriv)
{
  int n_rxn_spec = ( (int*) *rxn_ptr )[0];
  int dbl_index = ( (int*) *rxn_ptr )[1];
  int * spec_ids = &( ( (int*) *rxn_ptr )[2] );
  double * rxn_params = (double*) &( ( (int*) *rxn_ptr )[ dbl_index ] );

  // The actual operations performed depend on the reaction type, but only
  // involve reading from this reaction's parameters in rxn_ptr, reading the
  // value of any number of elements in the state array and adding values to
  // any number of elements in the deriv array.
  for ( int i_spec = 0; i_spec < n_rxn_spec; i_spec++ ) {
    atomicAdd( &( deriv[ spec_ids[ i_spec ] ] ),state[ spec_ids[ i_spec ] ] * rxn_params[ i_spec ] );
	//deriv[spec_ids[i_spec]]+= state[spec_ids[ i_spec]]*rxn_params[i_spec];
  }

  //*rxn_ptr = (void*) &( rxn_params[ n_rxn_spec ] );
}


/** Solve a block of reactions on a GPU
 *
 * @param start_rxn_param Relative starting address of the reaction parameters
 * @param end_rxn_param Relative ending address of the reaction parameters
 * @param dev_state State array on the device
 * @param dev_rxn_param Reaction parameters on the device
 * @param dev_deriv Derivative being calculated on the device
 * @param dev_temp_deriv Derivative for use on the GPU locally
 * param deriv Derivative on the host
 */
__global__ void solveRxnBlock( unsigned int * dev_start_rxn_param,
    unsigned int * dev_end_rxn_param, float * dev_state, void * dev_rxn_param,
    float * dev_deriv )
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int grid_size = gridDim.x * blockDim.x;

  //float shr_dev_deriv[ N_SPEC ];

  __shared__ float shr_dev_deriv[ N_SPEC ];


  //for ( int i_spec = threadIdx.x; i_spec < N_SPEC; i_spec += blockDim.x )
	//shr_dev_deriv[ i_spec ] = dev_deriv[ i_spec ];


  for ( int i_spec = threadIdx.x; i_spec < N_SPEC; i_spec += blockDim.x )
	shr_dev_deriv[ i_spec ] = 0.0;
  for ( int i_spec = index; i_spec < N_SPEC; i_spec += grid_size )
	dev_deriv[ i_spec ] = 0.0;

  __syncthreads();
  /*for ( void * rxn_ptr = (void*) &(((char*) dev_rxn_param)[ //Assings rxn_ptr start for every thread
        dev_start_rxn_param[ index ] ] ) ;
        rxn_ptr <= (void*) &(((char*) dev_rxn_param)[//This condition limit threads, doing it only the needed
        dev_end_rxn_param[ index ] ] ); ) {
*/
  	  void * rxn_ptr = (void*) &(((char*) dev_rxn_param)[dev_start_rxn_param[ index ]] );

 	  if (rxn_ptr){

	  	//devSolveRxn( &rxn_ptr, dev_state, shr_dev_deriv );
  	  devSolveRxn( &rxn_ptr, dev_state, shr_dev_deriv );

  }
  __syncthreads();

  for ( int i_spec = threadIdx.x; i_spec < N_SPEC; i_spec += blockDim.x )// Unnecesary creo. Solo se ejecuta una vez. Funciona porque asigna los indices segun el id del thread
  atomicAdd( &( dev_deriv[ i_spec ] ), shr_dev_deriv[ i_spec ] );

}

/** Solve the mechanism with GPUs
 *
 * @param n_gpu_block Number of GPU blocks to use
 * @param n_gpu_thread Number of GPU threads to use
 * @param state Pointer to the state array
 * @param dev_rxn_param Pointer to the reaction parameter block on the GPUs
 * @param rxn_param Pointer to the reaction parameter block
 * @param dev_deriv Pointer to the derivative array to calculate on the GPUs
 * @param dev_temp_deriv Pointer to a derivative array for use on the GPU only
 * @param deriv Pointer to the derivative array to calculate
 * @return Computational time (s)
 */
float solveWithGpus( int n_gpu_block, int n_gpu_thread, float * dev_state,
    void * dev_rxn_param, void * rxn_param, float * dev_deriv,
    float * deriv, size_t deriv_size )
{
  // Set the max number of GPUs to use
  dim3 dimGrid(n_gpu_block);
  dim3 dimBlock(n_gpu_thread);

  // Reaction parameters and GPU pointers would be copied during initialization
  // so they are not included in the computation time
  size_t size = MAX_N_GPU_BLOCK * MAX_N_GPU_THREAD * sizeof(unsigned int);
  unsigned int *start_rxn_param;
  unsigned int *end_rxn_param;
  HANDLE_ERROR( cudaHostAlloc( (void**) &start_rxn_param, size,
				cudaHostAllocWriteCombined |
					cudaHostAllocMapped ) );
  HANDLE_ERROR( cudaHostAlloc( (void**) &end_rxn_param, size,
				cudaHostAllocWriteCombined |
					cudaHostAllocMapped ) );

  int n_rxn_per_thread = ceil( N_RXN / ( (float) n_gpu_block * n_gpu_thread ) );
  void * rxn_ptr = rxn_param;
  int curr_rxn = 0;
  for ( int i_block = 0; i_block < n_gpu_block; i_block++ ) {
    for ( int i_thread = 0; i_thread < n_gpu_thread; i_thread++ ) {
      start_rxn_param[ i_block * n_gpu_thread + i_thread ] =
          (unsigned int) ( (char*) rxn_ptr - (char*) rxn_param );

      for ( int i_rxn = 0; i_rxn < n_rxn_per_thread && curr_rxn < N_RXN;
          i_rxn++, curr_rxn++) rxn_ptr = skipRxn( rxn_ptr );

      end_rxn_param[ i_block * n_gpu_thread + i_thread ] =
          (unsigned int) ( (char*) rxn_ptr - (char*) rxn_param - 1 );
    }
  }
  unsigned int * dev_start_rxn_param;
  unsigned int * dev_end_rxn_param;
  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_start_rxn_param,
					start_rxn_param, 0 ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_end_rxn_param,
					end_rxn_param, 0 ) );


  clock_t comp_begin = clock();
  // Solve the reactions in blocks on each GPU
  for (int i_repeat = 0; i_repeat < N_REPEAT; i_repeat++ ) {

    //for (int i_spec = 0; i_spec < N_SPEC; i_spec++ ) deriv[ i_spec ] = 0.0;

	  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_rxn_param, rxn_param, 0 ) );

    cudaMemcpy( dev_deriv, deriv, deriv_size, cudaMemcpyHostToDevice ); //0.02 secs more on copy, from 1.3x to 1.1x speedup
    solveRxnBlock<<< dimGrid, dimBlock >>>( dev_start_rxn_param, dev_end_rxn_param,
        dev_state, dev_rxn_param, dev_deriv );
    cudaMemcpy( deriv, dev_deriv, deriv_size, cudaMemcpyDeviceToHost );
  }
  //cudaMemcpy( deriv, dev_deriv, deriv_size, cudaMemcpyDeviceToHost );
  clock_t comp_end = clock();

  cudaFree( dev_start_rxn_param );
  cudaFree( dev_end_rxn_param );

  return ( (float) ( comp_end - comp_begin ) ) / CLOCKS_PER_SEC;
}

/** Solve the simulated mechanism with and without GPUs
 */
int main( int argc, char *argv[] ) {

  // Print device properties
  int deviceCount;
  HANDLE_ERROR( cudaGetDeviceCount(&deviceCount) );
  int device;
  for ( device = 0; device < deviceCount; device++ ) {
  	cudaDeviceProp deviceProp;
	HANDLE_ERROR( cudaGetDeviceProperties( &deviceProp, device ) );
	printf( "Device %d has compute capability %d.%d.",
		device, deviceProp.major, deviceProp.minor );
        if ( deviceProp.integrated == 1 ) printf( " Integrated device." );
        if ( deviceProp.canMapHostMemory == 1 ) printf( " Can map host memory." );
	printf( "\n" );
  }

  // Output results header
  printf( "\nnum_blocks\tnum_threads\tcomp_time" );
//  for ( int i_max_rxn = MIN_RXN_SIZE; i_max_rxn <= MAX_RXN_SIZE; i_max_rxn++ )
//    printf( "\t%d", i_max_rxn );
  printf("\n");

  // Allocate and initialize the state array
  size_t state_size = N_SPEC * sizeof(float);
  float * state;
  float * dev_state;
  HANDLE_ERROR( cudaHostAlloc( (void**) &state,
				state_size,
				cudaHostAllocWriteCombined |
					cudaHostAllocMapped ) );
  for ( int i_spec = 0; i_spec < N_SPEC; i_spec++ )
    state[ i_spec ] = 1.0e-3 * ( rand() % 2000 );
  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_state, state, 0 ) );

  // Get a set of random reactions
  void * rxn_param;
  void * dev_rxn_param;
  newReactionSet( &rxn_param, &dev_rxn_param );

  // Set up a derivative array to solve without GPUs for a reference
  float * deriv_no_gpus = (float*) malloc( N_SPEC * sizeof(float) );
  if ( deriv_no_gpus == NULL ) {
    printf( "\n\nError allocating derivative array without GPUs.\n");
    exit( 1 );
  }

  // Set up a derivative array to solve with GPUs
  size_t deriv_size = N_SPEC * sizeof(float);
  float * deriv;
  float * dev_deriv;
  HANDLE_ERROR( cudaHostAlloc( (void**) &deriv,
				deriv_size,
				cudaHostAllocMapped ) );
  for (int i_spec = 0; i_spec < N_SPEC; i_spec++ ) deriv[ i_spec ] = 0.0;
  HANDLE_ERROR( cudaHostGetDevicePointer( &dev_deriv, deriv, 0 ) );

  // Solve the system without GPUs
  float comp_time = 0.0;
  comp_time = solveNoGpus( state, rxn_param, deriv_no_gpus );

  // Output the non-GPU derivative array
  printDeriv( 0, 0, comp_time, deriv_no_gpus );

  // Solve the system with GPUs
  for ( int i_gpu_block = MAX_N_GPU_BLOCK; i_gpu_block <= MAX_N_GPU_BLOCK; i_gpu_block++ ) {
    for ( int i_gpu_thread = MAX_N_GPU_THREAD; i_gpu_thread <= MAX_N_GPU_THREAD; i_gpu_thread+=30 ) {

      // Solve
      comp_time = solveWithGpus( i_gpu_block, i_gpu_thread, dev_state, dev_rxn_param,
				 rxn_param, dev_deriv, deriv, deriv_size );

      // Output the GPU derivative array
      printDeriv( i_gpu_block, i_gpu_thread, comp_time, deriv );

      // Compare the results
      assertEqual( deriv_no_gpus, deriv );

    }
  }

  HANDLE_ERROR( cudaFreeHost( rxn_param ) );
  HANDLE_ERROR( cudaFreeHost( state ) );
  HANDLE_ERROR( cudaFreeHost( deriv ) );



}


