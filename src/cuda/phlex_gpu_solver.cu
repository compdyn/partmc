//Test
extern "C"{
#include "phlex_gpu_solver.h"
#include "rxns_gpu.h"
//#include "phlex_solver.h"
//}

const int N = 16; 
const int blocksize = 16; 

SolverDatagpu *sdgpu;
ModelDatagpu *mdgpu;

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

void solver_new_gpu_cu(SolverDatagpu *sd){

  /*sdgpu = (SolverDatagpu*)malloc(sizeof(SolverDatagpu));
    //SolverData *sd = (SolverData*) malloc(sizeof(SolverData));
  if (sdgpu==NULL) {
    printf("\n\nERROR allocating space for SolverDatagpu\n\n");
    exit(1);
  }*/

  /*
  sdgpu = (SolverDatagpu*)sd;
  sdgpu->curr_J_guess = false;

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

}

void solver_initialize_gpu_cu(void *sd) {



}

void solver_run_gpu_cu(void *sd) {
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
/*
  //TTODO: Put this loop at new_phlex_solver since data matrix don't differ in iterations
  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    //Reaction pointers
    start_rxn_param[i_rxn]=(unsigned int) ( (int*) rxn_ptr2 - (int*) rxn_param );
    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                model_data, deriv_data, (void*) rxn_data, time_step);
        break;
    }
  }
*/
  //solveRxnBlock<<< dimGrid, dimBlock >>>( model_data, deriv_data, (void*) rxn_data, time_step );

  ////cudaMemcpy( dev_deriv, deriv, deriv_size, cudaMemcpyHostToDevice ); //crec que no fa falta si ta link el punter de dades
  //    solveRxnBlock<<< dimGrid, dimBlock >>>( dev_start_rxn_param, dev_end_rxn_param,
  //        dev_state, dev_rxn_param, dev_deriv );
  //cudaMemcpy( deriv, dev_deriv, deriv_size, cudaMemcpyDeviceToHost );


}
/*
__global__ void solveRxnBlock( model_data, deriv_data, (void*) rxn_data, time_step ){


 // Loop through the reactions to determine the Jacobian elements used
  // advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_get_used_jac_elem(
                  model_data, (void*) rxn_data, jac_struct);
        break;
      case RXN_PDFITE_ACTIVITY:
        rxn_data = (int*) rxn_PDFiTE_activity_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
    }
  }



}*/

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
