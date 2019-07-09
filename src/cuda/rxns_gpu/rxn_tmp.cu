


extern "C" {

  //Works calling as a global

#define RXN_ARRHENIUS 1
#include "../rxns_gpu.h"

//TODO: Bug where calling an empty function an memcpy data rises an illegal instruction error.
  //a compiler bug in the CUDA 6.5 and 7.0 release toolkit for compute capability 3.0/3.5 devices

__global__ void rxn_gpu_tmp_arrhenius
          (
          //ModelDatagpu *model_data, double *state,
          //double *deriv, int *rxn_data, double *double_pointer_gpu,
          //double time_step, int n_rxn2

          ModelDatagpu *model_data, double *state, double *deriv,
          double time_step, int deriv_length, int n_rxn,
          int *int_pointer, double *double_pointer,
          unsigned int int_max_size, unsigned int double_max_size
          )
{

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ double deriv_data[MAX_SHARED_MEMORY_BLOCK_DOUBLE];

  if (threadIdx.x < deriv_length){ //This produces seg.fault for some large values seems
    deriv_data[index] = 0.0;
  }

  __syncthreads();

  if (index < n_rxn) {

    int *int_data = (int *) &(((int *) int_pointer)[index]);
    double *float_data = (double *) &(((double *) double_pointer)[index]);

    int rxn_type = int_data[0];
    int *rxn_data = (int *) &(int_data[n_rxn]);


    switch (rxn_type) {
      case RXN_ARRHENIUS :

        rxn_gpu_arrhenius_calc_deriv_contrib(
                model_data, state, deriv_data, (void *) rxn_data, float_data, time_step, deriv_length, n_rxn);


        //int n_rxn=n_rxn2;
        //double *state = model_data->state;//TODO: model_data->state[i] to calculate independent cells simultaneous


/*
        int_data =rxn_data;
        //int *int_data = (int*) rxn_data;
        //double *float_data = double_pointer;

        // Calculate the reaction rate
        double rate = float_data[6*n_rxn];
        for (int i_spec=0; i_spec<int_data[0]; i_spec++) rate *= state[int_data[(2 + i_spec)*n_rxn]-1];

        // Add contributions to the time derivative
        if (rate!=ZERO) {
          int i_dep_var = 0;
          for (int i_spec=0; i_spec<int_data[0]; i_spec++, i_dep_var++) {
            if (int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn] < 0) continue;
            //deriv[DERIV_ID_(i_dep_var)] -= rate;
            atomicAdd(&(deriv_data[int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn]]),-rate);
          }
          for (int i_spec=0; i_spec<int_data[1*n_rxn]; i_spec++, i_dep_var++) {
            if (int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn] < 0) continue;

            // Negative yields are allowed, but prevented from causing negative
            // concentrations that lead to solver failures
            if (-rate*float_data[(7 + i_spec)*n_rxn]*time_step <=state[int_data[(2 + int_data[0] + i_spec)*n_rxn]-1]) {
              //deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
              atomicAdd(&(deriv_data[int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn]]),
                        rate*float_data[(7 + i_spec)*n_rxn]);
            }
          }
        }
*/
        break;
    }
  }

  __syncthreads();

  if (threadIdx.x < deriv_length)
    deriv[index] = deriv_data[index];

}

}