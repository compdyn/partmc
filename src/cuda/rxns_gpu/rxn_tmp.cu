


extern "C" {


#include "../rxns_gpu.h"

//TODO: Bug where calling an empty function an memcpy data rises an illegal instruction error.
  //a compiler bug in the CUDA 6.5 and 7.0 release toolkit for compute capability 3.0/3.5 devices

__device__ void rxn_gpu_tmp_arrhenius(ModelDatagpu *model_data,
                                      double *deriv, int *rxn_data, double *double_pointer_gpu,
                                      double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  double *state = model_data->state;//TODO: model_data->state[i] to calculate independent domains simultaneous
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
        if (-rate*float_data[(7 + i_spec)*n_rxn]*time_step <=state[int_data[(2 + int_data[0] + i_spec)*n_rxn]-1]) {
          //deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
          atomicAdd(&(deriv[int_data[(2 + int_data[0] + int_data[1*n_rxn] + i_dep_var)*n_rxn]]),
                  rate*float_data[(7 + i_spec)*n_rxn]);
        }
    }
  }

}

}