//
// Created by Christian on 01/04/2020.
//

#ifndef CAMPGPU_CUDA_STRUCTS_H
#define CAMPGPU_CUDA_STRUCTS_H

typedef struct
{
  //Init variables ("public")
  int threads,blocks;
  int maxIt;
  int mattype;
  int nrows;
  int nnz;
  int n_cells;
  double tolmax;
  double* ddiag;

  // Intermediate variables ("private")
  double * dr0;
  double * dr0h;
  double * dn0;
  double * dp0;
  double * dt;
  double * ds;
  double * dAx2;
  double * dy;
  double * dz;

  // Matrix data ("input")
  double* A;
  int*    jA;
  int*    iA;

  //GPU pointers ("input")
  double* dA;
  int*    djA;
  int*    diA;
  double* dx;
  double* aux;
  double* daux;

  // ODE solver variables
  double* dewt;
  double* dacor;
  double* dacor_init;
  double* dtempv;
  double* dftemp;
  double* dzn;
  double* dcv_y;

#ifdef PMC_DEBUG_GPU
  int counterSendInit;
  int counterMatScaleAddI;
  int counterMatScaleAddISendA;
  int counterMatCopy;
  int counterNewtonIt;
  int counterLinSolSetup;
  int counterLinSolSolve;
  int countercvStep;
  int counterprecvStep;
  int counterDerivNewton;
  int counterBiConjGrad;
  int counterDerivSolve;
  int counterJac;

  double timeNewtonSendInit;
  double timeMatScaleAddI;
  double timeMatScaleAddISendA;
  double timeMatCopy;
  double timeNewtonIt;
  double timeLinSolSetup;
  double timeLinSolSolve;
  double timecvStep;
  double timeprecvStep;
  double timeDerivNewton;
  double timeBiConjGrad;
  double timeDerivSolve;
  double timeJac;

  cudaEvent_t startDerivNewton;
  cudaEvent_t startDerivSolve;
  cudaEvent_t startLinSolSetup;
  cudaEvent_t startLinSolSolve;
  cudaEvent_t startNewtonIt;
  cudaEvent_t startcvStep;
  cudaEvent_t startBiConjGrad;
  cudaEvent_t startJac;

  cudaEvent_t stopDerivNewton;
  cudaEvent_t stopDerivSolve;
  cudaEvent_t stopLinSolSetup;
  cudaEvent_t stopLinSolSolve;
  cudaEvent_t stopNewtonIt;
  cudaEvent_t stopcvStep;
  cudaEvent_t stopBiConjGrad;
  cudaEvent_t stopJac;



#endif

} itsolver;


#endif //CAMPGPU_CUDA_STRUCTS_H
