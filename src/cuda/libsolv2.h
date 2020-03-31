/* Copyright (C) 2020 Christian Guzman and Guillermo Oyarzun
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Basic GPU functions
 *
 */

#include<iostream>
#include<cuda.h>
//#include"/usr/local/cuda-10.1/include/cuda_runtime.h"
//#include"/usr/local/cuda-10.1/include/cuda_runtime_api.h"

//#include<cublas.h>

extern "C++" void cudaMallocDouble(double* &vector,int size);
extern "C++" void cudaMallocInt(int* &vector,int size);
extern "C++" void cudaMemcpyDToGpu(double* h_vect,double* d_vect,int size );
extern "C++" void cudaMemcpyIToGpu(int* h_vect,int* d_vect,int size );
extern "C++" void cudaMemcpyIToCpu(int* h_vect, int* d_vect,int size );
extern "C++" void cudaMemcpyDToCpu(double* h_vect, double* d_vect,int size );
extern "C++" void cudaFreeMem(void* vector);
extern "C++" void cudaGetLastErrorC();
extern "C++" void gpu_matScaleAddI(int nrows, double* dA, int* djA, int* diA, double alpha, int blocks, int threads);
extern "C++" void gpu_diagprecond(int nrows, double* dA, int* djA, int* diA, double* ddiag, int blocks, int threads);
extern "C++" void gpu_yequalsconst(double *dy, double constant, int nrows, int blocks, int threads);
extern "C++" void gpu_spmv(double* dx ,double* b, int nrows, double* dA, int *djA,int *diA,int mattype,int blocks,int  threads);
extern "C++" void gpu_axpby(double* dy ,double* dx, double a, double b, int nrows, int blocks, int threads);
extern "C++" void gpu_yequalsx(double *dy, double* dx, int nrows, int blocks, int threads);
extern "C++" double gpu_VWRMS_Norm(int n, double* vec1,double* vec2,double* h_temp,double* d_temp, int blocks,int threads);
extern "C++" double gpu_dotxy(double* vec1, double* vec2, double* h_temp, double* d_temp, int nrows, int blocks,int threads);
//extern "C++" double gpu_dotxy(double *dy, double* dx, int nrows);
extern "C++" void gpu_zaxpbypc(double* dz, double* dx ,double* dy, double a, double b, int nrows, int blocks, int threads);
extern "C++" void gpu_multxy(double* dz, double* dx ,double* dy, int nrows, int blocks, int threads);
//extern "C++" void gpu_zaxpby(double* dz, double* dx ,double* dy, double a, double b, int nrows, int blocks, int threads);
extern "C++" void gpu_zaxpby(double a, double* dx ,double b, double* dy, double* dz, int nrows, int blocks, int threads);
extern "C++" void gpu_axpy(double* dy, double* dx ,double a, int nrows, int blocks, int threads);
extern "C++" double gpu_VWRMS_Norm(int n, double* vec1,double* vec2,double* h_temp,double* d_temp, int blocks,int threads);
extern "C++" void gpu_scaley(double* dy, double a, int nrows, int blocks, int threads);

