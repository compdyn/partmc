/* Copyright (C) 2020 Christian Guzman and Guillermo Oyarzun
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Basic GPU functions
 *
 */

#include<iostream>
#include<cuda.h>
#include<cuda_runtime.h>
#include<cuda_runtime_api.h>

//#include "libsolv2.h"

//#include<cublas.h> //todo fix cublas not working
//#include<cublas_v2.h>

using namespace std;

//
//dAthreads
//
// Para reservar memoria Double e Int
extern "C++" void cudaMallocDouble(double* &vector,int size)
{        
	cudaMalloc((void**)&vector,size*sizeof(double));
}

extern "C++" void cudaMallocInt(int* &vector,int size)
{        
	cudaMalloc((void**)&vector,size*sizeof(int));
}

// Para copiar a CPU->GPU Double e Int
extern "C++" void cudaMemcpyDToGpu(double* h_vect,double* d_vect,int size )
{
  cudaMemcpy(d_vect,h_vect,size*sizeof(double),cudaMemcpyHostToDevice);
}

extern "C++" void cudaMemcpyIToGpu(int* h_vect,int* d_vect,int size )
{
		cudaMemcpy(d_vect,h_vect,size*sizeof(int),cudaMemcpyHostToDevice);
}

// Para copiar a GPU->CPU Double e Int
extern "C++" void cudaMemcpyIToCpu(int* h_vect, int* d_vect,int size )
{
		cudaMemcpy(h_vect,d_vect,size*sizeof(int),cudaMemcpyDeviceToHost);
}

extern "C++" void cudaMemcpyDToCpu(double* h_vect, double* d_vect,int size )
{
  cudaMemcpy(h_vect,d_vect,size*sizeof(double),cudaMemcpyDeviceToHost);
}

// Para liberar memoria
extern "C++" void cudaFreeMem(void* vector)
{
	cudaFree(vector);
}

extern "C++" void cudaGetLastErrorC(){
     cudaError_t error;
     error=cudaGetLastError();
     if(error!= cudaSuccess)
     {
       cout<<" ERROR INSIDE A CUDA FUNCTION: "<<error<<" "<<cudaGetErrorString(error)<<endl;
       exit(0);
     }
}


__global__ void cudamatScaleAddI(int nrows, double* dA, int* djA, int* diA, double alpha)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
            int jstart = diA[row];
            int jend   = diA[row+1];
            for(int j=jstart; j<jend; j++)
		{
	            if(djA[j]==row){
			dA[j] = 1.0 + alpha*dA[j];
		    }
		    else{
			dA[j] = alpha*dA[j];
		    }	   
		}
	}
        
}


// OJO: La matrix ya tiene que estar cargada en la GPU para que funcione
// Esta escrito basandose en formato CSR, pero funciona para CSC tambien
// dA  : array de coeficientes de la matriz en la GPU (tama�o nnz)
// djA : array de columnas de la matriz en la GPU ( tama�o nnz)
// diA : array de filas de la matriz en la GPU (tama�o nrows+1)
// alpha : coeficiente para escalar la matriz (ni idea de donde sale:P )
//extern C
extern "C++" void gpu_matScaleAddI(int nrows, double* dA, int* djA, int* diA, double alpha, int blocks, int threads)
{

   blocks = (nrows+threads-1)/threads;
   
   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1);

  cudamatScaleAddI<<<dimGrid,dimBlock>>>(nrows, dA, djA, diA, alpha);
}

// Diagonal precond
__global__ void cudadiagprecond(int nrows, double* dA, int* djA, int* diA, double* ddiag)
{
  int row= threadIdx.x + blockDim.x*blockIdx.x;
  if(row < nrows){
    int jstart=diA[row];
    int jend  =diA[row+1];
    for(int j=jstart;j<jend;j++){
      if(djA[j]==row){
        if(dA[j]!=0.0)
          ddiag[row]= 1.0/dA[j];
        else{
          ddiag[row]= 1.0;
        }
      }
    }
  }

}

extern "C++" void gpu_diagprecond(int nrows, double* dA, int* djA, int* diA, double* ddiag, int blocks, int threads)
{

  blocks = (nrows+threads-1)/threads;

  dim3 dimGrid(blocks,1,1);
  dim3 dimBlock(threads,1,1);

  cudadiagprecond<<<dimGrid,dimBlock>>>(nrows, dA, djA, diA, ddiag);
}

// y = constant
__global__ void cudasetconst(double* dy,double constant,int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dy[row]=constant;
	}
}

extern "C++" void gpu_yequalsconst(double *dy, double constant, int nrows, int blocks, int threads)
{
   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudasetconst<<<dimGrid,dimBlock>>>(dy,constant,nrows);

}


// x=A*b
__global__ void cudaSpmvCSR(double* dx, double* db, int nrows, double* dA, int* djA, int* diA)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
            int jstart = diA[row];
            int jend   = diA[row+1];
            double sum = 0.0;
            for(int j=jstart; j<jend; j++)
		{
			sum+= db[djA[j]]*dA[j];
		}
	    dx[row]=sum;
	}
 
}

__global__ void cudaSpmvCSC(double* dx, double* db, int nrows, double* dA, int* djA, int* diA)
{
	double mult;
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
            int jstart = diA[row];
            int jend   = diA[row+1];
            for(int j=jstart; j<jend; j++)
		{
			mult = db[row]*dA[j];
			atomicAdd(&(dx[djA[j]]),mult);
//			dx[djA[j]]+= db[row]*dA[j];
		}
	}
}

extern "C++" void gpu_spmv(double* dx ,double* db, int nrows, double* dA, int *djA,int *diA,int mattype,int blocks,int  threads)
{
   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 

   if(mattype==0) {
	cudaSpmvCSR<<<dimGrid,dimBlock>>>(dx, db, nrows, dA, djA, diA);		
   }
   else
   {
	cudasetconst<<<dimGrid,dimBlock>>>(dx,0.0, nrows);
	cudaSpmvCSC<<<dimGrid,dimBlock>>>(dx, db, nrows, dA, djA, diA);
   }
}

// y= a*x+ b*y
__global__ void cudaaxpby(double* dy,double* dx, double a, double b, int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dy[row]= a*dx[row] + b*dy[row];
	}
}

extern "C++" void gpu_axpby(double* dy ,double* dx, double a, double b, int nrows, int blocks, int threads)
{

   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudaaxpby<<<dimGrid,dimBlock>>>(dy,dx,a,b,nrows);
}

// y = x
__global__ void cudayequalsx(double* dy,double* dx,int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dy[row]=dx[row];
	}
}

extern "C++" void gpu_yequalsx(double *dy, double* dx, int nrows, int blocks, int threads)
{
   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudayequalsx<<<dimGrid,dimBlock>>>(dy,dx,nrows);

}

__global__ void cudadotxy(double *g_idata1, double *g_idata2, double *g_odata, unsigned int n)
{
  extern __shared__ double sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  double mySum = (i < n) ? g_idata1[i]*g_idata2[i] : 0;

  if (i + blockDim.x < n)
    mySum += g_idata1[i+blockDim.x]*g_idata2[i+blockDim.x];

  sdata[tid] = mySum;
  __syncthreads();

  for (unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s)
      sdata[tid] = mySum = mySum + sdata[tid + s];

    __syncthreads();
  }

  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

extern "C++" double gpu_dotxy(double* vec1, double* vec2, double* h_temp, double* d_temp, int nrows, int blocks,int threads)
{
  dim3 dimGrid(blocks,1,1);
  dim3 dimBlock(threads,1,1);

  cudadotxy<<<dimGrid,dimBlock,threads*sizeof(double)>>>(vec1,vec2,d_temp,nrows);
  cudaMemcpy(h_temp, d_temp, blocks * sizeof(double), cudaMemcpyDeviceToHost);

  double sum=0;
  for(int i=0;i<blocks;i++)
  {
    sum+=h_temp[i];
  }
  return sum;
}

/*
extern "C++" double gpu_dotxy(double *dy, double* dx, int nrows)
{
   double dot=0.0;
   cublasHandle_t hl;
   cublasCreate(&hl);

   cublasDdot(hl,nrows,dy,1,dx,1,&dot);

   cublasDestroy(hl);
   return dot;
}
*/

// z= a*z + x + b*y
__global__ void cudazaxpbypc(double* dz, double* dx,double* dy, double a, double b, int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dz[row]=a*dz[row]  + dx[row] + b*dy[row];
	}
}

extern "C++" void gpu_zaxpbypc(double* dz, double* dx ,double* dy, double a, double b, int nrows, int blocks, int threads)
{

   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudazaxpbypc<<<dimGrid,dimBlock>>>(dz,dx,dy,a,b,nrows);
}

// z= x*y
__global__ void cudamultxy(double* dz, double* dx,double* dy, int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dz[row]=dx[row]*dy[row];
	}
}

extern "C++" void gpu_multxy(double* dz, double* dx ,double* dy, int nrows, int blocks, int threads)
{

   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudamultxy<<<dimGrid,dimBlock>>>(dz,dx,dy,nrows);
}


// z= a*x + b*y
__global__ void cudazaxpby(double* dz, double* dx,double* dy, double a, double b, int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dz[row]=a*dx[row] + b*dy[row];
	}
}

extern "C++" void gpu_zaxpby(double a, double* dx ,double b, double* dy, double* dz, int nrows, int blocks, int threads)
{

   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudazaxpby<<<dimGrid,dimBlock>>>(dz,dx,dy,a,b,nrows);
}

// y= a*x + y
__global__ void cudaaxpy(double* dy,double* dx, double a, int nrows)
{
	int row= threadIdx.x + blockDim.x*blockIdx.x;
   	if(row < nrows){
		dy[row]=a*dx[row] + dy[row];
	}
}

extern "C++" void gpu_axpy(double* dy, double* dx ,double a, int nrows, int blocks, int threads)
{

   dim3 dimGrid(blocks,1,1);
   dim3 dimBlock(threads,1,1); 
   
   cudaaxpy<<<dimGrid,dimBlock>>>(dy,dx,a,nrows);
}

// sqrt(sum ( (x_i*y_i)^2)/n)
__global__ void cudaDVWRMS_Norm(double *g_idata1, double *g_idata2, double *g_odata, unsigned int n)
{
  extern __shared__ double sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  double mySum = (i < n) ? g_idata1[i]*g_idata1[i]*g_idata2[i]*g_idata2[i] : 0;

  if (i + blockDim.x < n)
    mySum += g_idata1[i+blockDim.x]*g_idata1[i+blockDim.x]*g_idata2[i+blockDim.x]*g_idata2[i+blockDim.x];

  sdata[tid] = mySum;
  __syncthreads();

  for (unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s)
      sdata[tid] = mySum = mySum + sdata[tid + s];

    __syncthreads();
  }

  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

extern "C++" double gpu_VWRMS_Norm(int n, double* vec1,double* vec2,double* h_temp,double* d_temp, int blocks,int threads)
{
  dim3 dimGrid(blocks,1,1);
  dim3 dimBlock(threads,1,1);

  cudaDVWRMS_Norm<<<dimGrid,dimBlock,threads*sizeof(double)>>>(vec1,vec2,d_temp,n);

  cudaMemcpy(h_temp, d_temp, blocks * sizeof(double), cudaMemcpyDeviceToHost);

  double sum=0;
  for(int i=0;i<blocks;i++)
  {
    sum+=h_temp[i];
  }
  return sqrt(sum/n);
}

// y=alpha*y
__global__ void cudascaley(double* dy, double a, int nrows)
{
  int row= threadIdx.x + blockDim.x*blockIdx.x;
  if(row < nrows){
    dy[row]=a*dy[row];
  }
}

extern "C++" void gpu_scaley(double* dy, double a, int nrows, int blocks, int threads)
{
  dim3 dimGrid(blocks,1,1);
  dim3 dimBlock(threads,1,1);

  cudascaley<<<dimGrid,dimBlock>>>(dy,a,nrows);
}

