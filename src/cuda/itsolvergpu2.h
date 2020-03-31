/* Copyright (C) 2020 Christian Guzman and Guillermo Oyarzun
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Iterative GPU solver
 *
 */

#ifndef ITSOLVERGPU_H
#define ITSOLVERGPU_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include"libsolv2.h"

extern "C" {
#include "../camp_solver.h"
}


using namespace std;

void createSolver(itsolver2 *bicg)
{
  //Init variables ("public")
  int nrows = bicg->nrows;
  int blocks = bicg->blocks;
  int threads = bicg->threads;
  int maxIt = bicg->maxIt;
  int mattype = bicg->mattype;
  double tolmax = bicg->tolmax;

  //Intermediate variables ("private")
  double ** dr0 = &bicg->dr0;
  double ** dr0h = &bicg->dr0h;
  double ** dn0 = &bicg->dn0;
  double ** dp0 = &bicg->dp0;
  double ** dt = &bicg->dt;
  double ** ds = &bicg->ds;
  double ** dAx2 = &bicg->dAx2;
  double ** dy = &bicg->dy;
  double ** dz = &bicg->dz;
  double ** aux = &bicg->aux;
  double ** daux = &bicg->daux;
  double ** ddiag = &bicg->ddiag;

  //Allocate
  cudaMalloc(dr0,bicg->nrows*sizeof(double));
  cudaMalloc(dr0h,bicg->nrows*sizeof(double));
  cudaMalloc(dn0,bicg->nrows*sizeof(double));
  cudaMalloc(dp0,bicg->nrows*sizeof(double));
  cudaMalloc(dt,bicg->nrows*sizeof(double));
  cudaMalloc(ds,bicg->nrows*sizeof(double));
  cudaMalloc(dAx2,bicg->nrows*sizeof(double));
  cudaMalloc(dy,bicg->nrows*sizeof(double));
  cudaMalloc(dz,bicg->nrows*sizeof(double));
  cudaMalloc(ddiag,bicg->nrows*sizeof(double));
  cudaMalloc(daux,bicg->nrows*sizeof(double));

  bicg->aux=(double*)malloc(sizeof(double)*bicg->blocks);

}

//Biconjugate gradient
void solveGPU2(itsolver2 *bicg, double *dA, int *djA, int *diA, double *dx, double *dtempv)
{
  //Init variables ("public")
  int nrows = bicg->nrows;
  int blocks = bicg->blocks;
  int threads = bicg->threads;
  int maxIt = bicg->maxIt;
  int mattype = bicg->mattype;
  double tolmax = bicg->tolmax;
  double * ddiag = bicg->ddiag;

  //Intermediate variables ("private")
  double * dr0 = bicg->dr0;
  double * dr0h = bicg->dr0h;
  double * dn0 = bicg->dn0;
  double * dp0 = bicg->dp0;
  double * dt = bicg->dt;
  double * ds = bicg->ds;
  //double * dAx;
  double * dAx2 = bicg->dAx2;
  double * dy = bicg->dy;
  double * dz = bicg->dz;
  //double * diag = bicg->diag;
  double * aux = bicg->aux;
  double * daux = bicg->daux;

  //Function private variables
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;

  gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x

  gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0        //y=ax+by

  gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0
  gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0

  for(int it=0;it<maxIt;it++){

    //rho1=gpu_dotxy(dr0,dr0h,nrows); //rho1 =<r0,r0h>
    rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);

    //  cout<<"rho1 "<<rho1<<endl;
    beta=(rho1/rho0)*(alpha/omega0);

    //    cout<<"rho1 "<<rho1<<" beta "<<beta<<endl;

    gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c

    gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag

    gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y

    //temp1=gpu_dotxy(dr0h,dn0,nrows);
    temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);

    alpha=rho1/temp1;

    //       cout<<"temp1 "<<temp1<<" alpha "<<alpha<<endl;

    gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);
    //gpu_zaxpby(ds,dr0,dn0,1.0,-1.0*alpha,nrows,blocks,threads);

    gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s

    gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);

    gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);

    //temp1=gpu_dotxy(dz,dAx2,nrows);
    temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    //temp2=gpu_dotxy(dAx2,dAx2,nrows);
    temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    omega0= temp1/temp2;

    gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x

    gpu_axpy(dx,dz,omega0,nrows,blocks,threads);

    gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);
    //gpu_zaxpby(dr0,ds,dt,1.0,-1.0*omega0,nrows,blocks,threads);

    //temp1=gpu_dotxy(dr0,dr0,nrows);
    temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=sqrt(temp1);

    //cout<<it<<": "<<temp1<<endl;

    if(temp1<tolmax){
      break;
    }
    rho0=rho1;
  }

  //cudaMemcpyDToCpu(x,dx,nrows);
}

/*
class itsolver{

   public:
	   itsolver(){};
	   itsolver(int, int, double*, int*, int*,int);
	   ~itsolver(){};

      // matrix data
      double* A;
      int*    jA;
      int*    iA;
      int     nrows;
      int     nnz;

      //GPU pointers
      double* dA;
      int*    djA;
      int*    diA;
      double * dx;
      double* ddiag;
      double* aux;
      double* daux;

      // Allocate ewt, acor, tempv, ftemp
      double* dewt;
      double* dacor;
      double* dacor_init;
      double* dtempv;
      double* dftemp;
      double* dzn;
      double* dcv_y;
      int threads,blocks;

      // Auxiliary scalars
      double tolmax;
      int maxIt;
      int mattype;

      int totits;
      double totres;
      // Auxiliary vectors
      double * dr0;
      double * dr0h;
      double * dn0;
      double * dp0;
      double * dt;
      double * ds;
      double * dAx;
      double * dAx2;
      double * dy;
      double * dz;
      double * diag;

      //subroutines
      void setUpSolver(double, double *, int, int, double *, int *, int *, int, int, double *, double *);
      void setUpGPU();
      void solveGPU();

};


itsolver::itsolver(int tnrows,int tnnz,double *tA, int *tjA, int *tiA, int tmattype)
{
}

void itsolver::setUpSolver(double reltol, double *ewt, int tnrows,int tnnz,double *tA, int *tjA, int *tiA, int tmattype, int qmax, double *dACamp, double *dftempCamp)
{
  //Linking Matrix data, later this data must be allocated in GPU
  nrows=tnrows;
  nnz=tnnz;
  mattype=tmattype;
  A=tA;
  jA=tjA;
  iA=tiA;

  int device=0;//Selected GPU
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  threads=prop.maxThreadsPerBlock;//1024; //128 set at max gpu
  blocks=(nrows+threads-1)/threads;

  // allocating matrix data to the GPU
  //cudaMallocDouble(dA,nnz);
  dA = dACamp;//set itsolver gpu pointer to jac pointer initialized at camp

  cudaMallocInt(djA,nnz);
  cudaMallocInt(diA,(nrows+1));

  // moving matrix data to the GPU

  //cudaMemcpyDToGpu(A,dA,nnz);

  cudaMemcpyIToGpu(jA,djA,nnz);
  cudaMemcpyIToGpu(iA,diA,(nrows+1));

  cudaMallocDouble(dr0,nrows);
  cudaMallocDouble(dr0h,nrows);
  cudaMallocDouble(dn0,nrows);
  cudaMallocDouble(dp0,nrows);
  cudaMallocDouble(dt,nrows);
  cudaMallocDouble(ds,nrows);
  cudaMallocDouble(dAx,nrows);
  cudaMallocDouble(dAx2,nrows);
  cudaMallocDouble(dy,nrows);
  cudaMallocDouble(dz,nrows);
  cudaMallocDouble(ddiag,nrows);
  cudaMallocDouble(dx,nrows);

  aux=(double*)malloc(sizeof(double)*blocks);
  cudaMallocDouble(daux,nrows);

  // Allocate ewt, acor, tempv, ftemp
  cudaMallocDouble(dewt,nrows);
  cudaMallocDouble(dacor,nrows);
  cudaMallocDouble(dtempv,nrows);
  //cudaMallocDouble(dacor_init,nrows);

  //cudaMallocDouble(dftemp,nrows);
  dftemp=dftempCamp;

  cudaMallocDouble(dzn,nrows*qmax);
  cudaMallocDouble(dcv_y,nrows);//if anything, set to zero to avoid seg fault errors

  cudaMemcpyDToGpu(ewt,dewt,nnz);
  cudaMemcpyDToGpu(ewt,dacor,nnz);
  //cudaMemcpyDToGpu(ewt,dtempv,nnz);
  cudaMemcpyDToGpu(ewt,dftemp,nnz);

  cudaMemset(dx, 0.0, nrows*sizeof(double));

  diag=(double*)malloc(sizeof(double)*nrows);

  maxIt=100;
  tolmax=reltol; //Set to CAMP selected accuracy (1e-8) //1e-10;//1e-6
}

void itsolver::setUpGPU()
{
  gpu_diagprecond(nrows,dA,djA,diA,ddiag,blocks,threads);
}

//Biconjugate gradient
void itsolver::solveGPU()
{
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;

  gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x

  gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0        //y=ax+by

  gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0
  gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0

  for(int it=0;it<maxIt;it++)
  {

    //rho1=gpu_dotxy(dr0,dr0h,nrows); //rho1 =<r0,r0h>
    rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);

    //  cout<<"rho1 "<<rho1<<endl;
    beta=(rho1/rho0)*(alpha/omega0);

    //    cout<<"rho1 "<<rho1<<" beta "<<beta<<endl;

    gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c

    gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag

    gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y

    //temp1=gpu_dotxy(dr0h,dn0,nrows);
    temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);

    alpha=rho1/temp1;

    //       cout<<"temp1 "<<temp1<<" alpha "<<alpha<<endl;

    gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);
    //gpu_zaxpby(ds,dr0,dn0,1.0,-1.0*alpha,nrows,blocks,threads);

    gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s

    gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);

    gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);

    //temp1=gpu_dotxy(dz,dAx2,nrows);
    temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    //temp2=gpu_dotxy(dAx2,dAx2,nrows);
    temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    omega0= temp1/temp2;

    gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x

    gpu_axpy(dx,dz,omega0,nrows,blocks,threads);

    gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);
    //gpu_zaxpby(dr0,ds,dt,1.0,-1.0*omega0,nrows,blocks,threads);

    //temp1=gpu_dotxy(dr0,dr0,nrows);
    temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=sqrt(temp1);

    //cout<<it<<": "<<temp1<<endl;

    if(temp1<tolmax){
      totits=it;
      totres=temp1;
      break;
    }
    rho0=rho1;
  }

  //cudaMemcpyDToCpu(x,dx,nrows);
}

*/

//void createSolver(itsolver2 *bicg);
//void solveGPU2(itsolver2 *bicg, double *dA, int *djA, int *diA, double *dx, double *dtempv);

#endif