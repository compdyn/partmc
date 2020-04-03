#include "itsolver_gpu.h"

void createSolver(itsolver *bicg)
{
  //Init variables ("public")
  int nrows = bicg->nrows;
  int blocks = bicg->blocks;

  //Auxiliary vectors ("private")
  double ** dr0 = &bicg->dr0;
  double ** dr0h = &bicg->dr0h;
  double ** dn0 = &bicg->dn0;
  double ** dp0 = &bicg->dp0;
  double ** dt = &bicg->dt;
  double ** ds = &bicg->ds;
  double ** dAx2 = &bicg->dAx2;
  double ** dy = &bicg->dy;
  double ** dz = &bicg->dz;
  double ** daux = &bicg->daux;
  double ** ddiag = &bicg->ddiag;

  //Allocate
  cudaMalloc(dr0,nrows*sizeof(double));
  cudaMalloc(dr0h,nrows*sizeof(double));
  cudaMalloc(dn0,nrows*sizeof(double));
  cudaMalloc(dp0,nrows*sizeof(double));
  cudaMalloc(dt,nrows*sizeof(double));
  cudaMalloc(ds,nrows*sizeof(double));
  cudaMalloc(dAx2,nrows*sizeof(double));
  cudaMalloc(dy,nrows*sizeof(double));
  cudaMalloc(dz,nrows*sizeof(double));
  cudaMalloc(ddiag,nrows*sizeof(double));
  cudaMalloc(daux,nrows*sizeof(double));
  bicg->aux=(double*)malloc(sizeof(double)*blocks);

}

__global__
void cudaSolveGPU(
        double *dA, int *djA, int *diA, double *dx, double *dtempv, //Input data
        int nrows, int blocks, int threads, int maxIt, int mattype,
        int n_cells, double tolmax, double *ddiag, //Init variables
        double *dr0, double *dr0h, double *dn0, double *dp0,
        double *dt, double *ds, double *dAx2, double *dy, double *dz,
        double *daux, // Auxiliary vectors
        double *aux_params
        //double *alpha, double *rho0, double* omega0, double *beta,
        //double *rho1, double *temp1, double *temp2 //Auxiliary parameters
        )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;
  //alpha  = 1.0;
  //rho0   = 1.0;
  //omega0 = 1.0;

  //int it=0;
  //do
  //{

  alpha = aux_params[0];
  rho0 = aux_params[1];
  omega0 = aux_params[2];
  beta = aux_params[3];
  rho1 = aux_params[4];
  temp1 = aux_params[5];
  temp2 = aux_params[6];

  //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);
  cudaDevicedotxy(dr0, dr0h, daux, nrows);
  //__syncthreads();
  cudaDevicereducey(daux, blocks);
  //__syncthreads();
  //Note more computations have been done on other blocks, but we only need block 0
  rho1=daux[blockIdx.x];//All threads access first value of daux, so it should be ok (plot twist: is not okay)
  beta=(rho1/rho0)*(alpha/omega0);


  //hum in some way WE NEED to use the same value beta value :/

  //gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c
  cudaDevicezaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows);   //z = ax + by + c
  //if (tid == 0)

  //Fails because there are no command to SYNB BLOCKS, so maybe one block update global memory
  //and other block tries to access global memory that is not updated because no sync
  //todo for some reason, threads has different beta values producing diff values in next function
  if (blockIdx.x == 0)//why block 1 is different that block 0 in daux[0] if daux[0] should be global
  {
    aux_params[0]=alpha;
    aux_params[1]=rho0;
    aux_params[2]=omega0;
    aux_params[3]=beta;//0.01;
    aux_params[4]=rho1;//rho1
    aux_params[5]=temp1;
    aux_params[6]=temp2;
  }
  /*

  //gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag
  cudaDevicemultxy(dy,ddiag,dp0,nrows);

  //gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y
  cudaDevicesetconst(dn0,0.0,nrows);
  cudaDeviceSpmvCSC(dn0,dy,nrows,dA,djA,diA);



  //temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);
  cudaDevicedotxy(dr0h, dn0, daux, nrows);
  //__syncthreads();
  cudaDevicereducey(daux, blocks);
  temp1=daux[blockIdx.x];
  alpha=rho1/temp1;



  //gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);
  cudaDevicezaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows);

  //gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s
  cudaDevicemultxy(dz,ddiag,ds,nrows); // precond z=diag*s


  //gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);
  cudaDevicesetconst(dt,0.0,nrows);
  cudaDeviceSpmvCSC(dt,dz,nrows,dA,djA,diA);

  //gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);
  cudaDevicemultxy(dAx2,ddiag,dt,nrows);

  //temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
  cudaDevicedotxy(dz, dAx2, daux, nrows);
  //__syncthreads();
  cudaDevicereducey(daux, blocks);
  temp1=daux[blockIdx.x];

  //temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
  cudaDevicedotxy(dAx2, dAx2, daux, nrows);
  //__syncthreads();
  cudaDevicereducey(daux, blocks);
  temp2=daux[blockIdx.x];
  omega0=temp1/temp2;

  //gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x
  cudaDeviceaxpy(dx,dy,alpha,nrows); // x=alpha*y +x

  //gpu_axpy(dx,dz,omega0,nrows,blocks,threads);
  cudaDeviceaxpy(dx,dz,omega0,nrows);

  //gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);
  cudaDevicezaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows);

  //temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
  cudaDevicedotxy(dr0, dr0, daux, nrows);
  //__syncthreads();
  cudaDevicereducey(daux, blocks);
  temp1=daux[blockIdx.x];
  temp1=sqrt(temp1);

  rho0=rho1;
*/
  //  it++;
  //} while(it<maxIt || temp1>tolmax);
}

//todo quit syncs with cpu/gpu
//Biconjugate gradient
void solveGPU(itsolver *bicg, double *dA, int *djA, int *diA, double *dx, double *dtempv)
{
  //Init variables ("public")
  int nrows_total = bicg->nrows;
  int blocks = bicg->blocks;
  int threads = bicg->threads;
  int maxIt = bicg->maxIt;
  int mattype = bicg->mattype;
  int n_cells = bicg->n_cells;
  double tolmax = bicg->tolmax;
  double *ddiag = bicg->ddiag;

  // Auxiliary vectors ("private")
  double *dr0 = bicg->dr0;
  double *dr0h = bicg->dr0h;
  double *dn0 = bicg->dn0;
  double *dp0 = bicg->dp0;
  double *dt = bicg->dt;
  double *ds = bicg->ds;
  double *dAx2 = bicg->dAx2;
  double *dy = bicg->dy;
  double *dz = bicg->dz;
  double *aux = bicg->aux;
  double *daux = bicg->daux;

  //Function private variables
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;
  int nrows = nrows_total;

  gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x

  gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0 //y=ax+by

  gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0

  gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0 //memset???
  gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  //Todo apply multicells (and remember threads should be power of two for dotxy reduce)
  //Each block will compute only a cell/group of cells
  /*
  int n_rows_cell = nrows_total/n_cells;
  int max_threads = bicg->threads;
  int cells_per_block = (max_threads+n_rows_cell-1)/n_rows_cell;
  //threads = n_rows_cell*cells_per_block;
  //blocks = (nrows_total+threads-1)/threads;
   */

  //dim3 dimGrid(blocks,1,1);
  //dim3 dimBlock(threads,1,1);

  int n_aux_params=7;
  double *aux_params;
  aux_params=(double*)malloc(n_aux_params*sizeof(double));
  double *daux_params;
  cudaMalloc(&daux_params,n_aux_params*sizeof(double));
  //cudaMemcpy(bicg->djA,bicg->jA,7*sizeof(double),cudaMemcpyHostToDevice);

  int redsize= sqrt(blocks) +1;
  redsize=pow(2,redsize);

  for(int it=0;it<maxIt;it++){

    aux_params[0]=alpha;
    aux_params[1]=rho0;
    aux_params[2]=omega0;
    aux_params[3]=beta;
    aux_params[4]=rho1;
    aux_params[5]=temp1;
    aux_params[6]=temp2;

    cudaMemcpy(daux_params,aux_params,n_aux_params*sizeof(double),cudaMemcpyHostToDevice);
//dimGrid,dimBlock
    cudaSolveGPU<<<blocks,threads,threads*sizeof(double)>>>
                                    (dA,djA,diA,dx,dtempv,nrows,redsize,threads,maxIt,mattype,n_cells,
                                        tolmax,ddiag,dr0,dr0h,dn0,dp0,dt,ds,dAx2,dy,dz,daux,
                                        daux_params
                                        //&alpha,&rho0,&omega0,&beta,&rho1,&temp1,&temp2
                                        );
    cudaDeviceSynchronize();
    cudaMemcpy(aux_params,daux_params,n_aux_params*sizeof(double),cudaMemcpyDeviceToHost);

    alpha = aux_params[0];
    rho0 = aux_params[1];
    omega0 = aux_params[2];
    beta = aux_params[3];
    rho1 = aux_params[4];
    temp1 = aux_params[5];
    temp2 = aux_params[6];



    /*

    rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);//rho1 =<r0,r0h>

    //  cout<<"rho1 "<<rho1<<endl;

    beta=(rho1/rho0)*(alpha/omega0);

    //    cout<<"rho1 "<<rho1<<" beta "<<beta<<endl;
*/
    //printf("rho1 %f", rho1);
    //printf("beta %f", beta);
    gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c

    gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag

    gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y

    temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);

    alpha=rho1/temp1;

    //       cout<<"temp1 "<<temp1<<" alpha "<<alpha<<endl;

    printf("rho1 %f", rho1);
    //printf("alpha %f", alpha);
    gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);

    gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s

    gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);

    gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);

    temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);

    omega0= temp1/temp2;

    gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x

    gpu_axpy(dx,dz,omega0,nrows,blocks,threads);

    gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);



    temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=sqrt(temp1);

    //cout<<it<<": "<<temp1<<endl;

    if(temp1<tolmax){
      break;
    }
    rho0=rho1;
  }

  cudaFreeMem(daux_params);

}

//todo free
//void free_gpu(itsolver *bicg)

 /*
void setUpSolver(itsolver *bicg, double reltol, double *ewt, int tnrows,int tnnz,double *tA, int *tjA, int *tiA, int tmattype, int qmax, double *dACamp, double *dftempCamp);
{

  bicg.tolmax=reltol;

}
*/