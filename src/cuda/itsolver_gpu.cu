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
//todo instead sending all in one kernel, divide in 2 or 4 kernels with streams and check if
//cuda reassigns better the resources
//todo put comments under some DEBUG flag (Like DEBUG_itsolver or something like this, it even
//doesn't need to be open for the user in CMAKE, it could be a defined value in the same file or something like that (ask matt)
//todo profiling del dot y ver cuanta occupancy me esta dando de shared memory porque me limita
//el numero de bloques que se ejecutan a la vez(solo se ejecutan a la vez en toda la function
// los bloques que "quepan" con la shared memory available: solution use cudastreams and launch instead
//of only 1 kernel use 2 or 4 to cubrir huecos (de memoria y eso), y tmb reducir la shared
//con una implementacion hibrida del dotxy

//*se puede probar la funcion en la cpu a ver que hace
__global__
void cudaSolveGPU(
        double *dA, int *djA, int *diA, double *dx, double *dtempv, //Input data
        int nrows, int blocks, int n_shr_empty, int maxIt, int mattype,
        int n_cells, double tolmax, double *ddiag, //Init variables
        double *dr0, double *dr0h, double *dn0, double *dp0,
        double *dt, double *ds, double *dAx2, double *dy, double *dz,
        double *daux // Auxiliary vectors
        //,double *aux_params
        //double *alpha, double *rho0, double* omega0, double *beta,
        //double *rho1, double *temp1, double *temp2 //Auxiliary parameters
)
{
  int tid = threadIdx.x;
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;

  //gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x
  //cudaDevicesetconst(dr0, 0.0, nrows);
  //cudaDeviceSpmvCSC(dr0,dx,nrows,dA,djA,diA);

  //gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0 //y=ax+by
  cudaDeviceaxpby(dr0,dtempv,1.0,-1.0,nrows);

  //gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0
  cudaDeviceyequalsx(dr0h,dr0,nrows);

  //gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0 //memset???
  //gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0
  cudaDevicesetconst(dn0, 0.0, nrows);
  cudaDevicesetconst(dp0, 0.0, nrows);

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  int it=0;
  do
  {

    /*alpha = aux_params[0];
    rho0 = aux_params[1];
    omega0 = aux_params[2];
    beta = aux_params[3];
    rho1 = aux_params[4];
    temp1 = aux_params[5];
    temp2 = aux_params[6];*/

    //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);
    __syncthreads();
    cudaDevicedotxy(dr0, dr0h, &rho1, nrows, n_shr_empty);//&rho1
    __syncthreads();//necessary to reduce accuracy error
    beta = (rho1 / rho0) * (alpha / omega0);

    //gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c
    cudaDevicezaxpbypc(dp0, dr0, dn0, beta, -1.0 * omega0 * beta, nrows);   //z = ax + by + c

    //gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag
    cudaDevicemultxy(dy, ddiag, dp0, nrows);

    //gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y
    cudaDevicesetconst(dn0, 0.0, nrows);
    cudaDeviceSpmvCSC(dn0, dy, nrows, dA, djA, diA);

    //temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);
    cudaDevicedotxy(dr0h, dn0, &temp1, nrows, n_shr_empty);
    __syncthreads();
    alpha = rho1 / temp1;

    //gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);
    cudaDevicezaxpby(1.0, dr0, -1.0 * alpha, dn0, ds, nrows);

    //gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s
    cudaDevicemultxy(dz, ddiag, ds, nrows); // precond z=diag*s

    //gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);
    cudaDevicesetconst(dt, 0.0, nrows);
    cudaDeviceSpmvCSC(dt, dz, nrows, dA, djA, diA);

    //gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);
    cudaDevicemultxy(dAx2, ddiag, dt, nrows);

    //temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    cudaDevicedotxy(dz, dAx2, &temp1, nrows, n_shr_empty);
    __syncthreads();

    //temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    cudaDevicedotxy(dAx2, dAx2, &temp2, nrows, n_shr_empty);
    __syncthreads();
    omega0 = temp1 / temp2;
    __syncthreads();
    //gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x
    cudaDeviceaxpy(dx, dy, alpha, nrows); // x=alpha*y +x
    __syncthreads();

    //gpu_axpy(dx,dz,omega0,nrows,blocks,threads);
    cudaDeviceaxpy(dx, dz, omega0, nrows);

    //gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);
    cudaDevicezaxpby(1.0, ds, -1.0 * omega0, dt, dr0, nrows);

    //temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    cudaDevicedotxy(dr0, dr0, &temp1, nrows, n_shr_empty);
    __syncthreads();
    temp1 = sqrt(temp1);

    rho0 = rho1;
/**/
    __syncthreads();
    /**/

    it++;
  } while(it<maxIt && temp1>tolmax);//while(it<maxIt && temp1>tolmax);//while(0);

  /*
  if (id == 0) //return aux variables if debugging
  {
    aux_params[0]=alpha;
    aux_params[1]=rho0;
    aux_params[2]=omega0;
    aux_params[3]=beta;//0.01;
    aux_params[4]=rho1;//rho1
    aux_params[5]=temp1;
    aux_params[6]=temp2;
  }
*/

}

//solveGPU_multi: Each block will compute only a cell/group of cells
//Algorithm: Biconjugate gradient
void solveGPU_multi(itsolver *bicg, double *dA, int *djA, int *diA, double *dx, double *dtempv)
{
  //Init variables ("public")
  int nrows = bicg->nrows;
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
//todo eliminate atomicadd in spmv through using CSR or something like that
  gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x
/*
  gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0 //y=ax+by

  gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0

  gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0 //memset???
  gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;
*/
  /*int n_aux_params=7;
  double *aux_params;
  aux_params=(double*)malloc(n_aux_params*sizeof(double));
  double *daux_params;
  cudaMalloc(&daux_params,n_aux_params*sizeof(double));*/
  //cudaMemcpy(bicg->djA,bicg->jA,7*sizeof(double),cudaMemcpyHostToDevice);

  int max_threads = bicg->threads;
  int size_cell = nrows/n_cells;
  int n_shr_empty = max_threads%size_cell; //+ size_cell*11 to simulate more shr_empty
  int active_threads = max_threads - n_shr_empty; //last multiple of size_cell before max_threads

  threads = bicg->threads;//active_threads;//bicg->threads;
  blocks = (nrows+active_threads-1)/active_threads; //blocks counting some threads are not working

  /*aux_params[0] = alpha;
  aux_params[1] = rho0;
  aux_params[2] = omega0;
  aux_params[3] = beta;
  aux_params[4] = rho1;
  aux_params[5] = temp1;
  aux_params[6] = temp2;

  cudaMemcpy(daux_params, aux_params, n_aux_params * sizeof(double), cudaMemcpyHostToDevice);*/
  cudaSolveGPU << < blocks, active_threads, bicg->threads * sizeof(double) >> >//threads
                                            (dA, djA, diA, dx, dtempv, nrows, blocks, n_shr_empty, maxIt, mattype, n_cells,
                                                    tolmax, ddiag, dr0, dr0h, dn0, dp0, dt, ds, dAx2, dy, dz, daux
                                                    //,daux_params
                                            );
  /*cudaDeviceSynchronize();
  cudaMemcpy(aux_params, daux_params, n_aux_params * sizeof(double), cudaMemcpyDeviceToHost);

  alpha = aux_params[0];
  rho0 = aux_params[1];
  omega0 = aux_params[2];
  beta = aux_params[3];
  rho1 = aux_params[4];
  temp1 = aux_params[5];
  temp2 = aux_params[6];*/
  //printf("temp1 %-le", temp1);
  //printf("rho1 %f", rho1);

  //cudaFreeMem(daux_params);

}

//Algorithm: Biconjugate gradient
void solveGPU(itsolver *bicg, double *dA, int *djA, int *diA, double *dx, double *dtempv)
{
  //Init variables ("public")
  int nrows = bicg->nrows;
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

  gpu_spmv(dr0,dx,nrows,dA,djA,diA,mattype,blocks,threads);  // r0= A*x

  gpu_axpby(dr0,dtempv,1.0,-1.0,nrows,blocks,threads); // r0=1.0*rhs+-1.0r0 //y=ax+by

  gpu_yequalsx(dr0h,dr0,nrows,blocks,threads);  //r0h=r0

  gpu_yequalsconst(dn0,0.0,nrows,blocks,threads);  //n0=0.0 //memset???
  gpu_yequalsconst(dp0,0.0,nrows,blocks,threads);  //p0=0.0

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  //printf("temp1 %-le", temp1);
  //printf("rho1 %f", rho1);


  //for(int it=0;it<maxIt;it++){
  int it=0;
  do {

    rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);//rho1 =<r0,r0h>
    //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,blocks, threads);//rho1 =<r0,r0h>
    beta=(rho1/rho0)*(alpha/omega0);

    //    cout<<"rho1 "<<rho1<<" beta "<<beta<<endl;

    gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c

    gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag

    gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y

    temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);
    //temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows, blocks, threads);

    alpha=rho1/temp1;

    //       cout<<"temp1 "<<temp1<<" alpha "<<alpha<<endl;

    gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);

    gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s

    gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);

    gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);

    temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    //temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,blocks, threads);

    temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    //temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,blocks, threads);

    omega0= temp1/temp2;

    gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x

    gpu_axpy(dx,dz,omega0,nrows,blocks,threads);

    gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);

    temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    //temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,blocks, threads);
    temp1=sqrt(temp1);

  //cout<<it<<": "<<temp1<<endl;

    rho0=rho1;

    //if(temp1<tolmax){
    //  break;}}

    it++;
  }while(it<maxIt && temp1>tolmax);

}

void free_itsolver(itsolver *bicg)
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

  cudaFree(dr0);
  cudaFree(dr0h);
  cudaFree(dn0);
  cudaFree(dp0);
  cudaFree(dt);
  cudaFree(ds);
  cudaFree(dAx2);
  cudaFree(dy);
  cudaFree(dz);
  cudaFree(ddiag);
  cudaFree(daux);
  free(bicg->aux);

}

 /*
void setUpSolver(itsolver *bicg, double reltol, double *ewt, int tnrows,int tnnz,double *tA, int *tjA, int *tiA, int tmattype, int qmax, double *dACamp, double *dftempCamp);
{

  bicg.tolmax=reltol;

}
*/