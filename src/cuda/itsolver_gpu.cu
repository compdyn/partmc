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

//todo profiling del dot y ver cuanta occupancy me esta dando de shared memory porque me limita
//el numero de bloques que se ejecutan a la vez(solo se ejecutan a la vez en toda la function
// los bloques que "quepan" con la shared memory available: solution use cudastreams and launch instead
//of only 1 kernel use 2 or 4 to cubrir huecos (de memoria y eso), y tmb reducir la shared
//con una implementacion hibrida del dotxy

//*se puede probar la funcion en la cpu a ver que hace
__global__
void cudaSolveGPU(
        double *dA, int *djA, int *diA, double *dx, double *dtempv, //Input data
        int nrows, int blocks, int active_threads, int maxIt, int mattype,
        int n_cells, double tolmax, double *ddiag, //Init variables
        double *dr0, double *dr0h, double *dn0, double *dp0,
        double *dt, double *ds, double *dAx2, double *dy, double *dz,
        double *daux, // Auxiliary vectors
        double *aux_params
        //double *alpha, double *rho0, double* omega0, double *beta,
        //double *rho1, double *temp1, double *temp2 //Auxiliary parameters
)
{
  int tid = threadIdx.x;
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  double alpha,rho0,omega0,beta,rho1,temp1,temp2;
  //adjust pointers by blockid
  /*int idle_threads = blockDim.x-active_threads;
  dr0 -= blockIdx.x*idle_threads;
  dr0h -= blockIdx.x*idle_threads;
  dn0 -= blockIdx.x*idle_threads;
  dp0 -= blockIdx.x*idle_threads;
  dt -= blockIdx.x*idle_threads;
  ds -= blockIdx.x*idle_threads;
  dAx2 -= blockIdx.x*idle_threads;
  dy -= blockIdx.x*idle_threads;
  dz -= blockIdx.x*idle_threads;
  daux -= blockIdx.x*idle_threads;

  dA -= blockIdx.x*idle_threads;
  djA -= blockIdx.x*idle_threads;
  diA -= blockIdx.x*idle_threads;
  dx -= blockIdx.x*idle_threads;
  dtempv -= blockIdx.x*idle_threads;
*/

  alpha  = 1.0;
  rho0   = 1.0;
  omega0 = 1.0;

  /*beta = 1.0;
  rho1 = 1.0;
  temp1 = 1.0;
  temp2 = 1.0;*/

  int it=0;
  do
  {

    //not_working_threads = blockDim.x

    /*alpha = aux_params[0];
    rho0 = aux_params[1];
    omega0 = aux_params[2];
    beta = aux_params[3];
    rho1 = aux_params[4];
    temp1 = aux_params[5];
    temp2 = aux_params[6];*/

    //FOR SOME REASON this works when active threads is 1024...
    //if(tid < active_threads) //If tid is inside active_threads
    if(1)
    {
      //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);
      //todo if thread==1023 then dont do nothing and at the init modify pointers to
      //-1*blockidx to adjust delay caused by 1023 do nothing
      __syncthreads();
      cudaDevicedotxy(dr0, dr0h, &rho1, nrows);//&rho1
      __syncthreads();
      //cudaDevicereducey(daux, blocks);
      //__syncthreads();
      //Note more computations have been done on other blocks, but we only need block 0
      //rho1 = daux[0];//blockIdx.x//All threads access first value of daux, so it should be ok (plot twist: is not okay)
      beta = (rho1 / rho0) * (alpha / omega0);
/**/
      //gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c
      cudaDevicezaxpbypc(dp0, dr0, dn0, beta, -1.0 * omega0 * beta, nrows);   //z = ax + by + c
      //if (tid == 0)

      //gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag
      cudaDevicemultxy(dy, ddiag, dp0, nrows);

      //gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y
      cudaDevicesetconst(dn0, 0.0, nrows);
      cudaDeviceSpmvCSC(dn0, dy, nrows, dA, djA, diA);

      //temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);
      cudaDevicedotxy(dr0h, dn0, &temp1, nrows);
      __syncthreads();
      //cudaDevicereducey(daux, blocks);
      //temp1 = daux[0];
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
      cudaDevicedotxy(dz, dAx2, &temp1, nrows);
      __syncthreads();
      //cudaDevicereducey(daux, blocks);
      //temp1 = daux[0];

      //temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
      cudaDevicedotxy(dAx2, dAx2, &temp2, nrows);
      __syncthreads();
      //cudaDevicereducey(daux, blocks);
      //temp2 = daux[0];
      omega0 = temp1 / temp2;

      //gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x
      cudaDeviceaxpy(dx, dy, alpha, nrows); // x=alpha*y +x

      //gpu_axpy(dx,dz,omega0,nrows,blocks,threads);
      cudaDeviceaxpy(dx, dz, omega0, nrows);

      //gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);
      cudaDevicezaxpby(1.0, ds, -1.0 * omega0, dt, dr0, nrows);

      //temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
      cudaDevicedotxy(dr0, dr0, &temp1, nrows);
      __syncthreads();
      //cudaDevicereducey(daux, blocks);
      //temp1 = daux[0];
      temp1 = sqrt(temp1);

      rho0 = rho1;

      __syncthreads();
      /**/
    }
    else {
      //giving problems AS ALWAYS

      cudaDevicedotxy(dr0, dr0h, &rho1, 0);//Assing sdata[tid] to zero

      /*
      cudaDevicedotxy(dr0h, dn0, &rho1, 0);
      cudaDevicedotxy(dz, dAx2, &rho1, 0);
      cudaDevicedotxy(dAx2, dAx2, &rho1, 0);
      cudaDevicedotxy(dr0, dr0, &rho1, 0);
      */
    }

    it++;
  } while(it<maxIt && temp1>tolmax);//while(it<maxIt && temp1>tolmax);

  //Notice needed sync with other blocks if no multicells technique is applied
  if (id == 0)//why block 1 is different that block 0 in daux[0] if daux[0] should be global
  {
    aux_params[0]=alpha;
    aux_params[1]=rho0;
    aux_params[2]=omega0;
    aux_params[3]=beta;//0.01;
    aux_params[4]=rho1;//rho1
    aux_params[5]=temp1;
    aux_params[6]=temp2;
  }

}



//todo quit syncs with cpu/gpu
//Biconjugate gradient
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

  double *x=(double*)malloc(sizeof(double)*nrows);

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

  //Todo apply multicells (and remember threads should be power of two for dotxy reduce)
  //todo if difficulty in convergence (slow reaction vs fast), use 1 cell each block
  //todo if size_cell > n_threads_max warning to execute in CPU
  //Each block will compute only a cell/group of cells

  //printf("multi threads %d", multi_threads);

  /*//don't work because send kernels of 1 block size is inefficient
  int size_cell = nrows/n_cells;
  int mod_threads_sizecell = threads%size_cell;
  int extra_rows = (threads%size_cell)*blocks; //Extra rows need per block
  int total_rows = nrows + extra_rows;
  int total_blocks = (total_rows+threads-1)/threads; //threads = max_threads_block
  int n_rows_block = threads-mod_threads_sizecell; //Rows computed each block
*/

  //nrows = n_rows_block;
  //blocks = total_blocks;

  //three conditions:
  //1-rows_block is less than threads_max
  //2-rows_block is divisible by size_cell
  //3-Is a power of two
  //or better say, it must be power2 BUT only compute n_rows_cell... how the hell I did this?
  //initialising multiple kernels with only 1 block? (Seems easy way but not totally like the idea...)

  int n_aux_params=7;
  double *aux_params;
  aux_params=(double*)malloc(n_aux_params*sizeof(double));
  double *daux_params;
  cudaMalloc(&daux_params,n_aux_params*sizeof(double));
  //cudaMemcpy(bicg->djA,bicg->jA,7*sizeof(double),cudaMemcpyHostToDevice);

  int max_threads = bicg->threads;
  int size_cell = nrows/n_cells;
  int idle_threads = max_threads%size_cell;
  int active_threads = max_threads - idle_threads; //last multiple of size_cell before max_threads
  //Recalculate n_blocks
  //int multi_blocks = (nrows+multi_threads-1)/multi_threads; //threads = max_threads_block

  threads = active_threads;//active_threads;//bicg->threads;
  blocks = (nrows+active_threads-1)/active_threads; //blocks counting some threads are not working
  int nrows2 = nrows;//nrows+idle_threads*(blocks-1);

  //int redsize= sqrt(blocks) +1;
  //redsize=pow(2,redsize);

  //cudaMemcpy(x,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
  //printf("dx0 1 %f %f\n", x[0], x[1]);

  /**/aux_params[0] = alpha;
  aux_params[1] = rho0;
  aux_params[2] = omega0;
  aux_params[3] = beta;
  aux_params[4] = rho1;
  aux_params[5] = temp1;
  aux_params[6] = temp2;

  cudaMemcpy(daux_params, aux_params, n_aux_params * sizeof(double), cudaMemcpyHostToDevice);
  cudaSolveGPU << < blocks, active_threads, threads * sizeof(double) >> >//threads
                                            (dA, djA, diA, dx, dtempv, nrows2, blocks, threads, maxIt, mattype, n_cells,
                                                    tolmax, ddiag, dr0, dr0h, dn0, dp0, dt, ds, dAx2, dy, dz, daux,
                                                    daux_params
                                            );
  cudaDeviceSynchronize();
  cudaMemcpy(aux_params, daux_params, n_aux_params * sizeof(double), cudaMemcpyDeviceToHost);

  alpha = aux_params[0];
  rho0 = aux_params[1];
  omega0 = aux_params[2];
  beta = aux_params[3];
  rho1 = aux_params[4];
  temp1 = aux_params[5];
  temp2 = aux_params[6];
  printf("temp1 %-le", temp1);

  //for(int it=0;it<maxIt;it++){
  //int it=0;
  //do {

    //todo maybe set dotxy and reducey to work with 1023 threads to check if same result
    //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,(blocks + 1) / 2, threads);//rho1 =<r0,r0h>
    //rho1=gpu_dotxy(dr0, dr0h, aux, daux, nrows,blocks, threads);//rho1 =<r0,r0h>
    //beta=(rho1/rho0)*(alpha/omega0);



    //printf("rho1 %f", rho1);
/**/
/*
    //Si cada uno lo hace con los threads de su bloque porque usa la mitad de bloques

    //    cout<<"rho1 "<<rho1<<" beta "<<beta<<endl;

    gpu_zaxpbypc(dp0,dr0,dn0,beta,-1.0*omega0*beta,nrows,blocks,threads);   //z = ax + by + c

    gpu_multxy(dy,ddiag,dp0,nrows,blocks,threads);  // precond y= p0*diag

    gpu_spmv(dn0,dy,nrows,dA,djA,diA,mattype,blocks,threads);  // n0= A*y

    //temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=gpu_dotxy(dr0h, dn0, aux, daux, nrows, blocks, threads);

    alpha=rho1/temp1;

    //       cout<<"temp1 "<<temp1<<" alpha "<<alpha<<endl;

    gpu_zaxpby(1.0,dr0,-1.0*alpha,dn0,ds,nrows,blocks,threads);

    gpu_multxy(dz,ddiag,ds,nrows,blocks,threads); // precond z=diag*s

    gpu_spmv(dt,dz,nrows,dA,djA,diA,mattype,blocks,threads);

    gpu_multxy(dAx2,ddiag,dt,nrows,blocks,threads);

    //temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=gpu_dotxy(dz, dAx2, aux, daux, nrows,blocks, threads);

    //temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp2=gpu_dotxy(dAx2, dAx2, aux, daux, nrows,blocks, threads);

    omega0= temp1/temp2;

    gpu_axpy(dx,dy,alpha,nrows,blocks,threads); // x=alpha*y +x

    gpu_axpy(dx,dz,omega0,nrows,blocks,threads);

    gpu_zaxpby(1.0,ds,-1.0*omega0,dt,dr0,nrows,blocks,threads);

    //temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,(blocks + 1) / 2, threads);
    temp1=gpu_dotxy(dr0, dr0, aux, daux, nrows,blocks, threads);
    temp1=sqrt(temp1);

    //cout<<it<<": "<<temp1<<endl;


    if(temp1<tolmax){
      break;
    }
    //rho0=rho1;
    // }
    it++;
  }while(it<maxIt && temp1>tolmax);
*/
  cudaFreeMem(daux_params);

  //cudaMemcpy(x,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
  //printf("dx0 1 %f %f\n", x[0], x[1]);

}


//todo free
//void free_gpu(itsolver *bicg)

 /*
void setUpSolver(itsolver *bicg, double reltol, double *ewt, int tnrows,int tnnz,double *tA, int *tjA, int *tiA, int tmattype, int qmax, double *dACamp, double *dftempCamp);
{

  bicg.tolmax=reltol;

}
*/

//todo recover when multi-cells in independent gpu blocks
/*
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
//}

/*  int n_aux_params=7;
  double *aux_params;
  aux_params=(double*)malloc(n_aux_params*sizeof(double));
  double *daux_params;
  cudaMalloc(&daux_params,n_aux_params*sizeof(double));
  //cudaMemcpy(bicg->djA,bicg->jA,7*sizeof(double),cudaMemcpyHostToDevice);

  int redsize= sqrt(blocks) +1;
  redsize=pow(2,redsize);
*/

/*aux_params[0]=alpha;
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

//printf("rho1 %f", rho1);
//cudaFreeMem(daux_params);
*/