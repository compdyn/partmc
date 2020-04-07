/* Copyright (C) 2020 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * ODE GPU solver
 *
 */

#include "itsolver_gpu.h"

extern "C" {
#include "cvode_gpu.h"
#include "../camp_solver.h"
}

//todo good error code checking:
//Each function at the end update his own cv_mem->function_status code and return the same status.
//IF at the end of the function none of the subfunctions returned error code, return success.
//IF not, return error code and print the stack of status code to see the failed functions.

//this means translate all this flags defined to a proper status memory array (and map the array
// with the names of each function, maybe define a structure or something).
#define CV_SUCCESS               0

#define DO_ERROR_TEST    +2
#define PREDICT_AGAIN    +3
#define CONV_FAIL        +4
#define TRY_AGAIN        +5
#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8
#define RHSFUNC_RECVR    +9

#define NUM_TESTS    5     /* number of error test quantities     */

/*=================================================================*/
/*             CVODE Private Constants                             */
/*=================================================================*/

//#define ZERO    RCONST(0.0)     /* real 0.0     */
//#define TINY    RCONST(1.0e-10) /* small number */
#define PT1     RCONST(0.1)     /* real 0.1     */
#define POINT2  RCONST(0.2)     /* real 0.2     */
#define FOURTH  RCONST(0.25)    /* real 0.25    */
//#define HALF    RCONST(0.5)     /* real 0.5     */
//#define ONE     RCONST(1.0)     /* real 1.0     */
#define TWO     RCONST(2.0)     /* real 2.0     */
#define THREE   RCONST(3.0)     /* real 3.0     */
#define FOUR    RCONST(4.0)     /* real 4.0     */
#define FIVE    RCONST(5.0)     /* real 5.0     */
#define TWELVE  RCONST(12.0)    /* real 12.0    */
#define HUNDRED RCONST(100.0)   /* real 100.0   */
#define PMC_TINY RCONST(1.0e-30) /* small number for PMC */

/*=================================================================*/
/*             CVODE Routine-Specific Constants                    */
/*=================================================================*/

/*
 * Control constants for lower-level functions used by cvStep
 * ----------------------------------------------------------
 *
 * cvHin return values:
 *    CV_SUCCESS
 *    CV_RHSFUNC_FAIL
 *    CV_TOO_CLOSE
 *
 * cvStep control constants:
 *    DO_ERROR_TEST
 *    PREDICT_AGAIN
 *
 * cvStep return values:
 *    CV_SUCCESS,
 *    CV_LSETUP_FAIL,  CV_LSOLVE_FAIL,
 *    CV_RHSFUNC_FAIL, CV_RTFUNC_FAIL
 *    CV_CONV_FAILURE, CV_ERR_FAILURE,
 *    CV_FIRST_RHSFUNC_ERR
 *
 * cvNls input nflag values:
 *    FIRST_CALL
 *    PREV_CONV_FAIL
 *    PREV_ERR_FAIL
 *
 * cvNls return values:
 *    CV_SUCCESS,
 *    CV_LSETUP_FAIL, CV_LSOLVE_FAIL, CV_RHSFUNC_FAIL,
 *    CONV_FAIL, RHSFUNC_RECVR
 *
 * cvNewtonIteration return values:
 *    CV_SUCCESS,
 *    CV_LSOLVE_FAIL, CV_RHSFUNC_FAIL
 *    CONV_FAIL, RHSFUNC_RECVR,
 *    TRY_AGAIN
 *
 */

#define DO_ERROR_TEST    +2
#define PREDICT_AGAIN    +3

#define CONV_FAIL        +4
#define TRY_AGAIN        +5

#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8

#define RHSFUNC_RECVR    +9

/*
 * Control constants for lower-level rootfinding functions
 * -------------------------------------------------------
 *
 * cvRcheck1 return values:
 *    CV_SUCCESS,
 *    CV_RTFUNC_FAIL,
 * cvRcheck2 return values:
 *    CV_SUCCESS
 *    CV_RTFUNC_FAIL,
 *    CLOSERT
 *    RTFOUND
 * cvRcheck3 return values:
 *    CV_SUCCESS
 *    CV_RTFUNC_FAIL,
 *    RTFOUND
 * cvRootfind return values:
 *    CV_SUCCESS
 *    CV_RTFUNC_FAIL,
 *    RTFOUND
 */

#define RTFOUND          +1
#define CLOSERT          +3

/*
 * Control constants for tolerances
 * --------------------------------
 */

#define CV_NN  0
#define CV_SS  1
#define CV_SV  2
#define CV_WF  3

/*
 * Algorithmic constants
 * ---------------------
 *
 * CVodeGetDky and cvStep
 *
 *    FUZZ_FACTOR
 *
 * cvHin
 *
 *    HLB_FACTOR
 *    HUB_FACTOR
 *    H_BIAS
 *    MAX_ITERS
 *
 * CVodeCreate
 *
 *   CORTES
 *
 * cvStep
 *
 *    THRESH
 *    ETAMX1
 *    ETAMX2
 *    ETAMX3
 *    ETAMXF
 *    ETAMIN
 *    ETACF
 *    ADDON
 *    BIAS1
 *    BIAS2
 *    BIAS3
 *    ONEPSM
 *
 *    SMALL_NST   nst > SMALL_NST => use ETAMX3
 *    MXNCF       max no. of convergence failures during one step try
 *    MXNEF       max no. of error test failures during one step try
 *    MXNEF1      max no. of error test failures before forcing a reduction of order
 *    SMALL_NEF   if an error failure occurs and SMALL_NEF <= nef <= MXNEF1, then
 *                reset eta =  SUNMIN(eta, ETAMXF)
 *    LONG_WAIT   number of steps to wait before considering an order change when
 *                q==1 and MXNEF1 error test failures have occurred
 *
 * cvNls
 *
 *    NLS_MAXCOR  maximum no. of corrector iterations for the nonlinear solver
 *    CRDOWN      constant used in the estimation of the convergence rate (crate)
 *                of the iterates for the nonlinear equation
 *    DGMAX       iter == CV_NEWTON, |gamma/gammap-1| > DGMAX => call lsetup
 *    RDIV        declare divergence if ratio del/delp > RDIV
 *    MSBP        max no. of steps between lsetup calls
 *
 */

#define FUZZ_FACTOR RCONST(100.0)

#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define MAX_ITERS  4000

#define CORTES RCONST(0.1)

#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0)
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ETACF  RCONST(0.25)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)
#define ONEPSM RCONST(1.000001)

#define SMALL_NST    10
#define MXNCF        10
#define MXNEF         7
#define MXNEF1        3
#define SMALL_NEF     2
#define LONG_WAIT    10

#define NLS_MAXCOR 3
#define CRDOWN RCONST(0.3)
#define DGMAX  RCONST(0.3)

#define RDIV      TWO
#define MSBP       20

#ifdef PMC_DEBUG_GPU
  int counterSendInit=0;
  int counterMatScaleAddI=0;
  int counterMatScaleAddISendA=0;
  int counterMatCopy=0;
  int counterNewtonIt=0;
  int counterLinSolSetup=0;
  int counterLinSolSolve=0;
  int counterJacCVODE=0;
  int countercvStep=0;
  int counterprecvStep=0;
  int counterDerivNewton=0;
  int counterBiConjGrad=0;
  int counterDerivSolve=0;

  double timeNewtonSendInit=PMC_TINY;
  double timeMatScaleAddI=PMC_TINY;
  double timeMatScaleAddISendA=PMC_TINY;
  double timeMatCopy=PMC_TINY;
  double timeNewtonIt=PMC_TINY;
  double timeLinSolSetup=PMC_TINY;
  double timeLinSolSolve=PMC_TINY;
  double timeJacCVODE=PMC_TINY;
  double timecvStep=PMC_TINY;
  double timeprecvStep=PMC_TINY;
  double timeDerivNewton=PMC_TINY;
  double timeBiConjGrad=PMC_TINY;
  double timeDerivSolve=PMC_TINY;
#endif

void allocSolverGPU(CVodeMem cv_mem, SolverData *sd)
{
  itsolver *bicg = &(sd->bicg);
  ModelData *md = &(sd->model_data);
  CVDlsMem cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  SUNMatrix J = cvdls_mem->A;

  double *ewt = NV_DATA_S(cv_mem->cv_ewt);
  double *tempv = NV_DATA_S(cv_mem->cv_tempv);

  //Init GPU ODE solver variables
  //Linking Matrix data, later this data must be allocated in GPU
  bicg->nnz=SM_NNZ_S(J);
  bicg->nrows=SM_NP_S(J);
  bicg->mattype=1; //CSC
  bicg->A=(double*)SM_DATA_S(J);
  bicg->jA=(int*)SM_INDEXVALS_S(J);
  bicg->iA=(int*)SM_INDEXPTRS_S(J);

  int device=0;//Selected GPU
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  bicg->threads=prop.maxThreadsPerBlock;//1024; //128 set at max gpu
  bicg->blocks=(bicg->nrows+bicg->threads-1)/bicg->threads;
  bicg->n_cells=md->n_cells;
  bicg->dA=md->jac_gpu_data;//set itsolver gpu pointer to jac pointer initialized at camp
  bicg->dftemp=md->deriv_gpu_data;

  // Allocating matrix data to the GPU
  cudaMalloc((void**)&bicg->djA,bicg->nnz*sizeof(int));
  cudaMalloc((void**)&bicg->diA,(bicg->nrows+1)*sizeof(int));

  cudaMalloc((void**)&bicg->dewt,bicg->nrows*sizeof(double));
  cudaMalloc((void**)&bicg->dacor,bicg->nrows*sizeof(double));
  cudaMalloc((void**)&bicg->dtempv,bicg->nrows*sizeof(double));
  cudaMalloc((void**)&bicg->dzn,bicg->nrows*(cv_mem->cv_qmax+1)*sizeof(double));
  cudaMalloc((void**)&bicg->dcv_y,bicg->nrows*sizeof(double));
  cudaMalloc((void**)&bicg->dx,bicg->nrows*sizeof(double));

  cudaMemcpy(bicg->djA,bicg->jA,bicg->nnz*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->diA,bicg->iA,(bicg->nrows+1)*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dewt,ewt,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dacor,ewt,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dftemp,ewt,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dx,tempv,bicg->nnz*sizeof(double),cudaMemcpyHostToDevice);

  bicg->maxIt=100;
  bicg->tolmax=1e-10; //cv_mem->cv_reltol CAMP selected accuracy (1e-8) //1e-10;//1e-6

  //Init Linear Solver variables
  createSolver(bicg);

  //Alloc struct gpu //not working for some reason (maybe in the future?)
/*  //printf("Size of struct %d \n", sizeof(itsolver));
  cudaMalloc((void**)&dbicg,sizeof(itsolver));
  //Copy struct gpu
  cudaMemcpy(dbicg,bicg,sizeof(itsolver),cudaMemcpyHostToDevice);
  cudaMemcpy(&(dbicg->dx),&(bicg->dx),sizeof(dbicg->dx),cudaMemcpyHostToDevice);
*/
}

int CVode_gpu2(void *cvode_mem, realtype tout, N_Vector yout,
          realtype *tret, int itask, SolverData *sd)
{
  CVodeMem cv_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ir, ier, irfndp;
  int ewtsetOK;
  realtype troundoff, tout_hin, rh, nrm;
  booleantype inactive_roots;

#ifdef PMC_DEBUG_GPU
  clock_t start;
  start=clock();
#endif

  /*
   * -------------------------------------
   * 1. Check and process inputs
   * -------------------------------------
   */

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "CVode", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if cvode_mem was allocated */
  if (cv_mem->cv_MallocDone == SUNFALSE) {
    cvProcessError(cv_mem, CV_NO_MALLOC, "CVODE", "CVode", MSGCV_NO_MALLOC);
    return(CV_NO_MALLOC);
  }

  /* Check for yout != NULL */
  if ((cv_mem->cv_y = yout) == NULL) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode", MSGCV_YOUT_NULL);
    return(CV_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode", MSGCV_TRET_NULL);
    return(CV_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != CV_NORMAL) && (itask != CV_ONE_STEP) ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode", MSGCV_BAD_ITASK);
    return(CV_ILL_INPUT);
  }

  if (itask == CV_NORMAL) cv_mem->cv_toutc = tout;
  cv_mem->cv_taskc = itask;

  /*
   * ----------------------------------------
   * 2. Initializations performed only at
   *    the first step (nst=0):
   *    - initial setup
   *    - initialize Nordsieck history array
   *    - compute initial step size
   *    - check for approach to tstop
   *    - check for approach to a root
   * ----------------------------------------
   */

  if (cv_mem->cv_nst == 0) {

    cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;

    ier = cvInitialSetup_gpu2(cv_mem);
    if (ier!= CV_SUCCESS) return(ier);

    /* Call f at (t0,y0), set zn[1] = y'(t0),
       set initial h (from H0 or cvHin), and scale zn[1] by h.
       Also check for zeros of root function g at and near t0.    */
    //retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_zn[0],
    //                      cv_mem->cv_zn[1], cv_mem->cv_user_data);
    retval = f(cv_mem->cv_tn, cv_mem->cv_zn[0], cv_mem->cv_zn[1], cv_mem->cv_user_data);

    N_VScale(ONE, cv_mem->cv_zn[0], yout);

    cv_mem->cv_nfe++;
    if (retval < 0) {
      cvProcessError(cv_mem, CV_RHSFUNC_FAIL, "CVODE", "CVode",
                     MSGCV_RHSFUNC_FAILED, cv_mem->cv_tn);
      return(CV_RHSFUNC_FAIL);
    }
    if (retval > 0) {
      cvProcessError(cv_mem, CV_FIRST_RHSFUNC_ERR, "CVODE", "CVode",
                     MSGCV_RHSFUNC_FIRST);
      return(CV_FIRST_RHSFUNC_ERR);
    }

    /* Test input tstop for legality. */

    if (cv_mem->cv_tstopset) {
      if ( (cv_mem->cv_tstop - cv_mem->cv_tn)*(tout - cv_mem->cv_tn) <= ZERO ) {
        cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode",
                       MSGCV_BAD_TSTOP, cv_mem->cv_tstop, cv_mem->cv_tn);
        return(CV_ILL_INPUT);
      }
    }

    /* Set initial h (from H0 or cvHin). */

    cv_mem->cv_h = cv_mem->cv_hin;
    if ( (cv_mem->cv_h != ZERO) && ((tout-cv_mem->cv_tn)*cv_mem->cv_h < ZERO) ) {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode", MSGCV_BAD_H0);
      return(CV_ILL_INPUT);
    }
    if (cv_mem->cv_h == ZERO) {
      tout_hin = tout;
      if ( cv_mem->cv_tstopset && (tout-cv_mem->cv_tn)*(tout-cv_mem->cv_tstop) > ZERO )
        tout_hin = cv_mem->cv_tstop;
      hflag = cvHin_gpu2(cv_mem, tout_hin); //set cv_y
      if (hflag != CV_SUCCESS) {
        istate = cvHandleFailure_gpu2(cv_mem, hflag);
        return(istate);
      }
    }
    rh = SUNRabs(cv_mem->cv_h)*cv_mem->cv_hmax_inv;
    if (rh > ONE) cv_mem->cv_h /= rh;
    if (SUNRabs(cv_mem->cv_h) < cv_mem->cv_hmin)
      cv_mem->cv_h *= cv_mem->cv_hmin/SUNRabs(cv_mem->cv_h);

    /* Check for approach to tstop */

    if (cv_mem->cv_tstopset) {
      if ( (cv_mem->cv_tn + cv_mem->cv_h - cv_mem->cv_tstop)*cv_mem->cv_h > ZERO )
        cv_mem->cv_h = (cv_mem->cv_tstop - cv_mem->cv_tn)*(ONE-FOUR*cv_mem->cv_uround);
    }

    /* Scale zn[1] by h.*/

    cv_mem->cv_hscale = cv_mem->cv_h;
    cv_mem->cv_h0u    = cv_mem->cv_h;
    cv_mem->cv_hprime = cv_mem->cv_h;

    N_VScale(cv_mem->cv_h, cv_mem->cv_zn[1], cv_mem->cv_zn[1]);

    /* Try to improve initial guess of zn[1] */
    if (cv_mem->cv_ghfun) {
      N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, cv_mem->cv_zn[1], cv_mem->cv_tempv1);
      cv_mem->cv_ghfun(cv_mem->cv_tn + cv_mem->cv_h, cv_mem->cv_h, cv_mem->cv_tempv1,
                       cv_mem->cv_zn[0], cv_mem->cv_zn[1], cv_mem->cv_user_data,
                       cv_mem->cv_tempv2, cv_mem->cv_acor_init);
    }

    /* Check for zeros of root function g at and near t0. */

    if (cv_mem->cv_nrtfn > 0) {

      retval = cvRcheck1_gpu2(cv_mem);

      if (retval == CV_RTFUNC_FAIL) {
        cvProcessError(cv_mem, CV_RTFUNC_FAIL, "CVODE", "cvRcheck1",
                       MSGCV_RTFUNC_FAILED, cv_mem->cv_tn);
        return(CV_RTFUNC_FAIL);
      }

    }

  } /* end of first call block */

  /*
   * ------------------------------------------------------
   * 3. At following steps, perform stop tests:
   *    - check for root in last step
   *    - check if we passed tstop
   *    - check if we passed tout (NORMAL mode)
   *    - check if current tn was returned (ONE_STEP mode)
   *    - check if we are close to tstop
   *      (adjust step size if needed)
   * -------------------------------------------------------
   */

  if (cv_mem->cv_nst > 0) {

    /* Estimate an infinitesimal time interval to be used as
       a roundoff for time quantities (based on current time
       and step size) */
    troundoff = FUZZ_FACTOR*cv_mem->cv_uround*(SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_h));

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If itask = CV_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */
    if (cv_mem->cv_nrtfn > 0) {

      irfndp = cv_mem->cv_irfnd;

      retval = cvRcheck2_gpu2(cv_mem);

      if (retval == CLOSERT) {
        cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "cvRcheck2",
                       MSGCV_CLOSE_ROOTS, cv_mem->cv_tlo);
        return(CV_ILL_INPUT);
      } else if (retval == CV_RTFUNC_FAIL) {
        cvProcessError(cv_mem, CV_RTFUNC_FAIL, "CVODE", "cvRcheck2",
                       MSGCV_RTFUNC_FAILED, cv_mem->cv_tlo);
        return(CV_RTFUNC_FAIL);
      } else if (retval == RTFOUND) {
        cv_mem->cv_tretlast = *tret = cv_mem->cv_tlo;
        return(CV_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( SUNRabs(cv_mem->cv_tn - cv_mem->cv_tretlast) > troundoff ) {

        retval = cvRcheck3_gpu2(cv_mem);

        if (retval == CV_SUCCESS) {     /* no root found */
          cv_mem->cv_irfnd = 0;
          if ((irfndp == 1) && (itask == CV_ONE_STEP)) {
            cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
            N_VScale(ONE, cv_mem->cv_zn[0], yout);
            return(CV_SUCCESS);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          cv_mem->cv_irfnd = 1;
          cv_mem->cv_tretlast = *tret = cv_mem->cv_tlo;
          return(CV_ROOT_RETURN);
        } else if (retval == CV_RTFUNC_FAIL) {  /* g failed */
          cvProcessError(cv_mem, CV_RTFUNC_FAIL, "CVODE", "cvRcheck3",
                         MSGCV_RTFUNC_FAILED, cv_mem->cv_tlo);
          return(CV_RTFUNC_FAIL);
        }

      }

    } /* end of root stop check */

    /* In CV_NORMAL mode, test if tout was reached */
    if ( (itask == CV_NORMAL) && ((cv_mem->cv_tn-tout)*cv_mem->cv_h >= ZERO) ) {
      cv_mem->cv_tretlast = *tret = tout;
      ier =  CVodeGetDky(cv_mem, tout, 0, yout);
      if (ier != CV_SUCCESS) {
        cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode",
                       MSGCV_BAD_TOUT, tout);
        return(CV_ILL_INPUT);
      }
      return(CV_SUCCESS);
    }

    /* In CV_ONE_STEP mode, test if tn was returned */
    if ( itask == CV_ONE_STEP &&
         SUNRabs(cv_mem->cv_tn - cv_mem->cv_tretlast) > troundoff ) {
      cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], yout);
      return(CV_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( cv_mem->cv_tstopset ) {

      if ( SUNRabs(cv_mem->cv_tn - cv_mem->cv_tstop) <= troundoff) {
        ier =  CVodeGetDky(cv_mem, cv_mem->cv_tstop, 0, yout);
        if (ier != CV_SUCCESS) {
          cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode",
                         MSGCV_BAD_TSTOP, cv_mem->cv_tstop, cv_mem->cv_tn);
          return(CV_ILL_INPUT);
        }
        cv_mem->cv_tretlast = *tret = cv_mem->cv_tstop;
        cv_mem->cv_tstopset = SUNFALSE;
        return(CV_TSTOP_RETURN);
      }

      /* If next step would overtake tstop, adjust stepsize */
      if ( (cv_mem->cv_tn + cv_mem->cv_hprime - cv_mem->cv_tstop)*cv_mem->cv_h > ZERO ) {
        cv_mem->cv_hprime = (cv_mem->cv_tstop - cv_mem->cv_tn)*(ONE-FOUR*cv_mem->cv_uround);
        cv_mem->cv_eta = cv_mem->cv_hprime/cv_mem->cv_h;
      }

    }

  } /* end stopping tests block */

  /*
   * --------------------------------------------------
   * 4. Looping point for internal steps
   *
   *    4.1. check for errors (too many steps, too much
   *         accuracy requested, step size too small)
   *    4.2. take a new step (call cvStep)
   *    4.3. stop on error
   *    4.4. perform stop tests:
   *         - check for root in last step
   *         - check if tout was passed
   *         - check if close to tstop
   *         - check if in ONE_STEP mode (must return)
   * --------------------------------------------------
   */

  // GPU initializations
  //set_data_gpu2(cv_mem, sd);

  nstloc = 0;
  for(;;) {

    cv_mem->cv_next_h = cv_mem->cv_h;
    cv_mem->cv_next_q = cv_mem->cv_q;

    /* Reset and check ewt */
    if (cv_mem->cv_nst > 0) {

      ewtsetOK = cv_mem->cv_efun(cv_mem->cv_zn[0], cv_mem->cv_ewt, cv_mem->cv_e_data);
      //set here copy of ewt to gpu

      if (ewtsetOK != 0) {

        if (cv_mem->cv_itol == CV_WF)
          cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode",
                         MSGCV_EWT_NOW_FAIL, cv_mem->cv_tn);
        else
          cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVode",
                         MSGCV_EWT_NOW_BAD, cv_mem->cv_tn);

        istate = CV_ILL_INPUT;
        cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
        N_VScale(ONE, cv_mem->cv_zn[0], yout);
        break;

      }
    }

    /* Check for too many steps */
    if ( (cv_mem->cv_mxstep>0) && (nstloc >= cv_mem->cv_mxstep) ) {
      cvProcessError(cv_mem, CV_TOO_MUCH_WORK, "CVODE", "CVode",
                     MSGCV_MAX_STEPS, cv_mem->cv_tn);
      istate = CV_TOO_MUCH_WORK;
      cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(cv_mem->cv_zn[0], cv_mem->cv_ewt);
    cv_mem->cv_tolsf = cv_mem->cv_uround * nrm;
    if (cv_mem->cv_tolsf > ONE) {
      cvProcessError(cv_mem, CV_TOO_MUCH_ACC, "CVODE", "CVode",
                     MSGCV_TOO_MUCH_ACC, cv_mem->cv_tn);
      istate = CV_TOO_MUCH_ACC;
      cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], yout);
      cv_mem->cv_tolsf *= TWO;
      break;
    } else {
      cv_mem->cv_tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (cv_mem->cv_tn + cv_mem->cv_h == cv_mem->cv_tn) {
      cv_mem->cv_nhnil++;
      if (cv_mem->cv_nhnil <= cv_mem->cv_mxhnil)
        cvProcessError(cv_mem, CV_WARNING, "CVODE", "CVode",
                       MSGCV_HNIL, cv_mem->cv_tn, cv_mem->cv_h);
      if (cv_mem->cv_nhnil == cv_mem->cv_mxhnil)
        cvProcessError(cv_mem, CV_WARNING, "CVODE", "CVode", MSGCV_HNIL_DONE);
    }

#ifdef PMC_DEBUG_GPU
    timeprecvStep+= clock() - start;
    counterprecvStep++;

    start=clock();
#endif
    /* Call cvStep to take a step */
    //kflag = cvStep(cv_mem);
    kflag = cvStep_gpu2(sd, cv_mem);

#ifdef PMC_DEBUG_GPU
    timecvStep+= clock() - start;
    countercvStep++;
#endif

    /* Process failed step cases, and exit loop */
    if (kflag != CV_SUCCESS) {
      istate = cvHandleFailure_gpu2(cv_mem, kflag);
      cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], yout);
      break;
    }

    nstloc++;

    /* Check for root in last step taken. */
    if (cv_mem->cv_nrtfn > 0) {

      retval = cvRcheck3_gpu2(cv_mem);

      if (retval == RTFOUND) {  /* A new root was found */
        cv_mem->cv_irfnd = 1;
        istate = CV_ROOT_RETURN;
        cv_mem->cv_tretlast = *tret = cv_mem->cv_tlo;
        break;
      } else if (retval == CV_RTFUNC_FAIL) { /* g failed */
        cvProcessError(cv_mem, CV_RTFUNC_FAIL, "CVODE", "cvRcheck3",
                       MSGCV_RTFUNC_FAILED, cv_mem->cv_tlo);
        istate = CV_RTFUNC_FAIL;
        break;
      }

      /* If we are at the end of the first step and we still have
       * some event functions that are inactive, issue a warning
       * as this may indicate a user error in the implementation
       * of the root function. */

      if (cv_mem->cv_nst==1) {
        inactive_roots = SUNFALSE;
        for (ir=0; ir<cv_mem->cv_nrtfn; ir++) {
          if (!cv_mem->cv_gactive[ir]) {
            inactive_roots = SUNTRUE;
            break;
          }
        }
        if ((cv_mem->cv_mxgnull > 0) && inactive_roots) {
          cvProcessError(cv_mem, CV_WARNING, "CVODES", "CVode",
                         MSGCV_INACTIVE_ROOTS);
        }
      }

    }

    /* In NORMAL mode, check if tout reached */
    if ( (itask == CV_NORMAL) &&  (cv_mem->cv_tn-tout)*cv_mem->cv_h >= ZERO ) {
      istate = CV_SUCCESS;
      cv_mem->cv_tretlast = *tret = tout;
      (void) CVodeGetDky(cv_mem, tout, 0, yout);
      cv_mem->cv_next_q = cv_mem->cv_qprime;
      cv_mem->cv_next_h = cv_mem->cv_hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( cv_mem->cv_tstopset ) {

      troundoff = FUZZ_FACTOR*cv_mem->cv_uround*(SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_h));
      if ( SUNRabs(cv_mem->cv_tn - cv_mem->cv_tstop) <= troundoff) {
        (void) CVodeGetDky(cv_mem, cv_mem->cv_tstop, 0, yout);
        cv_mem->cv_tretlast = *tret = cv_mem->cv_tstop;
        cv_mem->cv_tstopset = SUNFALSE;
        istate = CV_TSTOP_RETURN;
        break;
      }

      if ( (cv_mem->cv_tn + cv_mem->cv_hprime - cv_mem->cv_tstop)*cv_mem->cv_h > ZERO ) {
        cv_mem->cv_hprime = (cv_mem->cv_tstop - cv_mem->cv_tn)*(ONE-FOUR*cv_mem->cv_uround);
        cv_mem->cv_eta = cv_mem->cv_hprime/cv_mem->cv_h;
      }

    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == CV_ONE_STEP) {
      istate = CV_SUCCESS;
      cv_mem->cv_tretlast = *tret = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], yout);
      cv_mem->cv_next_q = cv_mem->cv_qprime;
      cv_mem->cv_next_h = cv_mem->cv_hprime;
      break;
    }

  } /* end looping for internal steps */

  free_gpu2(cv_mem);

  return(istate);
}


/*-----------------------------------------------------------------*/

/*
 * CVodeGetDky
 *
 * This routine computes the k-th derivative of the interpolating
 * polynomial at the time t and stores the result in the vector dky.
 * The formula is:
 *         q
 *  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] ,
 *        j=k
 * where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 * zn[j] is the j-th column of the Nordsieck history array.
 *
 * This function is called by CVode with k = 0 and t = tout, but
 * may also be called directly by the user.
 */

int CVodeGetDky_gpu2(void *cvode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  CVodeMem cv_mem;

  /* Check all inputs for legality */

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeGetDky", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (dky == NULL) {
    cvProcessError(cv_mem, CV_BAD_DKY, "CVODE", "CVodeGetDky", MSGCV_NULL_DKY);
    return(CV_BAD_DKY);
  }

  if ((k < 0) || (k > cv_mem->cv_q)) {
    cvProcessError(cv_mem, CV_BAD_K, "CVODE", "CVodeGetDky", MSGCV_BAD_K);
    return(CV_BAD_K);
  }

  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * cv_mem->cv_uround * (SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_hu));
  if (cv_mem->cv_hu < ZERO) tfuzz = -tfuzz;
  tp = cv_mem->cv_tn - cv_mem->cv_hu - tfuzz;
  tn1 = cv_mem->cv_tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    cvProcessError(cv_mem, CV_BAD_T, "CVODE", "CVodeGetDky", MSGCV_BAD_T,
                   t, cv_mem->cv_tn-cv_mem->cv_hu, cv_mem->cv_tn);
    return(CV_BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - cv_mem->cv_tn) / cv_mem->cv_h;
  for (j=cv_mem->cv_q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == cv_mem->cv_q) {
      N_VScale(c, cv_mem->cv_zn[cv_mem->cv_q], dky);
    } else {
      N_VLinearSum(c, cv_mem->cv_zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(CV_SUCCESS);
  r = SUNRpowerI(cv_mem->cv_h,-k);
  N_VScale(r, dky, dky);
  return(CV_SUCCESS);
}

/*
 * CVodeFree
 *
 * This routine frees the problem memory allocated by CVodeInit.
 * Such memory includes all the vectors allocated by cvAllocVectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree).
 */

void CVodeFree_gpu2(void **cvode_mem)
{
  CVodeMem cv_mem;

  if (*cvode_mem == NULL) return;

  cv_mem = (CVodeMem) (*cvode_mem);

  cvFreeVectors_gpu2(cv_mem);

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  if (cv_mem->cv_nrtfn > 0) {
    free(cv_mem->cv_glo); cv_mem->cv_glo = NULL;
    free(cv_mem->cv_ghi); cv_mem->cv_ghi = NULL;
    free(cv_mem->cv_grout); cv_mem->cv_grout = NULL;
    free(cv_mem->cv_iroots); cv_mem->cv_iroots = NULL;
    free(cv_mem->cv_rootdir); cv_mem->cv_rootdir = NULL;
    free(cv_mem->cv_gactive); cv_mem->cv_gactive = NULL;
  }

  free(*cvode_mem);
  *cvode_mem = NULL;
}

/*
 * =================================================================
 *  Private Functions Implementation
 * =================================================================
 */

/*
 * cvCheckNvector
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns SUNFALSE.
 */

booleantype cvCheckNvector_gpu2(N_Vector tmpl)
{
  if((tmpl->ops->nvclone     == NULL) ||
     (tmpl->ops->nvdestroy   == NULL) ||
     (tmpl->ops->nvlinearsum == NULL) ||
     (tmpl->ops->nvconst     == NULL) ||
     (tmpl->ops->nvprod      == NULL) ||
     (tmpl->ops->nvdiv       == NULL) ||
     (tmpl->ops->nvscale     == NULL) ||
     (tmpl->ops->nvabs       == NULL) ||
     (tmpl->ops->nvinv       == NULL) ||
     (tmpl->ops->nvaddconst  == NULL) ||
     (tmpl->ops->nvmaxnorm   == NULL) ||
     (tmpl->ops->nvwrmsnorm  == NULL) ||
     (tmpl->ops->nvmin       == NULL))
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

/*
 * cvAllocVectors
 *
 * This routine allocates the CVODE vectors ewt, acor, tempv, ftemp, and
 * zn[0], ..., zn[maxord].
 * If all memory allocations are successful, cvAllocVectors returns SUNTRUE.
 * Otherwise all allocated memory is freed and cvAllocVectors returns SUNFALSE.
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 */

booleantype cvAllocVectors_gpu2(CVodeMem cv_mem, N_Vector tmpl)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */

  cv_mem->cv_ewt = N_VClone(tmpl);
  if (cv_mem->cv_ewt == NULL) return(SUNFALSE);

  cv_mem->cv_acor = N_VClone(tmpl);
  if (cv_mem->cv_acor == NULL) {
    N_VDestroy(cv_mem->cv_ewt);
    return(SUNFALSE);
  }

  cv_mem->cv_tempv = N_VClone(tmpl);
  if (cv_mem->cv_tempv == NULL) {
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }
  N_VConst(ZERO, cv_mem->cv_tempv);

  cv_mem->cv_tempv1 = N_VClone(tmpl);
  if (cv_mem->cv_tempv1 == NULL) {
    N_VDestroy(cv_mem->cv_tempv);
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }
  N_VConst(ZERO, cv_mem->cv_tempv1);

  cv_mem->cv_tempv2 = N_VClone(tmpl);
  if (cv_mem->cv_tempv2 == NULL) {
    N_VDestroy(cv_mem->cv_tempv);
    N_VDestroy(cv_mem->cv_tempv1);
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }
  N_VConst(ZERO, cv_mem->cv_tempv2);
  cv_mem->cv_acor_init = N_VClone(tmpl);
  if (cv_mem->cv_acor_init == NULL) {
    N_VDestroy(cv_mem->cv_tempv);
    N_VDestroy(cv_mem->cv_tempv1);
    N_VDestroy(cv_mem->cv_tempv2);
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }

  cv_mem->cv_last_yn = N_VClone(tmpl);
  if (cv_mem->cv_last_yn == NULL) {
    N_VDestroy(cv_mem->cv_acor_init);
    N_VDestroy(cv_mem->cv_tempv);
    N_VDestroy(cv_mem->cv_tempv1);
    N_VDestroy(cv_mem->cv_tempv2);
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }

  cv_mem->cv_ftemp = N_VClone(tmpl);
  if (cv_mem->cv_ftemp == NULL) {
    N_VDestroy(cv_mem->cv_last_yn);
    N_VDestroy(cv_mem->cv_acor_init);
    N_VDestroy(cv_mem->cv_tempv);
    N_VDestroy(cv_mem->cv_tempv1);
    N_VDestroy(cv_mem->cv_tempv2);
    N_VDestroy(cv_mem->cv_ewt);
    N_VDestroy(cv_mem->cv_acor);
    return(SUNFALSE);
  }

  /* Allocate zn[0] ... zn[qmax] */

  for (j=0; j <= cv_mem->cv_qmax; j++) {
    cv_mem->cv_zn[j] = N_VClone(tmpl);
    if (cv_mem->cv_zn[j] == NULL) {
      N_VDestroy(cv_mem->cv_ewt);
      N_VDestroy(cv_mem->cv_acor);
      N_VDestroy(cv_mem->cv_acor_init);
      N_VDestroy(cv_mem->cv_last_yn);
      N_VDestroy(cv_mem->cv_tempv);
      N_VDestroy(cv_mem->cv_tempv1);
      N_VDestroy(cv_mem->cv_tempv2);
      N_VDestroy(cv_mem->cv_ftemp);
      for (i=0; i < j; i++) N_VDestroy(cv_mem->cv_zn[i]);
      return(SUNFALSE);
    }
  }

  /* Update solver workspace lengths  */
  cv_mem->cv_lrw += (cv_mem->cv_qmax + 5)*cv_mem->cv_lrw1;
  cv_mem->cv_liw += (cv_mem->cv_qmax + 5)*cv_mem->cv_liw1;

  /* Store the value of qmax used here */
  cv_mem->cv_qmax_alloc = cv_mem->cv_qmax;

  return(SUNTRUE);
}

/*
 * cvFreeVectors
 *
 * This routine frees the CVODE vectors allocated in cvAllocVectors.
 */

void cvFreeVectors_gpu2(CVodeMem cv_mem)
{
  int j, maxord;

  maxord = cv_mem->cv_qmax_alloc;

  N_VDestroy(cv_mem->cv_ewt);
  N_VDestroy(cv_mem->cv_acor);
  N_VDestroy(cv_mem->cv_tempv);
  N_VDestroy(cv_mem->cv_tempv1);
  N_VDestroy(cv_mem->cv_tempv2);
  N_VDestroy(cv_mem->cv_acor_init);
  N_VDestroy(cv_mem->cv_last_yn);
  N_VDestroy(cv_mem->cv_ftemp);
  for (j=0; j <= maxord; j++) N_VDestroy(cv_mem->cv_zn[j]);

  cv_mem->cv_lrw -= (maxord + 5)*cv_mem->cv_lrw1;
  cv_mem->cv_liw -= (maxord + 5)*cv_mem->cv_liw1;

  if (cv_mem->cv_VabstolMallocDone) {
    N_VDestroy(cv_mem->cv_Vabstol);
    cv_mem->cv_lrw -= cv_mem->cv_lrw1;
    cv_mem->cv_liw -= cv_mem->cv_liw1;
  }
}

/*
 * cvInitialSetup
 *
 * This routine performs input consistency checks at the first step.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
 */

int cvInitialSetup_gpu2(CVodeMem cv_mem)
{
  int ier;

  /* Did the user specify tolerances? */
  if (cv_mem->cv_itol == CV_NN) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "cvInitialSetup", MSGCV_NO_TOLS);
    return(CV_ILL_INPUT);
  }

  /* Set data for efun */
  if (cv_mem->cv_user_efun) cv_mem->cv_e_data = cv_mem->cv_user_data;
  else                      cv_mem->cv_e_data = cv_mem;

  /* Load initial error weights */
  ier = cv_mem->cv_efun(cv_mem->cv_zn[0], cv_mem->cv_ewt, cv_mem->cv_e_data);
  if (ier != 0) {
    if (cv_mem->cv_itol == CV_WF)
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "cvInitialSetup", MSGCV_EWT_FAIL);
    else
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "cvInitialSetup", MSGCV_BAD_EWT);
    return(CV_ILL_INPUT);
  }

  /* Check if lsolve function exists (if needed) and call linit function (if it exists) */
  if (cv_mem->cv_iter == CV_NEWTON) {
    if (cv_mem->cv_lsolve == NULL) {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "cvInitialSetup", MSGCV_LSOLVE_NULL);
      return(CV_ILL_INPUT);
    }
    if (cv_mem->cv_linit != NULL) {
      ier = cv_mem->cv_linit(cv_mem);
      if (ier != 0) {
        cvProcessError(cv_mem, CV_LINIT_FAIL, "CVODE", "cvInitialSetup", MSGCV_LINIT_FAIL);
        return(CV_LINIT_FAIL);
      }
    }
  }

  return(CV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * PRIVATE FUNCTIONS FOR CVODE
 * -----------------------------------------------------------------
 */

/*
 * cvHin
 *
 * This routine computes a tentative initial step size h0.
 * If tout is too close to tn (= t0), then cvHin returns CV_TOO_CLOSE
 * and h remains uninitialized. Note that here tout is either the value
 * passed to CVode at the first call or the value of tstop (if tstop is
 * enabled and it is closer to t0=tn than tout).
 * If the RHS function fails unrecoverably, cvHin returns CV_RHSFUNC_FAIL.
 * If the RHS function fails recoverably too many times and recovery is
 * not possible, cvHin returns CV_REPTD_RHSFUNC_ERR.
 * Otherwise, cvHin sets h to the chosen value h0 and returns CV_SUCCESS.
 *
 * The algorithm used seeks to find h0 as a solution of
 *       (WRMS norm of (h0^2 ydd / 2)) = 1,
 * where ydd = estimated second derivative of y.
 *
 * We start with an initial estimate equal to the geometric mean of the
 * lower and upper bounds on the step size.
 *
 * Loop up to MAX_ITERS times to find h0.
 * Stop if new and previous values differ by a factor < 2.
 * Stop if hnew/hg > 2 after one iteration, as this probably means
 * that the ydd value is bad because of cancellation error.
 *
 * For each new proposed hg, we allow MAX_ITERS attempts to
 * resolve a possible recoverable failure from f() by reducing
 * the proposed stepsize by a factor of 0.2. If a legal stepsize
 * still cannot be found, fall back on a previous value if possible,
 * or else return CV_REPTD_RHSFUNC_ERR.
 *
 * Finally, we apply a bias (0.5) and verify that h0 is within bounds.
 */
int cvHin_gpu2(CVodeMem cv_mem, realtype tout)
{
  int retval, sign, count1, count2;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  booleantype hgOK, hnewOK;

  /* If tout is too close to tn, give up */

  if ((tdiff = tout-cv_mem->cv_tn) == ZERO) return(CV_TOO_CLOSE);

  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = SUNRabs(tdiff);
  tround = cv_mem->cv_uround * SUNMAX(SUNRabs(cv_mem->cv_tn), SUNRabs(tout));

  if (tdist < TWO*tround) return(CV_TOO_CLOSE);

  /*
     Set lower and upper bounds on h0, and take geometric mean
     as first trial value.
     Exit with this value if the bounds cross each other.
  */

  hlb = HLB_FACTOR * tround;
  hub = cvUpperBoundH0_gpu2(cv_mem, tdist);

  hg  = SUNRsqrt(hlb*hub);

  if (hub < hlb) {
    if (sign == -1) cv_mem->cv_h = -hg;
    else            cv_mem->cv_h =  hg;
    return(CV_SUCCESS);
  }

  /* Outer loop */

  hnewOK = SUNFALSE;
  hs = hg;         /* safeguard against 'uninitialized variable' warning */

  for(count1 = 1; count1 <= MAX_ITERS; count1++) {

    /* Attempts to estimate ydd */

    hgOK = SUNFALSE;

    for (count2 = 1; count2 <= MAX_ITERS; count2++) {
      hgs = hg*sign;
      retval = cvYddNorm_gpu2(cv_mem, hgs, &yddnrm);
      /* If f() failed unrecoverably, give up */
      if (retval < 0) return(CV_RHSFUNC_FAIL);
      /* If successful, we can use ydd */
      if (retval == CV_SUCCESS) {hgOK = SUNTRUE; break;}
      /* f() failed recoverably; cut step size and test it again */
      hg *= POINT2;
    }

    /* If f() failed recoverably MAX_ITERS times */

    if (!hgOK) {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) return(CV_REPTD_RHSFUNC_ERR);
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* If the stopping criteria was met, or if this is the last pass, stop */
    if ( (hnewOK) || (count1 == MAX_ITERS))  {hnew = hg; break;}

    /* Propose new step size */
    hnew = (yddnrm*hub*hub > TWO) ? SUNRsqrt(TWO/yddnrm) : SUNRsqrt(hg*hub);
    hrat = hnew/hg;

    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO)) {
      hnewOK = SUNTRUE;
    }

    /* After one pass, if ydd seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO)) {
      hnew = hg;
      hnewOK = SUNTRUE;
    }

    /* Send this value back through f() */
    hg = hnew;

  }

  /* Apply bounds, bias factor, and attach sign */

  h0 = H_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  cv_mem->cv_h = h0;

  return(CV_SUCCESS);
}

/*
 * cvUpperBoundH0
 *
 * This routine sets an upper bound on abs(h0) based on
 * tdist = tn - t0 and the values of y[i]/y'[i].
 */

realtype cvUpperBoundH0_gpu2(CVodeMem cv_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  /*
   * Bound based on |y0|/|y0'| -- allow at most an increase of
   * HUB_FACTOR in y0 (based on a forward Euler step). The weight
   * factor is used as a safeguard against zero components in y0.
   */

  temp1 = cv_mem->cv_tempv;
  temp2 = cv_mem->cv_acor;

  N_VAbs(cv_mem->cv_zn[0], temp2);
  cv_mem->cv_efun(cv_mem->cv_zn[0], temp1, cv_mem->cv_e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(cv_mem->cv_zn[1], temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /*
   * bound based on tdist -- allow at most a step of magnitude
   * HUB_FACTOR * tdist
   */

  hub = HUB_FACTOR*tdist;

  /* Use the smaller of the two */

  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}

/*
 * cvYddNorm
 *
 * This routine computes an estimate of the second derivative of y
 * using a difference quotient, and returns its WRMS norm.
 */

int cvYddNorm_gpu2(CVodeMem cv_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  N_VLinearSum(hg, cv_mem->cv_zn[1], ONE, cv_mem->cv_zn[0], cv_mem->cv_y);
  //retval = cv_mem->cv_f(cv_mem->cv_tn+hg, cv_mem->cv_y,
  //                      cv_mem->cv_tempv, cv_mem->cv_user_data);
  retval = f(cv_mem->cv_tn+hg, cv_mem->cv_y, cv_mem->cv_tempv, cv_mem->cv_user_data);
  cv_mem->cv_nfe++;
  if (retval < 0) return(CV_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VLinearSum(ONE, cv_mem->cv_tempv, -ONE, cv_mem->cv_zn[1], cv_mem->cv_tempv);
  N_VScale(ONE/hg, cv_mem->cv_tempv, cv_mem->cv_tempv);

  *yddnrm = N_VWrmsNorm(cv_mem->cv_tempv, cv_mem->cv_ewt);

  return(CV_SUCCESS);
}




/*
 * -----------------------------------------------------------------
 * Functions for rootfinding
 * -----------------------------------------------------------------
 */

/*
 * cvRcheck1
 *
 * This routine completes the initialization of rootfinding memory
 * information, and checks whether g has a zero both at and very near
 * the initial point of the IVP.
 *
 * This routine returns an int equal to:
 *  CV_RTFUNC_FAIL < 0 if the g function failed, or
 *  CV_SUCCESS     = 0 otherwise.
 */

int cvRcheck1_gpu2(CVodeMem cv_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  for (i = 0; i < cv_mem->cv_nrtfn; i++) cv_mem->cv_iroots[i] = 0;
  cv_mem->cv_tlo = cv_mem->cv_tn;
  cv_mem->cv_ttol = (SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_h)) *
                    cv_mem->cv_uround*HUNDRED;

  /* Evaluate g at initial t and check for zero values. */
  retval = cv_mem->cv_gfun(cv_mem->cv_tlo, cv_mem->cv_zn[0],
                           cv_mem->cv_glo, cv_mem->cv_user_data);
  cv_mem->cv_nge = 1;
  if (retval != 0) return(CV_RTFUNC_FAIL);

  zroot = SUNFALSE;
  for (i = 0; i < cv_mem->cv_nrtfn; i++) {
    if (SUNRabs(cv_mem->cv_glo[i]) == ZERO) {
      zroot = SUNTRUE;
      cv_mem->cv_gactive[i] = SUNFALSE;
    }
  }
  if (!zroot) return(CV_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = SUNMAX(cv_mem->cv_ttol/SUNRabs(cv_mem->cv_h), PT1);
  smallh = hratio*cv_mem->cv_h;
  tplus = cv_mem->cv_tlo + smallh;
  N_VLinearSum(ONE, cv_mem->cv_zn[0], hratio,
               cv_mem->cv_zn[1], cv_mem->cv_y);
  retval = cv_mem->cv_gfun(tplus, cv_mem->cv_y,
                           cv_mem->cv_ghi, cv_mem->cv_user_data);
  cv_mem->cv_nge++;
  if (retval != 0) return(CV_RTFUNC_FAIL);

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < cv_mem->cv_nrtfn; i++) {
    if (!cv_mem->cv_gactive[i] && SUNRabs(cv_mem->cv_ghi[i]) != ZERO) {
      cv_mem->cv_gactive[i] = SUNTRUE;
      cv_mem->cv_glo[i] = cv_mem->cv_ghi[i];
    }
  }
  return(CV_SUCCESS);
}

/*
 * cvRcheck2
 *
 * This routine checks for exact zeros of g at the last root found,
 * if the last return was a root.  It then checks for a close pair of
 * zeros (an error condition), and for a new root at a nearby point.
 * The array glo = g(tlo) at the left endpoint of the search interval
 * is adjusted if necessary to assure that all g_i are nonzero
 * there, before returning to do a root search in the interval.
 *
 * On entry, tlo = tretlast is the last value of tret returned by
 * CVode.  This may be the previous tn, the previous tout value,
 * or the last root location.
 *
 * This routine returns an int equal to:
 *     CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *     CLOSERT         = 3 if a close pair of zeros was found, or
 *     RTFOUND         = 1 if a new zero of g was found near tlo, or
 *     CV_SUCCESS      = 0 otherwise.
 */

int cvRcheck2_gpu2(CVodeMem cv_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  if (cv_mem->cv_irfnd == 0) return(CV_SUCCESS);

  (void) CVodeGetDky(cv_mem, cv_mem->cv_tlo, 0, cv_mem->cv_y);
  retval = cv_mem->cv_gfun(cv_mem->cv_tlo, cv_mem->cv_y,
                           cv_mem->cv_glo, cv_mem->cv_user_data);
  cv_mem->cv_nge++;
  if (retval != 0) return(CV_RTFUNC_FAIL);

  zroot = SUNFALSE;
  for (i = 0; i < cv_mem->cv_nrtfn; i++) cv_mem->cv_iroots[i] = 0;
  for (i = 0; i < cv_mem->cv_nrtfn; i++) {
    if (!cv_mem->cv_gactive[i]) continue;
    if (SUNRabs(cv_mem->cv_glo[i]) == ZERO) {
      zroot = SUNTRUE;
      cv_mem->cv_iroots[i] = 1;
    }
  }
  if (!zroot) return(CV_SUCCESS);

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  cv_mem->cv_ttol = (SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_h)) *
                    cv_mem->cv_uround * HUNDRED;
  smallh = (cv_mem->cv_h > ZERO) ? cv_mem->cv_ttol : -cv_mem->cv_ttol;
  tplus = cv_mem->cv_tlo + smallh;
  if ( (tplus - cv_mem->cv_tn)*cv_mem->cv_h >= ZERO) {
    hratio = smallh/cv_mem->cv_h;
    N_VLinearSum(ONE, cv_mem->cv_y, hratio, cv_mem->cv_zn[1], cv_mem->cv_y);
  } else {
    (void) CVodeGetDky(cv_mem, tplus, 0, cv_mem->cv_y);
  }
  retval = cv_mem->cv_gfun(tplus, cv_mem->cv_y,
                           cv_mem->cv_ghi, cv_mem->cv_user_data);
  cv_mem->cv_nge++;
  if (retval != 0) return(CV_RTFUNC_FAIL);

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = SUNFALSE;
  for (i = 0; i < cv_mem->cv_nrtfn; i++) {
    if (!cv_mem->cv_gactive[i]) continue;
    if (SUNRabs(cv_mem->cv_ghi[i]) == ZERO) {
      if (cv_mem->cv_iroots[i] == 1) return(CLOSERT);
      zroot = SUNTRUE;
      cv_mem->cv_iroots[i] = 1;
    } else {
      if (cv_mem->cv_iroots[i] == 1)
        cv_mem->cv_glo[i] = cv_mem->cv_ghi[i];
    }
  }
  if (zroot) return(RTFOUND);
  return(CV_SUCCESS);
}

/*
 * cvRcheck3
 *
 * This routine interfaces to cvRootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * This routine returns an int equal to:
 *     CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *     RTFOUND         = 1 if a root of g was found, or
 *     CV_SUCCESS      = 0 otherwise.
 */

int cvRcheck3_gpu2(CVodeMem cv_mem)
{
  int i, ier, retval;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (cv_mem->cv_taskc == CV_ONE_STEP) {
    cv_mem->cv_thi = cv_mem->cv_tn;
    N_VScale(ONE, cv_mem->cv_zn[0], cv_mem->cv_y);
  }
  if (cv_mem->cv_taskc == CV_NORMAL) {
    if ( (cv_mem->cv_toutc - cv_mem->cv_tn)*cv_mem->cv_h >= ZERO) {
      cv_mem->cv_thi = cv_mem->cv_tn;
      N_VScale(ONE, cv_mem->cv_zn[0], cv_mem->cv_y);
    } else {
      cv_mem->cv_thi = cv_mem->cv_toutc;
      (void) CVodeGetDky(cv_mem, cv_mem->cv_thi, 0, cv_mem->cv_y);
    }
  }

  /* Set ghi = g(thi) and call cvRootfind to search (tlo,thi) for roots. */
  retval = cv_mem->cv_gfun(cv_mem->cv_thi, cv_mem->cv_y,
                           cv_mem->cv_ghi, cv_mem->cv_user_data);
  cv_mem->cv_nge++;
  if (retval != 0) return(CV_RTFUNC_FAIL);

  cv_mem->cv_ttol = (SUNRabs(cv_mem->cv_tn) + SUNRabs(cv_mem->cv_h)) *
                    cv_mem->cv_uround * HUNDRED;
  ier = cvRootfind_gpu2(cv_mem);
  if (ier == CV_RTFUNC_FAIL) return(CV_RTFUNC_FAIL);
  for(i=0; i<cv_mem->cv_nrtfn; i++) {
    if(!cv_mem->cv_gactive[i] && cv_mem->cv_grout[i] != ZERO)
      cv_mem->cv_gactive[i] = SUNTRUE;
  }
  cv_mem->cv_tlo = cv_mem->cv_trout;
  for (i = 0; i < cv_mem->cv_nrtfn; i++)
    cv_mem->cv_glo[i] = cv_mem->cv_grout[i];

  /* If no root found, return CV_SUCCESS. */
  if (ier == CV_SUCCESS) return(CV_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) CVodeGetDky(cv_mem, cv_mem->cv_trout, 0, cv_mem->cv_y);
  return(RTFOUND);
}

/*
 * cvRootfind
 *
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 * Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 * Defined Output Points for Solutions of ODEs, Sandia National
 * Laboratory Report SAND80-0180, February 1980.
 *
 * This routine uses the following parameters for communication:
 *
 * nrtfn    = number of functions g_i, or number of components of
 *            the vector-valued function g(t).  Input only.
 *
 * gfun     = user-defined function for g(t).  Its form is
 *            (void) gfun(t, y, gt, user_data)
 *
 * rootdir  = in array specifying the direction of zero-crossings.
 *            If rootdir[i] > 0, search for roots of g_i only if
 *            g_i is increasing; if rootdir[i] < 0, search for
 *            roots of g_i only if g_i is decreasing; otherwise
 *            always search for roots of g_i.
 *
 * gactive  = array specifying whether a component of g should
 *            or should not be monitored. gactive[i] is initially
 *            set to SUNTRUE for all i=0,...,nrtfn-1, but it may be
 *            reset to SUNFALSE if at the first step g[i] is 0.0
 *            both at the I.C. and at a small perturbation of them.
 *            gactive[i] is then set back on SUNTRUE only after the
 *            corresponding g function moves away from 0.0.
 *
 * nge      = cumulative counter for gfun calls.
 *
 * ttol     = a convergence tolerance for trout.  Input only.
 *            When a root at trout is found, it is located only to
 *            within a tolerance of ttol.  Typically, ttol should
 *            be set to a value on the order of
 *               100 * UROUND * max (SUNRabs(tlo), SUNRabs(thi))
 *            where UROUND is the unit roundoff of the machine.
 *
 * tlo, thi = endpoints of the interval in which roots are sought.
 *            On input, these must be distinct, but tlo - thi may
 *            be of either sign.  The direction of integration is
 *            assumed to be from tlo to thi.  On return, tlo and thi
 *            are the endpoints of the final relevant interval.
 *
 * glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
 *            and g(thi) respectively.  Input and output.  On input,
 *            none of the glo[i] should be zero.
 *
 * trout    = root location, if a root was found, or thi if not.
 *            Output only.  If a root was found other than an exact
 *            zero of g, trout is the endpoint thi of the final
 *            interval bracketing the root, with size at most ttol.
 *
 * grout    = array of length nrtfn containing g(trout) on return.
 *
 * iroots   = int array of length nrtfn with root information.
 *            Output only.  If a root was found, iroots indicates
 *            which components g_i have a root at trout.  For
 *            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
 *            and g_i is increasing, iroots[i] = -1 if g_i has a
 *            root and g_i is decreasing, and iroots[i] = 0 if g_i
 *            has no roots or g_i varies in the direction opposite
 *            to that indicated by rootdir[i].
 *
 * This routine returns an int equal to:
 *      CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *      RTFOUND         = 1 if a root of g was found, or
 *      CV_SUCCESS      = 0 otherwise.
 */

int cvRootfind_gpu2(CVodeMem cv_mem)
{
  realtype alph, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, retval, imax, side, sideprev;
  booleantype zroot, sgnchg;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot = SUNFALSE;
  sgnchg = SUNFALSE;
  for (i = 0;  i < cv_mem->cv_nrtfn; i++) {
    if(!cv_mem->cv_gactive[i]) continue;
    if (SUNRabs(cv_mem->cv_ghi[i]) == ZERO) {
      if(cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) {
        zroot = SUNTRUE;
      }
    } else {
      if ( (cv_mem->cv_glo[i]*cv_mem->cv_ghi[i] < ZERO) &&
           (cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) ) {
        gfrac = SUNRabs(cv_mem->cv_ghi[i]/(cv_mem->cv_ghi[i] - cv_mem->cv_glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = SUNTRUE;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  /* If no sign change was found, reset trout and grout.  Then return
     CV_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */
  if (!sgnchg) {
    cv_mem->cv_trout = cv_mem->cv_thi;
    for (i = 0; i < cv_mem->cv_nrtfn; i++) cv_mem->cv_grout[i] = cv_mem->cv_ghi[i];
    if (!zroot) return(CV_SUCCESS);
    for (i = 0; i < cv_mem->cv_nrtfn; i++) {
      cv_mem->cv_iroots[i] = 0;
      if(!cv_mem->cv_gactive[i]) continue;
      if ( (SUNRabs(cv_mem->cv_ghi[i]) == ZERO) &&
           (cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) )
        cv_mem->cv_iroots[i] = cv_mem->cv_glo[i] > 0 ? -1 : 1;
    }
    return(RTFOUND);
  }

  /* Initialize alph to avoid compiler warning */
  alph = ONE;

  /* A sign change was found.  Loop to locate nearest root. */

  side = 0;  sideprev = -1;
  for(;;) {                                    /* Looping point */

    /* If interval size is already less than tolerance ttol, break. */
    if (SUNRabs(cv_mem->cv_thi - cv_mem->cv_tlo) <= cv_mem->cv_ttol) break;

    /* Set weight alph.
       On the first two passes, set alph = 1.  Thereafter, reset alph
       according to the side (low vs high) of the subinterval in which
       the sign change was found in the previous two passes.
       If the sides were opposite, set alph = 1.
       If the sides were the same, then double alph (if high side),
       or halve alph (if low side).
       The next guess tmid is the secant method value if alph = 1, but
       is closer to tlo if alph < 1, and closer to thi if alph > 1.    */

    if (sideprev == side) {
      alph = (side == 2) ? alph*TWO : alph*HALF;
    } else {
      alph = ONE;
    }

    /* Set next root approximation tmid and get g(tmid).
       If tmid is too close to tlo or thi, adjust it inward,
       by a fractional distance that is between 0.1 and 0.5.  */
    tmid = cv_mem->cv_thi - (cv_mem->cv_thi - cv_mem->cv_tlo) *
                            cv_mem->cv_ghi[imax] / (cv_mem->cv_ghi[imax] - alph*cv_mem->cv_glo[imax]);
    if (SUNRabs(tmid - cv_mem->cv_tlo) < HALF*cv_mem->cv_ttol) {
      fracint = SUNRabs(cv_mem->cv_thi - cv_mem->cv_tlo)/cv_mem->cv_ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = cv_mem->cv_tlo + fracsub*(cv_mem->cv_thi - cv_mem->cv_tlo);
    }
    if (SUNRabs(cv_mem->cv_thi - tmid) < HALF*cv_mem->cv_ttol) {
      fracint = SUNRabs(cv_mem->cv_thi - cv_mem->cv_tlo)/cv_mem->cv_ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = cv_mem->cv_thi - fracsub*(cv_mem->cv_thi - cv_mem->cv_tlo);
    }

    (void) CVodeGetDky(cv_mem, tmid, 0, cv_mem->cv_y);
    retval = cv_mem->cv_gfun(tmid, cv_mem->cv_y, cv_mem->cv_grout,
                             cv_mem->cv_user_data);
    cv_mem->cv_nge++;
    if (retval != 0) return(CV_RTFUNC_FAIL);

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */
    maxfrac = ZERO;
    zroot = SUNFALSE;
    sgnchg = SUNFALSE;
    sideprev = side;
    for (i = 0;  i < cv_mem->cv_nrtfn; i++) {
      if(!cv_mem->cv_gactive[i]) continue;
      if (SUNRabs(cv_mem->cv_grout[i]) == ZERO) {
        if(cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) zroot = SUNTRUE;
      } else {
        if ( (cv_mem->cv_glo[i]*cv_mem->cv_grout[i] < ZERO) &&
             (cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) ) {
          gfrac = SUNRabs(cv_mem->cv_grout[i]/(cv_mem->cv_grout[i] - cv_mem->cv_glo[i]));
          if (gfrac > maxfrac) {
            sgnchg = SUNTRUE;
            maxfrac = gfrac;
            imax = i;
          }
        }
      }
    }
    if (sgnchg) {
      /* Sign change found in (tlo,tmid); replace thi with tmid. */
      cv_mem->cv_thi = tmid;
      for (i = 0; i < cv_mem->cv_nrtfn; i++)
        cv_mem->cv_ghi[i] = cv_mem->cv_grout[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (SUNRabs(cv_mem->cv_thi - cv_mem->cv_tlo) <= cv_mem->cv_ttol) break;
      continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      cv_mem->cv_thi = tmid;
      for (i = 0; i < cv_mem->cv_nrtfn; i++)
        cv_mem->cv_ghi[i] = cv_mem->cv_grout[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    cv_mem->cv_tlo = tmid;
    for (i = 0; i < cv_mem->cv_nrtfn; i++)
      cv_mem->cv_glo[i] = cv_mem->cv_grout[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (SUNRabs(cv_mem->cv_thi - cv_mem->cv_tlo) <= cv_mem->cv_ttol) break;

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  cv_mem->cv_trout = cv_mem->cv_thi;
  for (i = 0; i < cv_mem->cv_nrtfn; i++) {
    cv_mem->cv_grout[i] = cv_mem->cv_ghi[i];
    cv_mem->cv_iroots[i] = 0;
    if(!cv_mem->cv_gactive[i]) continue;
    if ( (SUNRabs(cv_mem->cv_ghi[i]) == ZERO) &&
         (cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) )
      cv_mem->cv_iroots[i] = cv_mem->cv_glo[i] > 0 ? -1 : 1;
    if ( (cv_mem->cv_glo[i]*cv_mem->cv_ghi[i] < ZERO) &&
         (cv_mem->cv_rootdir[i]*cv_mem->cv_glo[i] <= ZERO) )
      cv_mem->cv_iroots[i] = cv_mem->cv_glo[i] > 0 ? -1 : 1;
  }
  return(RTFOUND);
}

//malloc and copy
void set_data_gpu2(CVodeMem cv_mem, SolverData *sd)
{
  /*itsolver *bicg = &(sd->bicg);
  ModelData *md = &(sd->model_data);
  CVDlsMem cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  SUNMatrix J = cvdls_mem->A;

  //Init GPU ODE solver variables
  bicg->nnz=SM_NNZ_S(J);
  bicg->nrows=SM_NP_S(J);
  bicg->mattype=1; //CSC
  bicg->A=(double*)SM_DATA_S(J);
  bicg->jA=(int*)SM_INDEXVALS_S(J);
  bicg->iA=(int*)SM_INDEXPTRS_S(J);
  bicg->dA=md->jac_gpu_data;//set itsolver gpu pointer to jac pointer initialized at camp
  bicg->dftemp=md->deriv_gpu_data;

  // Setting data
  bicg->maxIt=100;
  bicg->tolmax=cv_mem->cv_reltol; //Set to CAMP selected accuracy (1e-8) //1e-10;//1e-6
*/
}

/*
 * cvStep
 *
 * This routine performs one internal cvode step, from tn to tn + h.
 * It calls other routines to do all the work.
 *
 * The main operations done here are as follows:
 * - preliminary adjustments if a new step size was chosen;
 * - prediction of the Nordsieck history array zn at tn + h;
 * - setting of multistep method coefficients and test quantities;
 * - solution of the nonlinear system;
 * - testing the local error;
 * - updating zn and other state data if successful;
 * - resetting stepsize and order for the next step.
 * - if SLDET is on, check for stability, reduce order if necessary.
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 */
int cvStep_gpu2(SolverData *sd, CVodeMem cv_mem)
{
  itsolver *bicg = &(sd->bicg);
  realtype saved_t, dsm;
  int ncf, nef;
  int nflag, kflag, eflag;

  double *ewt = NV_DATA_S(cv_mem->cv_ewt);

  cudaMemcpy(bicg->dewt,ewt,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);

  saved_t = cv_mem->cv_tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  if ((cv_mem->cv_nst > 0) && (cv_mem->cv_hprime != cv_mem->cv_h))
    cvAdjustParams_gpu2(cv_mem);

  /* Looping point for attempts to take a step */
  for(;;) {

    cvPredict_gpu2(cv_mem);

    cvSet_gpu2(cv_mem);

    //nflag = cvNls(cv_mem, nflag);
#ifdef PMC_DEBUG_GPU
    clock_t start=clock();
#endif

    nflag = cvNlsNewton_gpu2(sd, cv_mem, nflag);

#ifdef PMC_DEBUG_GPU
    timeNewtonIt+= clock() - start;
    counterNewtonIt++;
#endif
    kflag = cvHandleNFlag_gpu2(cv_mem, &nflag, saved_t, &ncf);

    /* Go back in loop if we need to predict again (nflag=PREV_CONV_FAIL)*/
    if (kflag == PREDICT_AGAIN) continue;

    /* Return if nonlinear solve failed and recovery not possible. */
    if (kflag != DO_ERROR_TEST) return(kflag);

    /* Perform error test (nflag=CV_SUCCESS) */
    eflag = cvDoErrorTest_gpu2(cv_mem, &nflag, saved_t, &nef, &dsm);

    /* Go back in loop if we need to predict again (nflag=PREV_ERR_FAIL) */
    if (eflag == TRY_AGAIN)  continue;

    /* Return if error test failed and recovery not possible. */
    if (eflag != CV_SUCCESS) return(eflag);

    /* Error test passed (eflag=CV_SUCCESS), break from loop */
    break;

  }

  /* Nonlinear system solve and error test were both successful.
     Update data, and consider change of step and/or order.       */

  cvCompleteStep_gpu2(cv_mem);

  cvPrepareNextStep_gpu2(cv_mem, dsm);//use tq calculated in cvset and tempv calc in cvnewton

  /* If Stablilty Limit Detection is turned on, call stability limit
     detection routine for possible order reduction. */

  if (cv_mem->cv_sldeton) cvBDFStab_gpu2(cv_mem);

  cv_mem->cv_etamax = (cv_mem->cv_nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /*  Finally, we rescale the acor array to be the
      estimated local error vector. */

  N_VScale(cv_mem->cv_tq[2], cv_mem->cv_acor, cv_mem->cv_acor);
  return(CV_SUCCESS);

}

/*
 * cvAdjustParams
 *
 * This routine is called when a change in step size was decided upon,
 * and it handles the required adjustments to the history array zn.
 * If there is to be a change in order, we call cvAdjustOrder and reset
 * q, L = q+1, and qwait.  Then in any case, we call cvRescale, which
 * resets h and rescales the Nordsieck array.
 */

void cvAdjustParams_gpu2(CVodeMem cv_mem)
{
  if (cv_mem->cv_qprime != cv_mem->cv_q) {
    //cvAdjustOrder(cv_mem, cv_mem->cv_qprime-cv_mem->cv_q);

    int deltaq = cv_mem->cv_qprime-cv_mem->cv_q;
    switch(deltaq) {
      case 1:
        cvIncreaseBDF_gpu2(cv_mem);
        break;
      case -1:
        cvDecreaseBDF_gpu2(cv_mem);
        break;
    }

    cv_mem->cv_q = cv_mem->cv_qprime;
    cv_mem->cv_L = cv_mem->cv_q+1;
    cv_mem->cv_qwait = cv_mem->cv_L;
  }
  cvRescale_gpu2(cv_mem);
}

/*
 * cvIncreaseBDF
 *
 * This routine adjusts the history array on an increase in the
 * order q in the case that lmm == CV_BDF.
 * A new column zn[q+1] is set equal to a multiple of the saved
 * vector (= acor) in zn[indx_acor].  Then each zn[j] is adjusted by
 * a multiple of zn[q+1].  The coefficients in the adjustment are the
 * coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 * where xi_j = [t_n - t_(n-j)]/h.
 */

void cvIncreaseBDF_gpu2(CVodeMem cv_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;

  for (i=0; i <= cv_mem->cv_qmax; i++) cv_mem->cv_l[i] = ZERO;
  cv_mem->cv_l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = cv_mem->cv_hscale;
  if (cv_mem->cv_q > 1) {
    for (j=1; j < cv_mem->cv_q; j++) {
      hsum += cv_mem->cv_tau[j+1];
      xi = hsum / cv_mem->cv_hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--)
        cv_mem->cv_l[i] = cv_mem->cv_l[i]*xiold + cv_mem->cv_l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, cv_mem->cv_zn[cv_mem->cv_indx_acor],
           cv_mem->cv_zn[cv_mem->cv_L]);
  for (j=2; j <= cv_mem->cv_q; j++)
    N_VLinearSum(cv_mem->cv_l[j], cv_mem->cv_zn[cv_mem->cv_L], ONE,
                 cv_mem->cv_zn[j], cv_mem->cv_zn[j]);
}

/*
 * cvDecreaseBDF
 *
 * This routine adjusts the history array on a decrease in the
 * order q in the case that lmm == CV_BDF.
 * Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 * in the adjustment are the coefficients of the polynomial
 *   x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.
 */

void cvDecreaseBDF_gpu2(CVodeMem cv_mem)
{
  realtype hsum, xi;
  int i, j;

  for (i=0; i <= cv_mem->cv_qmax; i++) cv_mem->cv_l[i] = ZERO;
  cv_mem->cv_l[2] = ONE;
  hsum = ZERO;
  for (j=1; j <= cv_mem->cv_q-2; j++) {
    hsum += cv_mem->cv_tau[j];
    xi = hsum /cv_mem->cv_hscale;
    for (i=j+2; i >= 2; i--)
      cv_mem->cv_l[i] = cv_mem->cv_l[i]*xi + cv_mem->cv_l[i-1];
  }

  for (j=2; j < cv_mem->cv_q; j++)
    N_VLinearSum(-cv_mem->cv_l[j], cv_mem->cv_zn[cv_mem->cv_q],
                 ONE, cv_mem->cv_zn[j], cv_mem->cv_zn[j]);
}

/*
 * cvRescale
 *
 * This routine rescales the Nordsieck array by multiplying the
 * jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
 * h is rescaled by eta, and hscale is reset to h.
 */

void cvRescale_gpu2(CVodeMem cv_mem)
{
  int j;
  realtype factor;

  factor = cv_mem->cv_eta;
  for (j=1; j <= cv_mem->cv_q; j++) {
    N_VScale(factor, cv_mem->cv_zn[j], cv_mem->cv_zn[j]);
    factor *= cv_mem->cv_eta;
  }
  cv_mem->cv_h = cv_mem->cv_hscale * cv_mem->cv_eta;
  cv_mem->cv_next_h = cv_mem->cv_h;
  cv_mem->cv_hscale = cv_mem->cv_h;
  cv_mem->cv_nscon = 0;
}

/*
 * cvPredict
 *
 * This routine advances tn by the tentative step size h, and computes
 * the predicted array z_n(0), which is overwritten on zn.  The
 * prediction of zn is done by repeated additions.
 * If tstop is enabled, it is possible for tn + h to be past tstop by roundoff,
 * and in that case, we reset tn (after incrementing by h) to tstop.
 */

void cvPredict_gpu2(CVodeMem cv_mem)
{
  int j, k;

  cv_mem->cv_tn += cv_mem->cv_h;
  if (cv_mem->cv_tstopset) {
    if ((cv_mem->cv_tn - cv_mem->cv_tstop)*cv_mem->cv_h > ZERO)
      cv_mem->cv_tn = cv_mem->cv_tstop;
  }
  N_VScale(ONE, cv_mem->cv_zn[0], cv_mem->cv_last_yn);
  for (k = 1; k <= cv_mem->cv_q; k++)
    for (j = cv_mem->cv_q; j >= k; j--)
      N_VLinearSum(ONE, cv_mem->cv_zn[j-1], ONE,
                   cv_mem->cv_zn[j], cv_mem->cv_zn[j-1]);

}

/*
 * cvSet
 *
 * This routine is a high level routine which calls cvSetAdams or
 * cvSetBDF to set the polynomial l, the test quantity array tq,
 * and the related variables  rl1, gamma, and gamrat.
 *
 * The array tq is loaded with constants used in the control of estimated
 * local errors and in the nonlinear convergence test.  Specifically, while
 * running at order q, the components of tq are as follows:
 *   tq[1] = a coefficient used to get the est. local error at order q-1
 *   tq[2] = a coefficient used to get the est. local error at order q
 *   tq[3] = a coefficient used to get the est. local error at order q+1
 *   tq[4] = constant used in nonlinear iteration convergence test
 *   tq[5] = coefficient used to get the order q+2 derivative vector used in
 *           the est. local error at order q+1
 */

void cvSet_gpu2(CVodeMem cv_mem)
{

  cvSetBDF_gpu2(cv_mem);
  cv_mem->cv_rl1 = ONE / cv_mem->cv_l[1];
  cv_mem->cv_gamma = cv_mem->cv_h * cv_mem->cv_rl1;
  if (cv_mem->cv_nst == 0) cv_mem->cv_gammap = cv_mem->cv_gamma;
  cv_mem->cv_gamrat = (cv_mem->cv_nst > 0) ?
                      cv_mem->cv_gamma / cv_mem->cv_gammap : ONE;  /* protect x / x != 1.0 */
}

/*
 * cvSetBDF
 *
 * This routine computes the coefficients l and tq in the case
 * lmm == CV_BDF.  cvSetBDF calls cvSetTqBDF to set the test
 * quantity array tq.
 *
 * The components of the array l are the coefficients of a
 * polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
 *                                 q-1
 * Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
 *                                 i=1
 *  xi_i = [t_n - t_(n-i)] / h.
 *
 * The array tq is set to test quantities used in the convergence
 * test, the error test, and the selection of h at a new order.
 */

void cvSetBDF_gpu2(CVodeMem cv_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;

  cv_mem->cv_l[0] = cv_mem->cv_l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= cv_mem->cv_q; i++) cv_mem->cv_l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = cv_mem->cv_h;
  if (cv_mem->cv_q > 1) {
    for (j=2; j < cv_mem->cv_q; j++) {
      hsum += cv_mem->cv_tau[j-1];
      xi_inv = cv_mem->cv_h / hsum;
      alpha0 -= ONE / j;
      for (i=j; i >= 1; i--) cv_mem->cv_l[i] += cv_mem->cv_l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }

    /* j = q */
    alpha0 -= ONE / cv_mem->cv_q;
    xistar_inv = -cv_mem->cv_l[1] - alpha0;
    hsum += cv_mem->cv_tau[cv_mem->cv_q-1];
    xi_inv = cv_mem->cv_h / hsum;
    alpha0_hat = -cv_mem->cv_l[1] - xi_inv;
    for (i=cv_mem->cv_q; i >= 1; i--)
      cv_mem->cv_l[i] += cv_mem->cv_l[i-1]*xistar_inv;
  }

  cvSetTqBDF_gpu2(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/*
 * cvSetTqBDF
 *
 * This routine sets the test quantity array tq in the case
 * lmm == CV_BDF.
 */

void cvSetTqBDF_gpu2(CVodeMem cv_mem, realtype hsum, realtype alpha0,
                       realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, Cpinv, Cppinv;

  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + cv_mem->cv_q * A1;
  cv_mem->cv_tq[2] = SUNRabs(A1 / (alpha0 * A2));
  cv_mem->cv_tq[5] = SUNRabs(A2 * xistar_inv / (cv_mem->cv_l[cv_mem->cv_q] * xi_inv));
  if (cv_mem->cv_qwait == 1) {
    if (cv_mem->cv_q > 1) {
      C = xistar_inv / cv_mem->cv_l[cv_mem->cv_q];
      A3 = alpha0 + ONE / cv_mem->cv_q;
      A4 = alpha0_hat + xi_inv;
      Cpinv = (ONE - A4 + A3) / A3;
      cv_mem->cv_tq[1] = SUNRabs(C * Cpinv);
    }
    else cv_mem->cv_tq[1] = ONE;
    hsum += cv_mem->cv_tau[cv_mem->cv_q];
    xi_inv = cv_mem->cv_h / hsum;
    A5 = alpha0 - (ONE / (cv_mem->cv_q+1));
    A6 = alpha0_hat - xi_inv;
    Cppinv = (ONE - A6 + A5) / A2;
    cv_mem->cv_tq[3] = SUNRabs(Cppinv / (xi_inv * (cv_mem->cv_q+2) * A5));
  }
  cv_mem->cv_tq[4] = cv_mem->cv_nlscoef / cv_mem->cv_tq[2];
}


/*
 * cvHandleNFlag
 *
 * This routine takes action on the return value nflag = *nflagPtr
 * returned by cvNls, as follows:
 *
 * If cvNls succeeded in solving the nonlinear system, then
 * cvHandleNFlag returns the constant DO_ERROR_TEST, which tells cvStep
 * to perform the error test.
 *
 * If the nonlinear system was not solved successfully, then ncfn and
 * ncf = *ncfPtr are incremented and Nordsieck array zn is restored.
 *
 * If the solution of the nonlinear system failed due to an
 * unrecoverable failure by setup, we return the value CV_LSETUP_FAIL.
 *
 * If it failed due to an unrecoverable failure in solve, then we return
 * the value CV_LSOLVE_FAIL.
 *
 * If it failed due to an unrecoverable failure in rhs, then we return
 * the value CV_RHSFUNC_FAIL.
 *
 * Otherwise, a recoverable failure occurred when solving the
 * nonlinear system (cvNls returned nflag == CONV_FAIL or RHSFUNC_RECVR).
 * In this case, if ncf is now equal to maxncf or |h| = hmin,
 * we return the value CV_CONV_FAILURE (if nflag=CONV_FAIL) or
 * CV_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
 * If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
 * PREDICT_AGAIN, telling cvStep to reattempt the step.
 *
 */

int cvHandleNFlag_gpu2(CVodeMem cv_mem, int *nflagPtr, realtype saved_t,
                         int *ncfPtr)
{
  int nflag;

  nflag = *nflagPtr;

  if (nflag == CV_SUCCESS) return(DO_ERROR_TEST);

  /* The nonlinear soln. failed; increment ncfn and restore zn */
  cv_mem->cv_ncfn++;
  cvRestore_gpu2(cv_mem, saved_t);

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == CV_LSETUP_FAIL)  return(CV_LSETUP_FAIL);
  if (nflag == CV_LSOLVE_FAIL)  return(CV_LSOLVE_FAIL);
  if (nflag == CV_RHSFUNC_FAIL) return(CV_RHSFUNC_FAIL);

  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */

  (*ncfPtr)++;
  cv_mem->cv_etamax = ONE;

  /* If we had maxncf failures or |h| = hmin,
     return CV_CONV_FAILURE or CV_REPTD_RHSFUNC_ERR. */

  if ((SUNRabs(cv_mem->cv_h) <= cv_mem->cv_hmin*ONEPSM) ||
      (*ncfPtr == cv_mem->cv_maxncf)) {
    if (nflag == CONV_FAIL)     return(CV_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(CV_REPTD_RHSFUNC_ERR);
  }

  /* Reduce step size; return to reattempt the step */

  cv_mem->cv_eta = SUNMAX(ETACF, cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h));
  *nflagPtr = PREV_CONV_FAIL;
  cvRescale_gpu2(cv_mem);

  return(PREDICT_AGAIN);
}

/*
 * cvRestore
 *
 * This routine restores the value of tn to saved_t and undoes the
 * prediction.  After execution of cvRestore, the Nordsieck array zn has
 * the same values as before the call to cvPredict.
 */

void cvRestore_gpu2(CVodeMem cv_mem, realtype saved_t)
{
  int j, k;

  cv_mem->cv_tn = saved_t;
  for (k = 1; k <= cv_mem->cv_q; k++)
    for (j = cv_mem->cv_q; j >= k; j--)
      N_VLinearSum(ONE, cv_mem->cv_zn[j-1], -ONE,
                   cv_mem->cv_zn[j], cv_mem->cv_zn[j-1]);
  N_VScale(ONE, cv_mem->cv_last_yn, cv_mem->cv_zn[0]);
}

/*
 * cvDoErrorTest
 *
 * This routine performs the local error test.
 * The weighted local error norm dsm is loaded into *dsmPtr, and
 * the test dsm ?<= 1 is made.
 *
 * If the test passes, cvDoErrorTest returns CV_SUCCESS.
 *
 * If the test fails, we undo the step just taken (call cvRestore) and
 *
 *   - if maxnef error test failures have occurred or if SUNRabs(h) = hmin,
 *     we return CV_ERR_FAILURE.
 *
 *   - if more than MXNEF1 error test failures have occurred, an order
 *     reduction is forced. If already at order 1, restart by reloading
 *     zn from scratch. If f() fails we return either CV_RHSFUNC_FAIL
 *     or CV_UNREC_RHSFUNC_ERR (no recovery is possible at this stage).
 *
 *   - otherwise, set *nflagPtr to PREV_ERR_FAIL, and return TRY_AGAIN.
 *
 */
booleantype cvDoErrorTest_gpu2(CVodeMem cv_mem, int *nflagPtr,
                                 realtype saved_t, int *nefPtr, realtype *dsmPtr)
{
  realtype dsm;
  realtype min_val;
  int retval;

  /* Find the minimum concentration and if it's small and negative, make it
   * positive */
  N_VLinearSum(cv_mem->cv_l[0], cv_mem->cv_acor, ONE, cv_mem->cv_zn[0],
               cv_mem->cv_ftemp);
  min_val = N_VMin(cv_mem->cv_ftemp);
  if (min_val < ZERO && min_val > -PMC_TINY) {
    N_VAbs(cv_mem->cv_ftemp, cv_mem->cv_ftemp);
    N_VLinearSum(-cv_mem->cv_l[0], cv_mem->cv_acor, ONE, cv_mem->cv_ftemp,
                 cv_mem->cv_zn[0]);
    min_val = ZERO;
  }

  dsm = cv_mem->cv_acnrm * cv_mem->cv_tq[2];

  /* If est. local error norm dsm passes test and there are no negative values,
   * return CV_SUCCESS */
  *dsmPtr = dsm;
  if (dsm <= ONE && min_val >= ZERO) return(CV_SUCCESS);

  /* Test failed; increment counters, set nflag, and restore zn array */
  (*nefPtr)++;
  cv_mem->cv_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  cvRestore_gpu2(cv_mem, saved_t);

  /* At maxnef failures or |h| = hmin, return CV_ERR_FAILURE */
  if ((SUNRabs(cv_mem->cv_h) <= cv_mem->cv_hmin*ONEPSM) ||
      (*nefPtr == cv_mem->cv_maxnef)) return(CV_ERR_FAILURE);

  /* Set etamax = 1 to prevent step size increase at end of this step */
  cv_mem->cv_etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    cv_mem->cv_eta = ONE / (SUNRpowerR(BIAS2*dsm,ONE/cv_mem->cv_L) + ADDON);
    cv_mem->cv_eta = SUNMAX(ETAMIN, SUNMAX(cv_mem->cv_eta,
                                           cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h)));
    if (*nefPtr >= SMALL_NEF) cv_mem->cv_eta = SUNMIN(cv_mem->cv_eta, ETAMXF);
    cvRescale_gpu2(cv_mem);
    return(TRY_AGAIN);
  }

  /* After MXNEF1 failures, force an order reduction and retry step */
  if (cv_mem->cv_q > 1) {
    cv_mem->cv_eta = SUNMAX(ETAMIN, cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h));

    //cvAdjustOrder_gpu2(cv_mem,-1);
    cvDecreaseBDF_gpu2(cv_mem);

    cv_mem->cv_L = cv_mem->cv_q;
    cv_mem->cv_q--;
    cv_mem->cv_qwait = cv_mem->cv_L;
    cvRescale_gpu2(cv_mem);
    return(TRY_AGAIN);
  }

  /* If already at order 1, restart: reload zn from scratch */

  cv_mem->cv_eta = SUNMAX(ETAMIN, cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h));
  cv_mem->cv_h *= cv_mem->cv_eta;
  cv_mem->cv_next_h = cv_mem->cv_h;
  cv_mem->cv_hscale = cv_mem->cv_h;
  cv_mem->cv_qwait = LONG_WAIT;
  cv_mem->cv_nscon = 0;


  //retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_zn[0],
  //                      cv_mem->cv_tempv, cv_mem->cv_user_data);
  retval = f(cv_mem->cv_tn, cv_mem->cv_zn[0],cv_mem->cv_tempv, cv_mem->cv_user_data);

  cv_mem->cv_nfe++;
  if (retval < 0)  return(CV_RHSFUNC_FAIL);
  if (retval > 0)  return(CV_UNREC_RHSFUNC_ERR);

  N_VScale(cv_mem->cv_h, cv_mem->cv_tempv, cv_mem->cv_zn[1]);

  return(TRY_AGAIN);
}

/*
 * -----------------------------------------------------------------
 * Functions called after succesful step
 * -----------------------------------------------------------------
 */

/*
 * cvCompleteStep
 *
 * This routine performs various update operations when the solution
 * to the nonlinear system has passed the local error test.
 * We increment the step counter nst, record the values hu and qu,
 * update the tau array, and apply the corrections to the zn array.
 * The tau[i] are the last q values of h, with tau[1] the most recent.
 * The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 * we save acor and cv_mem->cv_tq[5] for a possible order increase.
 */
void cvCompleteStep_gpu2(CVodeMem cv_mem)
{
  int i, j;

  cv_mem->cv_nst++;
  cv_mem->cv_nscon++;
  cv_mem->cv_hu = cv_mem->cv_h;
  cv_mem->cv_qu = cv_mem->cv_q;

  for (i=cv_mem->cv_q; i >= 2; i--)  cv_mem->cv_tau[i] = cv_mem->cv_tau[i-1];
  if ((cv_mem->cv_q==1) && (cv_mem->cv_nst > 1))
    cv_mem->cv_tau[2] = cv_mem->cv_tau[1];
  cv_mem->cv_tau[1] = cv_mem->cv_h;

  /* Apply correction to column j of zn: l_j * Delta_n */
  for (j=0; j <= cv_mem->cv_q; j++)
    N_VLinearSum(cv_mem->cv_l[j], cv_mem->cv_acor, ONE,
                 cv_mem->cv_zn[j], cv_mem->cv_zn[j]);
  cv_mem->cv_qwait--;
  if ((cv_mem->cv_qwait == 1) && (cv_mem->cv_q != cv_mem->cv_qmax)) {
    N_VScale(ONE, cv_mem->cv_acor, cv_mem->cv_zn[cv_mem->cv_qmax]);
    cv_mem->cv_saved_tq5 = cv_mem->cv_tq[5];
    cv_mem->cv_indx_acor = cv_mem->cv_qmax;
  }
}


/*
 * cvPrepareNextStep
 *
 * This routine handles the setting of stepsize and order for the
 * next step -- hprime and qprime.  Along with hprime, it sets the
 * ratio eta = hprime/h.  It also updates other state variables
 * related to a change of step size or order.
 */
void cvPrepareNextStep_gpu2(CVodeMem cv_mem, realtype dsm)
{
  /* If etamax = 1, defer step size or order changes */
  if (cv_mem->cv_etamax == ONE) {
    cv_mem->cv_qwait = SUNMAX(cv_mem->cv_qwait, 2);
    cv_mem->cv_qprime = cv_mem->cv_q;
    cv_mem->cv_hprime = cv_mem->cv_h;
    cv_mem->cv_eta = ONE;
    return;
  }

  /* etaq is the ratio of new to old h at the current order */
  cv_mem->cv_etaq = ONE /(SUNRpowerR(BIAS2*dsm,ONE/cv_mem->cv_L) + ADDON);

  /* If no order change, adjust eta and acor in cvSetEta and return */
  if (cv_mem->cv_qwait != 0) {
    cv_mem->cv_eta = cv_mem->cv_etaq;
    cv_mem->cv_qprime = cv_mem->cv_q;
    cvSetEta_gpu2(cv_mem);
    return;
  }

  /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are
     the ratios of new to old h at orders q-1 and q+1, respectively.
     cvChooseEta selects the largest; cvSetEta adjusts eta and acor */
  cv_mem->cv_qwait = 2;

  //cv_mem->cv_etaqm1 = cvComputeEtaqm1_gpu2(cv_mem);
  //compute cv_etaqm1
  realtype ddn;
  cv_mem->cv_etaqm1 = ZERO;
  if (cv_mem->cv_q > 1) {
    ddn = N_VWrmsNorm(cv_mem->cv_zn[cv_mem->cv_q], cv_mem->cv_ewt) * cv_mem->cv_tq[1];
    cv_mem->cv_etaqm1 = ONE/(SUNRpowerR(BIAS1*ddn, ONE/cv_mem->cv_q) + ADDON);
  }

  //cv_mem->cv_etaqp1 = cvComputeEtaqp1_gpu2(cv_mem);
  //compute cv_etaqp1
  realtype dup, cquot;
  cv_mem->cv_etaqp1 = ZERO;
  if (cv_mem->cv_q != cv_mem->cv_qmax && cv_mem->cv_saved_tq5 != ZERO) {
    //if (cv_mem->cv_saved_tq5 != ZERO) return(cv_mem->cv_etaqp1);
    cquot = (cv_mem->cv_tq[5] / cv_mem->cv_saved_tq5) *
            SUNRpowerI(cv_mem->cv_h/cv_mem->cv_tau[2], cv_mem->cv_L);
    N_VLinearSum(-cquot, cv_mem->cv_zn[cv_mem->cv_qmax], ONE,
                 cv_mem->cv_acor, cv_mem->cv_tempv);
    dup = N_VWrmsNorm(cv_mem->cv_tempv, cv_mem->cv_ewt) * cv_mem->cv_tq[3];
    cv_mem->cv_etaqp1 = ONE / (SUNRpowerR(BIAS3*dup, ONE/(cv_mem->cv_L+1)) + ADDON);
  }

  cvChooseEta_gpu2(cv_mem);
  cvSetEta_gpu2(cv_mem);
}

/*
 * cvSetEta
 *
 * This routine adjusts the value of eta according to the various
 * heuristic limits and the optional input hmax.
 */

void cvSetEta_gpu2(CVodeMem cv_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (cv_mem->cv_eta < THRESH) {
    cv_mem->cv_eta = ONE;
    cv_mem->cv_hprime = cv_mem->cv_h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    cv_mem->cv_eta = SUNMIN(cv_mem->cv_eta, cv_mem->cv_etamax);
    cv_mem->cv_eta /= SUNMAX(ONE, SUNRabs(cv_mem->cv_h)*cv_mem->cv_hmax_inv*cv_mem->cv_eta);
    cv_mem->cv_hprime = cv_mem->cv_h * cv_mem->cv_eta;
    if (cv_mem->cv_qprime < cv_mem->cv_q) cv_mem->cv_nscon = 0;
  }
}

/*
 * cvChooseEta
 * Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 * q - 1, q, or q + 1, respectively), this routine chooses the
 * maximum eta value, sets eta to that value, and sets qprime to the
 * corresponding value of q.  If there is a tie, the preference
 * order is to (1) keep the same order, then (2) decrease the order,
 * and finally (3) increase the order.  If the maximum eta value
 * is below the threshhold THRESH, the order is kept unchanged and
 * eta is set to 1.
 */
void cvChooseEta_gpu2(CVodeMem cv_mem)
{
  realtype etam;

  etam = SUNMAX(cv_mem->cv_etaqm1, SUNMAX(cv_mem->cv_etaq, cv_mem->cv_etaqp1));

  if (etam < THRESH) {
    cv_mem->cv_eta = ONE;
    cv_mem->cv_qprime = cv_mem->cv_q;
    return;
  }

  if (etam == cv_mem->cv_etaq) {

    cv_mem->cv_eta = cv_mem->cv_etaq;
    cv_mem->cv_qprime = cv_mem->cv_q;

  } else if (etam == cv_mem->cv_etaqm1) {

    cv_mem->cv_eta = cv_mem->cv_etaqm1;
    cv_mem->cv_qprime = cv_mem->cv_q - 1;

  } else {

    cv_mem->cv_eta = cv_mem->cv_etaqp1;
    cv_mem->cv_qprime = cv_mem->cv_q + 1;

    if (cv_mem->cv_lmm == CV_BDF) {
      /*
       * Store Delta_n in zn[qmax] to be used in order increase
       *
       * This happens at the last step of order q before an increase
       * to order q+1, so it represents Delta_n in the ELTE at q+1
       */

      N_VScale(ONE, cv_mem->cv_acor, cv_mem->cv_zn[cv_mem->cv_qmax]);

    }
  }
}

/*
 * -----------------------------------------------------------------
 * Functions for BDF Stability Limit Detection
 * -----------------------------------------------------------------
 */

/*
 * cvBDFStab
 *
 * This routine handles the BDF Stability Limit Detection Algorithm
 * STALD.  It is called if lmm = CV_BDF and the SLDET option is on.
 * If the order is 3 or more, the required norm data is saved.
 * If a decision to reduce order has not already been made, and
 * enough data has been saved, cvSLdet is called.  If it signals
 * a stability limit violation, the order is reduced, and the step
 * size is reset accordingly.
 */

void cvBDFStab_gpu2(CVodeMem cv_mem)
{
  int i,k, ldflag, factorial;
  realtype sq, sqm1, sqm2;

  /* If order is 3 or greater, then save scaled derivative data,
     push old data down in i, then add current values to top.    */

  if (cv_mem->cv_q >= 3) {
    for (k = 1; k <= 3; k++)
      for (i = 5; i >= 2; i--)
        cv_mem->cv_ssdat[i][k] = cv_mem->cv_ssdat[i-1][k];
    factorial = 1;
    for (i = 1; i <= cv_mem->cv_q-1; i++) factorial *= i;
    sq = factorial * cv_mem->cv_q * (cv_mem->cv_q+1) *
         cv_mem->cv_acnrm / SUNMAX(cv_mem->cv_tq[5],TINY);
    sqm1 = factorial * cv_mem->cv_q *
           N_VWrmsNorm(cv_mem->cv_zn[cv_mem->cv_q], cv_mem->cv_ewt);
    sqm2 = factorial * N_VWrmsNorm(cv_mem->cv_zn[cv_mem->cv_q-1], cv_mem->cv_ewt);
    cv_mem->cv_ssdat[1][1] = sqm2*sqm2;
    cv_mem->cv_ssdat[1][2] = sqm1*sqm1;
    cv_mem->cv_ssdat[1][3] = sq*sq;
  }


  if (cv_mem->cv_qprime >= cv_mem->cv_q) {

    /* If order is 3 or greater, and enough ssdat has been saved,
       nscon >= q+5, then call stability limit detection routine.  */

    if ( (cv_mem->cv_q >= 3) && (cv_mem->cv_nscon >= cv_mem->cv_q+5) ) {
      ldflag = cvSLdet_gpu2(cv_mem);
      if (ldflag > 3) {
        /* A stability limit violation is indicated by
           a return flag of 4, 5, or 6.
           Reduce new order.                     */
        cv_mem->cv_qprime = cv_mem->cv_q-1;
        cv_mem->cv_eta = cv_mem->cv_etaqm1;
        cv_mem->cv_eta = SUNMIN(cv_mem->cv_eta,cv_mem->cv_etamax);
        cv_mem->cv_eta = cv_mem->cv_eta /
                         SUNMAX(ONE,SUNRabs(cv_mem->cv_h)*cv_mem->cv_hmax_inv*cv_mem->cv_eta);
        cv_mem->cv_hprime = cv_mem->cv_h*cv_mem->cv_eta;
        cv_mem->cv_nor = cv_mem->cv_nor + 1;
      }
    }
  }
  else {
    /* Otherwise, let order increase happen, and
       reset stability limit counter, nscon.     */
    cv_mem->cv_nscon = 0;
  }
}

/*
 * cvSLdet
 *
 * This routine detects stability limitation using stored scaled
 * derivatives data. cvSLdet returns the magnitude of the
 * dominate characteristic root, rr. The presence of a stability
 * limit is indicated by rr > "something a little less then 1.0",
 * and a positive kflag. This routine should only be called if
 * order is greater than or equal to 3, and data has been collected
 * for 5 time steps.
 *
 * Returned values:
 *    kflag = 1 -> Found stable characteristic root, normal matrix case
 *    kflag = 2 -> Found stable characteristic root, quartic solution
 *    kflag = 3 -> Found stable characteristic root, quartic solution,
 *                 with Newton correction
 *    kflag = 4 -> Found stability violation, normal matrix case
 *    kflag = 5 -> Found stability violation, quartic solution
 *    kflag = 6 -> Found stability violation, quartic solution,
 *                 with Newton correction
 *
 *    kflag < 0 -> No stability limitation,
 *                 or could not compute limitation.
 *
 *    kflag = -1 -> Min/max ratio of ssdat too small.
 *    kflag = -2 -> For normal matrix case, vmax > vrrt2*vrrt2
 *    kflag = -3 -> For normal matrix case, The three ratios
 *                  are inconsistent.
 *    kflag = -4 -> Small coefficient prevents elimination of quartics.
 *    kflag = -5 -> R value from quartics not consistent.
 *    kflag = -6 -> No corrected root passes test on qk values
 *    kflag = -7 -> Trouble solving for sigsq.
 *    kflag = -8 -> Trouble solving for B, or R via B.
 *    kflag = -9 -> R via sigsq[k] disagrees with R from data.
 */

int cvSLdet_gpu2(CVodeMem cv_mem)
{
  int i, k, j, it, kmin = 0, kflag = 0;
  realtype rat[5][4], rav[4], qkr[4], sigsq[4], smax[4], ssmax[4];
  realtype drr[4], rrc[4],sqmx[4], qjk[4][4], vrat[5], qc[6][4], qco[6][4];
  realtype rr, rrcut, vrrtol, vrrt2, sqtol, rrtol;
  realtype smink, smaxk, sumrat, sumrsq, vmin, vmax, drrmax, adrr;
  realtype tem, sqmax, saqk, qp, s, sqmaxk, saqj, sqmin;
  realtype rsa, rsb, rsc, rsd, rd1a, rd1b, rd1c;
  realtype rd2a, rd2b, rd3a, cest1, corr1;
  realtype ratp, ratm, qfac1, qfac2, bb, rrb;

  /* The following are cutoffs and tolerances used by this routine */

  rrcut  = RCONST(0.98);
  vrrtol = RCONST(1.0e-4);
  vrrt2  = RCONST(5.0e-4);
  sqtol  = RCONST(1.0e-3);
  rrtol  = RCONST(1.0e-2);

  rr = ZERO;

  /*  Index k corresponds to the degree of the interpolating polynomial. */
  /*      k = 1 -> q-1          */
  /*      k = 2 -> q            */
  /*      k = 3 -> q+1          */

  /*  Index i is a backward-in-time index, i = 1 -> current time, */
  /*      i = 2 -> previous step, etc    */

  /* get maxima, minima, and variances, and form quartic coefficients  */

  for (k=1; k<=3; k++) {
    smink = cv_mem->cv_ssdat[1][k];
    smaxk = ZERO;

    for (i=1; i<=5; i++) {
      smink = SUNMIN(smink,cv_mem->cv_ssdat[i][k]);
      smaxk = SUNMAX(smaxk,cv_mem->cv_ssdat[i][k]);
    }

    if (smink < TINY*smaxk) {
      kflag = -1;
      return(kflag);
    }
    smax[k] = smaxk;
    ssmax[k] = smaxk*smaxk;

    sumrat = ZERO;
    sumrsq = ZERO;
    for (i=1; i<=4; i++) {
      rat[i][k] = cv_mem->cv_ssdat[i][k]/cv_mem->cv_ssdat[i+1][k];
      sumrat = sumrat + rat[i][k];
      sumrsq = sumrsq + rat[i][k]*rat[i][k];
    }
    rav[k] = FOURTH*sumrat;
    vrat[k] = SUNRabs(FOURTH*sumrsq - rav[k]*rav[k]);

    qc[5][k] = cv_mem->cv_ssdat[1][k] * cv_mem->cv_ssdat[3][k] -
               cv_mem->cv_ssdat[2][k] * cv_mem->cv_ssdat[2][k];
    qc[4][k] = cv_mem->cv_ssdat[2][k] * cv_mem->cv_ssdat[3][k] -
               cv_mem->cv_ssdat[1][k] * cv_mem->cv_ssdat[4][k];
    qc[3][k] = ZERO;
    qc[2][k] = cv_mem->cv_ssdat[2][k] * cv_mem->cv_ssdat[5][k] -
               cv_mem->cv_ssdat[3][k] * cv_mem->cv_ssdat[4][k];
    qc[1][k] = cv_mem->cv_ssdat[4][k] * cv_mem->cv_ssdat[4][k] -
               cv_mem->cv_ssdat[3][k] * cv_mem->cv_ssdat[5][k];

    for (i=1; i<=5; i++)
      qco[i][k] = qc[i][k];

  }                            /* End of k loop */

  /* Isolate normal or nearly-normal matrix case. The three quartics will
     have a common or nearly-common root in this case.
     Return a kflag = 1 if this procedure works. If the three roots
     differ more than vrrt2, return error kflag = -3.    */

  vmin = SUNMIN(vrat[1],SUNMIN(vrat[2],vrat[3]));
  vmax = SUNMAX(vrat[1],SUNMAX(vrat[2],vrat[3]));

  if (vmin < vrrtol*vrrtol) {

    if (vmax > vrrt2*vrrt2) {
      kflag = -2;
      return(kflag);
    } else {
      rr = (rav[1] + rav[2] + rav[3])/THREE;
      drrmax = ZERO;
      for (k = 1;k<=3;k++) {
        adrr = SUNRabs(rav[k] - rr);
        drrmax = SUNMAX(drrmax, adrr);
      }
      if (drrmax > vrrt2) {kflag = -3; return(kflag);}

      kflag = 1;

      /*  can compute charactistic root, drop to next section   */
    }

  } else {

    /* use the quartics to get rr. */

    if (SUNRabs(qco[1][1]) < TINY*ssmax[1]) {
      kflag = -4;
      return(kflag);
    }

    tem = qco[1][2]/qco[1][1];
    for (i=2; i<=5; i++) {
      qco[i][2] = qco[i][2] - tem*qco[i][1];
    }

    qco[1][2] = ZERO;
    tem = qco[1][3]/qco[1][1];
    for (i=2; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][1];
    }
    qco[1][3] = ZERO;

    if (SUNRabs(qco[2][2]) < TINY*ssmax[2]) {
      kflag = -4;
      return(kflag);
    }

    tem = qco[2][3]/qco[2][2];
    for (i=3; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][2];
    }

    if (SUNRabs(qco[4][3]) < TINY*ssmax[3]) {
      kflag = -4;
      return(kflag);
    }

    rr = -qco[5][3]/qco[4][3];

    if (rr < TINY || rr > HUNDRED) {
      kflag = -5;
      return(kflag);
    }

    for (k=1; k<=3; k++)
      qkr[k] = qc[5][k] + rr*(qc[4][k] + rr*rr*(qc[2][k] + rr*qc[1][k]));

    sqmax = ZERO;
    for (k=1; k<=3; k++) {
      saqk = SUNRabs(qkr[k])/ssmax[k];
      if (saqk > sqmax) sqmax = saqk;
    }

    if (sqmax < sqtol) {
      kflag = 2;

      /*  can compute charactistic root, drop to "given rr,etc"   */

    } else {

      /* do Newton corrections to improve rr.  */

      for (it=1; it<=3; it++) {
        for (k=1; k<=3; k++) {
          qp = qc[4][k] + rr*rr*(THREE*qc[2][k] + rr*FOUR*qc[1][k]);
          drr[k] = ZERO;
          if (SUNRabs(qp) > TINY*ssmax[k]) drr[k] = -qkr[k]/qp;
          rrc[k] = rr + drr[k];
        }

        for (k=1; k<=3; k++) {
          s = rrc[k];
          sqmaxk = ZERO;
          for (j=1; j<=3; j++) {
            qjk[j][k] = qc[5][j] + s*(qc[4][j] + s*s*(qc[2][j] + s*qc[1][j]));
            saqj = SUNRabs(qjk[j][k])/ssmax[j];
            if (saqj > sqmaxk) sqmaxk = saqj;
          }
          sqmx[k] = sqmaxk;
        }

        sqmin = sqmx[1] + ONE;
        for (k=1; k<=3; k++) {
          if (sqmx[k] < sqmin) {
            kmin = k;
            sqmin = sqmx[k];
          }
        }
        rr = rrc[kmin];

        if (sqmin < sqtol) {
          kflag = 3;
          /*  can compute charactistic root   */
          /*  break out of Newton correction loop and drop to "given rr,etc" */
          break;
        } else {
          for (j=1; j<=3; j++) {
            qkr[j] = qjk[j][kmin];
          }
        }
      } /*  end of Newton correction loop  */

      if (sqmin > sqtol) {
        kflag = -6;
        return(kflag);
      }
    } /*  end of if (sqmax < sqtol) else   */
  } /*  end of if (vmin < vrrtol*vrrtol) else, quartics to get rr. */

  /* given rr, find sigsq[k] and verify rr.  */
  /* All positive kflag drop to this section  */

  for (k=1; k<=3; k++) {
    rsa = cv_mem->cv_ssdat[1][k];
    rsb = cv_mem->cv_ssdat[2][k]*rr;
    rsc = cv_mem->cv_ssdat[3][k]*rr*rr;
    rsd = cv_mem->cv_ssdat[4][k]*rr*rr*rr;
    rd1a = rsa - rsb;
    rd1b = rsb - rsc;
    rd1c = rsc - rsd;
    rd2a = rd1a - rd1b;
    rd2b = rd1b - rd1c;
    rd3a = rd2a - rd2b;

    if (SUNRabs(rd1b) < TINY*smax[k]) {
      kflag = -7;
      return(kflag);
    }

    cest1 = -rd3a/rd1b;
    if (cest1 < TINY || cest1 > FOUR) {
      kflag = -7;
      return(kflag);
    }
    corr1 = (rd2b/cest1)/(rr*rr);
    sigsq[k] = cv_mem->cv_ssdat[3][k] + corr1;
  }

  if (sigsq[2] < TINY) {
    kflag = -8;
    return(kflag);
  }

  ratp = sigsq[3]/sigsq[2];
  ratm = sigsq[1]/sigsq[2];
  qfac1 = FOURTH*(cv_mem->cv_q*cv_mem->cv_q - ONE);
  qfac2 = TWO/(cv_mem->cv_q - ONE);
  bb = ratp*ratm - ONE - qfac1*ratp;
  tem = ONE - qfac2*bb;

  if (SUNRabs(tem) < TINY) {
    kflag = -8;
    return(kflag);
  }

  rrb = ONE/tem;

  if (SUNRabs(rrb - rr) > rrtol) {
    kflag = -9;
    return(kflag);
  }

  /* Check to see if rr is above cutoff rrcut  */
  if (rr > rrcut) {
    if (kflag == 1) kflag = 4;
    if (kflag == 2) kflag = 5;
    if (kflag == 3) kflag = 6;
  }

  /* All positive kflag returned at this point  */

  return(kflag);

}

/*
 * cvEwtSetSV
 *
 * This routine sets ewt as decribed above in the case tol_type = CV_SV.
 * It tests for non-positive components before inverting. cvEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */
int cvEwtSetSV_gpu2(CVodeMem cv_mem, N_Vector cv_ewt, N_Vector weight)
{
  N_VAbs(cv_ewt, cv_mem->cv_tempv);
  N_VLinearSum(cv_mem->cv_reltol, cv_mem->cv_tempv, ONE,
               cv_mem->cv_Vabstol, cv_mem->cv_tempv);
  if (N_VMin(cv_mem->cv_tempv) <= ZERO) return(-1);
  N_VInv(cv_mem->cv_tempv, weight);
  return(0);
}

int cvNlsNewton_gpu2(SolverData *sd, CVodeMem cv_mem, int nflag)
{
  itsolver *bicg = &(sd->bicg);
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval, ier;
  booleantype callSetup;

#ifdef PMC_DEBUG_GPU
  clock_t start;
  start=clock();
#endif

  //double *acor_init = NV_DATA_S(cv_mem->cv_acor_init); //user-supplied value to improve guesses for zn(0)
  double *acor = NV_DATA_S(cv_mem->cv_acor);
  double *cv_y = NV_DATA_S(cv_mem->cv_y);
  double *tempv = NV_DATA_S(cv_mem->cv_tempv);
  double *ftemp = NV_DATA_S(cv_mem->cv_ftemp);

  int znUsedOnNewtonIt=2;//Only used zn[0] and zn[1] //0.01s
  for(int i=0; i<znUsedOnNewtonIt; i++){//cv_qmax+1
    double *zn = NV_DATA_S(cv_mem->cv_zn[i]);
    cudaMemcpy((i*bicg->nrows+bicg->dzn),zn,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  }

  cudaMemcpy(bicg->dacor,acor,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dtempv,tempv,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(bicg->dftemp,ftemp,bicg->nrows*sizeof(double),cudaMemcpyHostToDevice);

#ifdef PMC_DEBUG_GPU
  timeNewtonSendInit+= clock() - start;
  counterSendInit++;
#endif

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
             CV_NO_FAILURES : CV_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (cv_mem->cv_lsetup) {
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
                (cv_mem->cv_nst == 0) ||
                (cv_mem->cv_nst >= cv_mem->cv_nstlp + MSBP) ||
                (SUNRabs(cv_mem->cv_gamrat-ONE) > DGMAX);
  } else {
    cv_mem->cv_crate = ONE;
    callSetup = SUNFALSE;
  }

  /* Call a user-supplied function to improve guesses for zn(0), if one exists */
  //if not, set to zero
  //N_VConst(ZERO, cv_mem->cv_acor_init);
  //cudaMemcpyDToGpu(acor_init, bicg->dacor_init, bicg->nrows);

  if (cv_mem->cv_ghfun) {
    //todo use this ghfun
    /*N_VScale(cv_mem->cv_rl1, cv_mem->cv_zn[1], cv_mem->cv_ftemp);
    cv_mem->cv_ghfun(cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_zn[0],
                     N_VLinearSum(ONE, cv_mem->cv_zn[0], -ONE, cv_mem->cv_last_yn, cv_mem->cv_ftemp);
    retval = cv_mem->cv_ghfun(cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_zn[0],
                              cv_mem->cv_last_yn, cv_mem->cv_ftemp, cv_mem->cv_user_data,
                              cv_mem->cv_tempv, cv_mem->cv_acor_init);
    cv_mem->cv_tempv1, cv_mem->cv_acor_init);
    if (retval<0) return(RHSFUNC_RECVR);*/
  }

  //remove temps, not used in jac
  vtemp1 = cv_mem->cv_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = cv_mem->cv_acor;  /* rename y as vtemp2 for readability     */
  vtemp3 = cv_mem->cv_acor;  /* rename tempv as vtemp3 for readability */

  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call cvNewtonIteration for the Newton iteration itself.      */
  for(;;) {

    /* Load prediction into y vector */
    //N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, cv_mem->cv_acor_init, cv_mem->cv_y);
    //gpu_zaxpby(1.0, bicg->dzn, 1.0, bicg->dacor_init, bicg->dcv_y, bicg->nrows, bicg->blocks, bicg->threads);
    gpu_yequalsx(bicg->dcv_y,bicg->dzn, bicg->nrows, bicg->blocks, bicg->threads);//Consider acor_init=0


    cudaMemcpy(cv_y,bicg->dcv_y,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
    //retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_y,
    //                      cv_mem->cv_ftemp, cv_mem->cv_user_data);
    //int f(realtype t, N_Vector y, N_Vector deriv, void *solver_data)

#ifdef PMC_DEBUG_GPU
    start=clock();
#endif

    retval = f(cv_mem->cv_tn, cv_mem->cv_y, cv_mem->cv_ftemp, cv_mem->cv_user_data);
    if (retval < 0) return(CV_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

#ifdef PMC_DEBUG_GPU
    timeDerivNewton+= clock() - start;
    counterDerivNewton++;
#endif
    //cudaMemcpyDToGpu(ftemp, bicg->dftemp, bicg->nrows);

    cv_mem->cv_nfe++;
    if (callSetup)
    {

#ifdef PMC_DEBUG_GPU
      start=clock();
#endif

      ier = linsolsetup_gpu2(sd, cv_mem, convfail, vtemp1, vtemp2, vtemp3);

#ifdef PMC_DEBUG_GPU
      timeLinSolSetup+= clock() - start;
      counterLinSolSetup++;
#endif

      cv_mem->cv_nsetups++;
      callSetup = SUNFALSE;
      cv_mem->cv_gamrat = cv_mem->cv_crate = ONE;
      cv_mem->cv_gammap = cv_mem->cv_gamma;
      cv_mem->cv_nstlp = cv_mem->cv_nst;
      // Return if lsetup failed
      if (ier < 0) return(CV_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    // Set acor to the initial guess for adjustments to the y vector
    //N_VScale(ONE, cv_mem->cv_acor_init, cv_mem->cv_acor);
    //gpu_yequalsx(bicg->dacor, bicg->dacor_init, bicg->nrows, bicg->blocks, bicg->threads);
    cudaMemset(bicg->dacor, 0.0, bicg->nrows*sizeof(double));

#ifdef PMC_DEBUG_GPU
    start=clock();
#endif

    // Do the Newton iteration
    //ier = cvNewtonIteration(cv_mem);
    ier = linsolsolve_gpu2(sd, cv_mem);

#ifdef PMC_DEBUG_GPU
    timeLinSolSolve+= clock() - start;
    counterLinSolSolve++;
#endif
    // If there is a convergence failure and the Jacobian-related
    //   data appears not to be current, loop again with a call to lsetup
    //   in which convfail=CV_FAIL_BAD_J.  Otherwise return.
    if (ier != TRY_AGAIN) return(ier);

    callSetup = SUNTRUE;
    convfail = CV_FAIL_BAD_J;

  }
}

//todo check for possible changes in jac indices(indexptrs) in sunmatscaleaddI
//todo connect jac and deriv gpu from camp
//int linsolsetup_gpu2(CVodeMem cv_mem)
int linsolsetup_gpu2(SolverData *sd, CVodeMem cv_mem,int convfail,N_Vector vtemp1,N_Vector vtemp2,N_Vector vtemp3)
{
  itsolver *bicg = &(sd->bicg);
  booleantype jbad, jok;
  realtype dgamma;
  CVDlsMem cvdls_mem;
  int retval = 0;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) ||
         (cv_mem->cv_nst > cvdls_mem->nstlj + CVD_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;

  // If jok = SUNTRUE, use saved copy of J
  if (jok) {
    //if (0) {
    cv_mem->cv_jcur = SUNFALSE;
    retval = SUNMatCopy(cvdls_mem->savedJ, cvdls_mem->A);

    /* If jok = SUNFALSE, reset J, call jac routine for new J value and save a copy */
  } else {
    cvdls_mem->nje++;
    cvdls_mem->nstlj = cv_mem->cv_nst;
    cv_mem->cv_jcur = SUNTRUE;
    //retval = SUNMatZero(cvdls_mem->A);//we already set this to zero in our calc_jac

    //clock_t start = clock();

    //int Jac(realtype t, N_Vector y, N_Vector deriv, SUNMatrix J, void *solver_data,
    //        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    //retval = cvdls_mem->jac(cv_mem->cv_tn, cv_mem->cv_y,cv_mem->cv_ftemp, cvdls_mem->A,
    //                        cvdls_mem->J_data, vtemp1, vtemp2, vtemp3);

    retval = Jac(cv_mem->cv_tn, cv_mem->cv_y,cv_mem->cv_ftemp, cvdls_mem->A,cvdls_mem->J_data, vtemp1, vtemp2, vtemp3);

    /*
    timeJacCVODE+= clock() - start5;
    counterJacCVODE++;
    */

    if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVDLS",
                     "cvDlsSetup",  MSGD_JACFUNC_FAILED);
      cvdls_mem->last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      cvdls_mem->last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }

#ifdef PMC_DEBUG_GPU
    clock_t start = clock();
#endif

    retval = SUNMatCopy(cvdls_mem->A, cvdls_mem->savedJ);

#ifdef PMC_DEBUG_GPU
    timeMatCopy+= clock() - start;
    counterMatCopy++;
#endif

    if (retval) {
      cvProcessError(cv_mem, CVDLS_SUNMAT_FAIL, "CVDLS",
                     "cvDlsSetup",  MSGD_MATCOPY_FAILED);
      cvdls_mem->last_flag = CVDLS_SUNMAT_FAIL;
      return(-1);
    }

  }

  clock_t start = clock();

  cudaMemcpy(bicg->dA,bicg->A,bicg->nnz*sizeof(double),cudaMemcpyHostToDevice);

#ifdef PMC_DEBUG_GPU
  timeMatScaleAddISendA+= clock() - start;
  counterMatScaleAddISendA++;

  clock_t start2 = clock();
#endif

  gpu_matScaleAddI(bicg->nrows,bicg->dA,bicg->djA,bicg->diA,-cv_mem->cv_gamma,bicg->blocks,bicg->threads);

#ifdef PMC_DEBUG_GPU
  timeMatScaleAddI+= clock() - start2;
  counterMatScaleAddI++;
#endif

  gpu_diagprecond(bicg->nrows,bicg->dA,bicg->djA,bicg->diA,bicg->ddiag,bicg->blocks,bicg->threads); //Setup linear solver

  //return(cvdls_mem->last_flag);
  return retval;
}

//translating to cv new iteration
int linsolsolve_gpu2(SolverData *sd, CVodeMem cv_mem)
{
  itsolver *bicg = &(sd->bicg);
  int m, retval;
  realtype del, delp, dcon;
#ifdef PMC_DEBUG_GPU
  clock_t start;
#endif

  cv_mem->cv_mnewt = m = 0;

  //Delp = del from last iter (reduce iterations)
  del = delp = 0.0;

  double *acor = NV_DATA_S(cv_mem->cv_acor);
  double *cv_y = NV_DATA_S(cv_mem->cv_y);
  double *tempv = NV_DATA_S(cv_mem->cv_tempv);
  CVDlsMem cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  double *x = NV_DATA_S(cvdls_mem->x);

  // Looping point for Newton iteration
  for(;;) {

    // Evaluate the residual of the nonlinear system
    gpu_zaxpby(cv_mem->cv_rl1, (bicg->dzn+1*bicg->nrows), 1.0, bicg->dacor, bicg->dtempv, bicg->nrows, bicg->blocks, bicg->threads);
    gpu_zaxpby(cv_mem->cv_gamma, bicg->dftemp, -1.0, bicg->dtempv, bicg->dtempv, bicg->nrows, bicg->blocks, bicg->threads);
    //N_VLinearSum(cv_mem->cv_rl1, cv_mem->cv_zn[1], ONE,
    //             cv_mem->cv_acor, cv_mem->cv_tempv);
    //N_VLinearSum(cv_mem->cv_gamma, cv_mem->cv_ftemp, -ONE,
    //             cv_mem->cv_tempv, cv_mem->cv_tempv);

#ifdef PMC_DEBUG_GPU
    start=clock();
#endif

    // Call the lsolve function
    //todo use cuda get event timer to profile a gpu function with multiple calls to gpu without device_syncronhize
    //solveGPU(bicg,bicg->dA,bicg->djA,bicg->diA,bicg->dx,bicg->dtempv);
    solveGPU_multi(bicg,bicg->dA,bicg->djA,bicg->diA,bicg->dx,bicg->dtempv);

    //todo improve checking of good results (under flag debug and check all elements)
    //cudaMemcpy(x,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
    //printf("dx3_4 %f %f,", x[3], x[4]); //seems working

#ifdef PMC_DEBUG_GPU
    timeBiConjGrad+= clock() - start;
    counterBiConjGrad++;
#endif

    // scale the correction to account for change in gamma
    //todo add this part (from matt merge)
    /*if ((cv_mem->cv_lmm == CV_BDF) && (cv_mem->cv_gamrat != 1.0)) {
      //N_VScale(TWO/(ONE + cv_mem->cv_gamrat), b, b);
      gpu_scaley(bicg->dx, 2.0 / (1.0 + cv_mem->cv_gamrat), bicg->nrows, bicg->blocks, bicg->threads);
    }

    if (cv_mem->cv_ghfun) {
      N_VLinearSum(ONE, cv_mem->cv_y, ONE, b, cv_mem->cv_ftemp);
      retval = cv_mem->cv_ghfun(cv_mem->cv_tn, ZERO, cv_mem->cv_ftemp,
                                cv_mem->cv_y, b, cv_mem->cv_user_data,
                                cv_mem->cv_tempv1, cv_mem->cv_tempv2);
      if (retval==1) {
        //SUNDIALS_DEBUG_PRINT_FULL("Received updated adjustment from guess helper");
      } else if (retval<0) {
        if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup))
          return(TRY_AGAIN);
        else
          return(RHSFUNC_RECVR);
      }
    }
    // Check for negative concentrations
    N_VLinearSum(ONE, cv_mem->cv_y, ONE, b, cv_mem->cv_ftemp);
    if (N_VMin(cv_mem->cv_ftemp) < -PMC_TINY) {
      return(CONV_FAIL);
    }
*/


    // Get WRMS norm of correction
    //del = vwrms_gpu2(b, cv_mem->cv_ewt);
    del = gpu_VWRMS_Norm(bicg->nrows, bicg->dx, bicg->dewt, bicg->aux, bicg->daux, (bicg->blocks + 1) / 2, bicg->threads);

    //add correction to acor and y
    // z= a*x + b*y
    gpu_zaxpby(1.0, bicg->dacor, 1.0, bicg->dx, bicg->dacor, bicg->nrows, bicg->blocks, bicg->threads);
    gpu_zaxpby(1.0, bicg->dzn, 1.0, bicg->dacor, bicg->dcv_y, bicg->nrows, bicg->blocks, bicg->threads);

    //(T a, const T *X, T b, const T *Y, T *Z, I n)
    //Z[i] = a*X[i] + b*Y[i];
    //N_VLinearSum(ONE, cv_mem->cv_acor, ONE, b, cv_mem->cv_acor);
    //N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, cv_mem->cv_acor, cv_mem->cv_y);

    // Test for convergence.  If m > 0, an estimate of the convergence
    // rate constant is stored in crate, and used in the test.
    if (m > 0) {
      cv_mem->cv_crate = SUNMAX(0.3 * cv_mem->cv_crate, del / delp);
    }

    dcon = del * SUNMIN(1.0, cv_mem->cv_crate) / cv_mem->cv_tq[4];

    if (dcon <= 1.0) {
      //cv_mem->cv_acnrm = N_VWrmsNorm(cv_mem->cv_acor, cv_mem->cv_ewt);
      cv_mem->cv_acnrm = gpu_VWRMS_Norm(bicg->nrows, bicg->dacor, bicg->dewt, bicg->aux,
                                        bicg->daux, (bicg->blocks + 1) / 2, bicg->threads);
      cv_mem->cv_jcur = SUNFALSE;

      cudaMemcpy(acor,bicg->dacor,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(tempv,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      return (CV_SUCCESS);
    }
    cv_mem->cv_mnewt = ++m;

    // Stop at maxcor iterations or if iter. seems to be diverging.
    //     If still not converged and Jacobian data is not current,
    //     signal to try the solution again
    if ((m == cv_mem->cv_maxcor) || ((m >= 2) && (del > RDIV * delp))) {
      cudaMemcpy(acor,bicg->dacor,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(tempv,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup)) {
        return (TRY_AGAIN);
      } else {
        return (CONV_FAIL);
      }
    }

#ifdef PMC_DEBUG_GPU
    start=clock();
#endif

    // Save norm of correction, evaluate f, and loop again
    delp = del;

    cudaMemcpy(cv_y,bicg->dcv_y,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
    //retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_y,
    //                      cv_mem->cv_ftemp, cv_mem->cv_user_data);
    //int f(realtype t, N_Vector y, N_Vector deriv, void *solver_data)

    retval = f(cv_mem->cv_tn, cv_mem->cv_y, cv_mem->cv_ftemp, cv_mem->cv_user_data);

#ifdef PMC_DEBUG_GPU
    timeDerivSolve+= clock() - start;
    counterDerivSolve++;
#endif
    //cudaMemcpyDToGpu(ftemp, bicg->dftemp, bicg->nrows);

    //N_VLinearSum(ONE, cv_mem->cv_y, -ONE, cv_mem->cv_zn[0], cv_mem->cv_acor);
    // z= a*x + b*y
    gpu_zaxpby(1.0, bicg->dcv_y, -1.0, bicg->dzn, bicg->dacor, bicg->nrows, bicg->blocks, bicg->threads);

    if (retval < 0){
      cudaMemcpy(acor,bicg->dacor,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(tempv,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      return(CV_RHSFUNC_FAIL);
    }
    if (retval > 0) {
      cudaMemcpy(acor,bicg->dacor,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(tempv,bicg->dx,bicg->nrows*sizeof(double),cudaMemcpyDeviceToHost);
      if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup))
        return(TRY_AGAIN);
      else
        return(RHSFUNC_RECVR);
    }

    cv_mem->cv_nfe++;

  }
  //return 0;
}

/*
 * cvHandleFailure
 *
 * This routine prints error messages for all cases of failure by
 * cvHin and cvStep.
 * It returns to CVode the value that CVode is to return to the user.
 */

int cvHandleFailure_gpu2(CVodeMem cv_mem, int flag)
{

  /* Set vector of  absolute weighted local errors */
  /*
  N_VProd(acor, ewt, tempv);
  N_VAbs(tempv, tempv);
  */

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
    case CV_ERR_FAILURE:
      cvProcessError(cv_mem, CV_ERR_FAILURE, "CVODE", "CVode", MSGCV_ERR_FAILS,
                     cv_mem->cv_tn, cv_mem->cv_h);
      break;
    case CV_CONV_FAILURE:
      cvProcessError(cv_mem, CV_CONV_FAILURE, "CVODE", "CVode", MSGCV_CONV_FAILS,
                     cv_mem->cv_tn, cv_mem->cv_h);
      break;
    case CV_LSETUP_FAIL:
      cvProcessError(cv_mem, CV_LSETUP_FAIL, "CVODE", "CVode", MSGCV_SETUP_FAILED,
                     cv_mem->cv_tn);
      break;
    case CV_LSOLVE_FAIL:
      cvProcessError(cv_mem, CV_LSOLVE_FAIL, "CVODE", "CVode", MSGCV_SOLVE_FAILED,
                     cv_mem->cv_tn);
      break;
    case CV_RHSFUNC_FAIL:
      cvProcessError(cv_mem, CV_RHSFUNC_FAIL, "CVODE", "CVode", MSGCV_RHSFUNC_FAILED,
                     cv_mem->cv_tn);
      break;
    case CV_UNREC_RHSFUNC_ERR:
      cvProcessError(cv_mem, CV_UNREC_RHSFUNC_ERR, "CVODE", "CVode", MSGCV_RHSFUNC_UNREC,
                     cv_mem->cv_tn);
      break;
    case CV_REPTD_RHSFUNC_ERR:
      cvProcessError(cv_mem, CV_REPTD_RHSFUNC_ERR, "CVODE", "CVode", MSGCV_RHSFUNC_REPTD,
                     cv_mem->cv_tn);
      break;
    case CV_RTFUNC_FAIL:
      cvProcessError(cv_mem, CV_RTFUNC_FAIL, "CVODE", "CVode", MSGCV_RTFUNC_FAILED,
                     cv_mem->cv_tn);
      break;
    case CV_TOO_CLOSE:
      cvProcessError(cv_mem, CV_TOO_CLOSE, "CVODE", "CVode", MSGCV_TOO_CLOSE);
      break;
    default:
      return(CV_SUCCESS);
  }

  return(flag);
}


void free_gpu2(CVodeMem cv_mem)
{

  //HANDLE_ERROR(cudaFree(cv_mem->indexvals_gpu2));
  //HANDLE_ERROR(cudaFree(cv_mem->indexptrs_gpu2));
  //HANDLE_ERROR(cudaFree(cv_mem->jac_data_gpu2));

  //In principle, C++ guarantee destroy the classes when they go out of scope, so don't need to call destructor
  //bicg->~itsolver(){};

}

void printSolverCounters()
{
#ifdef PMC_DEBUG_GPU
  //printf("timeNewtonSendInit %lf, counterSendInit %d\n",timeNewtonSendInit/CLOCKS_PER_SEC,counterSendInit);
  //printf("timeMatScaleAddI %lf, counterMatScaleAddI %d\n",timeMatScaleAddI/CLOCKS_PER_SEC,counterMatScaleAddI);
  //printf("timeMatScaleAddISendA %lf, counterMatScaleAddISendA %d\n",timeMatScaleAddISendA/CLOCKS_PER_SEC,counterMatScaleAddISendA);
  //printf("timeMatCopy %lf, counterMatCopy %d\n",timeMatCopy/CLOCKS_PER_SEC,counterMatCopy);
  printf("timeDerivNewton %lf, counterDerivNewton %d\n",timeDerivNewton/CLOCKS_PER_SEC,counterDerivNewton);
  printf("timeDerivSolve %lf, counterDerivSolve %d\n",timeDerivSolve/CLOCKS_PER_SEC,counterDerivSolve);
  printf("timeBiConjGrad %lf, counterBiConjGrad %d\n",timeBiConjGrad/CLOCKS_PER_SEC,counterBiConjGrad);
  printf("timeLinSolSetup %lf, counterLinSolSetup %d\n",timeLinSolSetup/CLOCKS_PER_SEC,counterLinSolSetup);
  printf("timeLinSolSolve %lf, counterLinSolSolve %d\n",timeLinSolSolve/CLOCKS_PER_SEC,counterLinSolSolve);
  printf("timeNewtonIt %lf, counterNewtonIt %d\n",timeNewtonIt/CLOCKS_PER_SEC,counterNewtonIt);
  printf("timecvStep %lf, countercvStep %d\n",timecvStep/CLOCKS_PER_SEC,countercvStep);
  printf("timeprecvStep %lf, counterBiConjGrad %d\n",timeBiConjGrad/CLOCKS_PER_SEC,counterBiConjGrad);

  //printf("timeJacCVODE %lf, counterJacCVODE %d\n",timeJacCVODE/CLOCKS_PER_SEC,counterJacCVODE);
  //printf("timeprecvStep %lf, counterBiConjGrad %d\n",timeBiConjGrad/CLOCKS_PER_SEC,counterBiConjGrad);
#endif
}


/*
//cpu - for testing
void sun(int nrows, realtype* dA, sunindextype* djA, sunindextype* diA, double gamma) {

  for (int j=0; j < nrows; j++){
    for (int i = diA[j]; i < diA[j + 1]; i++)
      if (djA[i] == j) {
        dA[i] = 1.0 + gamma * dA[i];
      } else {
        dA[i] = gamma * dA[i];
      }
  }
}

__global__ void matscaleaddi_global(int nrows, double* dA, sunindextype* djA, sunindextype* diA, double gamma) {

  int j = threadIdx.x + blockDim.x*blockIdx.x;
  for (int i = diA[j]; i < diA[j + 1]; i++)
    if (djA[i] == j) {
      dA[i] = 1.0 + gamma * dA[i];
    } else {
      dA[i] = gamma * dA[i];
    }
}

 int matscaleaddi_gpu2(double gamma, SUNMatrix J, CVodeMem cv_mem)
{


  //Note that don't need to pass again J since pointer A it's the same than J and
  //its updated automatically!

  //bicg->matScaleAddI(gamma);

  cudaMemcpyDToGpu(bicg->A,bicg->dA,bicg->nnz);
  gpu_matScaleAddI(bicg->nrows,bicg->dA,bicg->djA,bicg->diA,gamma,bicg->blocks,bicg->threads);

  return 0;
}

*/