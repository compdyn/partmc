/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * This is the c ODE solver for the chemistry module 
 * It is currently set up to use the SUNDIALS BDF method, Newton
 * iteration with the CVDENSE sparse linear solver, and a user-supplied
 * Jacobian routine. 
 *
 * It uses a scalar relative tolerance and a vector absolute tolerance.
 *
*/
/** \file
 * \brief Interface to c solvers for chemistry
*/

#include <stdio.h>
#include <stdlib.h>

/* Header files with a description of contents used */

#if defined(PMC_USE_SUNDIALS)
#include <cvodes/cvodes.h>              /* prototypes for CVODE fcts., consts.  */
#include <cvodes/cvodes_spils.h>        /* CVSpils interface                    */
#include <nvector/nvector_serial.h>     /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_sparse.h> /* sparse SUNMatrix                     */
#include <sundials/sundials_dense.h>    /* SUNDIALS dense matrix                */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_klu.h>    /* KLU SUNLinearSolver                  */
#include <cvodes/cvodes_direct.h>       /* CVDls interface                      */
#include <sundials/sundials_types.h>    /* definition of type realtype          */
#include <sundials/sundials_math.h>     /* SUNDIALS math function macros        */
#endif

/// Flag to disallow negative concentrations
#define PMC_INTEGRATION_CHECK_NEGATIVE 0
/// Maximum number of convergence failures
#define PMC_INTEGRATION_MAX_CONV_FAILS 20

/// Result code indicating successful completion.
#define PMC_INTEGRATION_SUCCESS                0
/// Result code indicating no available integration routine
#define PMC_INTEGRATION_NO_AVAIL_SOLVER        12

// SUNDIALS result codes
#if defined(PMC_USE_SUNDIALS)
/// Result code indicating failure to allocate \c y vector.
#define PMC_INTEGRATION_INIT_Y                 1
/// Result code indicating failure to allocate \c abstol vector.
#define PMC_INTEGRATION_INIT_ABSTOL            2
/// Result code indicating failure to create the solver.
#define PMC_INTEGRATION_INIT_CVODE_MEM         3
/// Result code indicating failure to initialize the solver.
#define PMC_INTEGRATION_INIT_CVODE             4
/// Result code indicating failure to set tolerances.
#define PMC_INTEGRATION_SVTOL                  5
/// Result code indicating failure to set maximum steps.
#define PMC_INTEGRATION_SET_MAX_STEPS          6
/// Result code indicating failure of the solver.
#define PMC_INTEGRATION_FAIL                   7
/// Result code indicating failure to initialize sparse Jacobian
#define PMC_INTEGRATION_SPARSE_JAC             8
/// Result code indicating failure to set Jacobian function
#define PMC_INTEGRATION_JAC_FUNC               9
/// Result code indicating failure to set user data
#define PMC_INTEGRATION_SET_USER_DATA          10
/// Result code indicating SUNDIALS realtype is not set to double precision
#define PMC_INTEGRATION_WRONG_PRECISION        11
/// Result code indicating failure to get the KLU linear solver
#define PMC_INTEGRATION_KLU_LINEAR_SOLVER      13
/// Result code indicating failure to set the linear solver
#define PMC_INTEGRATION_SET_LINEAR_SOLVER      14
/// Result code indicating failure to set the maximum number of convergence failures
#define PMC_INTEGRATION_SET_MAX_CONV_FAILS     15
/// Result code indicating failure to set the SPGMR solver
#define PMC_INTEGRATION_SPGMR_LINEAR_SOLVER    16
/// Result code indicating failure to set the preconditioner functions
#define PMC_INTEGRATION_SET_PRECONDITIONER     17
#endif

/* Type : UserData
 * contains preconditioner blocks, pivot arrays, and void pointer to send to fortran 
 * functions.
 */
typedef struct {
  DlsMat P, Jbd;
  sunindextype *pivot;
  void *sysdata;
  SUNMatrix jac_init;
} *UserData;

/* Fortran support subroutines */
#ifndef DOXYGEN_SKIP_DOC
// Fortran subroutine to calcualte f(t,y) for the system
void deriv_func(int n_eqn_c, double curr_time_c, double *state_c_p,
                    double *f_c_p, void *sysdata_c_p);

// Fortran subroutine to calculate J(t,y) = df/dy for the system
void jac_func(int n_eqn_c, double curr_time_c, double *state_c_p,
                      double *jac_c_p, void *sysdata_c_p);
#endif

/* Functions Called by the Solver */

/// \brief f routine. Compute function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/// \brief Jacobian routine. Compute J(t,y) = df/dy.
static int Jac(realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/// Preconditioner
static int Precond(realtype tn, N_Vector c, N_Vector fc, booleantype jok,
		   booleantype *jcurPtr, realtype gamma, void *user_data);

/// Solver for the preconditioner
static int PSolve(realtype tn, N_Vector c, N_Vector fc, N_Vector r, N_Vector z,
		  realtype gamma, realtype delta, int lr, void *user_data);

#ifndef DOXYGEN_SKIP_DOC
static int test_f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int test_Jac(realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

/// \brief Private function to print species concentrations
static void integrate_print_output(realtype t, realtype y1, realtype y2, realtype y3);
/// \brief Private function to print final statistics
static void integrate_print_stats(void *cvode_mem, int n_eqn, N_Vector y_nv);

#if defined(PMC_USE_SUNDIALS)
/// \brief Private function to check SUNDIALS return values
static int integrate_SUNDIALS_check_flag(void *flagvalue, char *funcname, int opt);
#endif

/// \brief Main ODE solver
/**
 *----------------------------
 * Main ODE solver function
 *----------------------------
 * \param neq The number of equations
 * \param y A pointer to the aq. chemistry state array (atm or M)
 * \param abstol A pointer to array of absolute tolerances
 *             for each species (same units as in x)
 * \param reltol The relative tolerance for all species (unitless)
 * \param max_steps The maximum number of integration steps (unitless)
 * \param t_initial The start time (s)
 * \param t_final The end time (s)
 * \param sysdata A pointer to system data for use in f(t,y) and J(t,y) functions
 * \param n_jac_elem The exact number of non-zero elements in the sparse Jacobian matrix
 * \param n_jac_col_elem An array containing the number of non-zero elements in each column
 *          of the Jacobian matrix. Their sum should equal n_jac_elem.
 * \param jac_row_ids An array containing the row id of each element in the Jacobian matrix
 *          data arranged as a flat array. The array should contain n_jac_elem elements.
 */
int integrate_solver (int neq, double *y, double *abstol, double reltol,
                  int max_steps, double t_initial, double t_final, void *sysdata,
		  int n_jac_elem, int *jac_col_ptrs, int *jac_row_ids)
{
    
#if defined(PMC_USE_SUNDIALS)
#if defined(SUNDIALS_DOUBLE_PRECISION)
    /* SUNDIALS VARIABLES */
    realtype reltol_rt; // relative tolerance
    realtype t_rt;      // current time
    N_Vector y_nv;      // state variable (atm or M)
    N_Vector abstol_nv; // absolulte tolerances (atm or M)
    SUNMatrix A;
    SUNLinearSolver LS;
    UserData user_data; // sysdata and preconditioner data

    void *cvode_mem;    // CVODE solver parameters
    int flag;

    // Set up user data
    user_data = (UserData) malloc(sizeof *user_data);
#ifdef KRYLOV
    user_data->P = NewDenseMat(neq, neq);
    user_data->Jbd = NewDenseMat(neq, neq);
    user_data->pivot = newIndexArray(neq);
#endif
    user_data->sysdata = sysdata;
    user_data->jac_init = NULL;

    y_nv = abstol_nv = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;
    
    // set up N_Vectors and point them to input vectors
    // !!! SUNDIALS realtype must be set to double precision floating point at install !!!
    y_nv = N_VMake_Serial(neq, (realtype*)y);
    if (integrate_SUNDIALS_check_flag((void *)y, "N_VNew_Serial", 0)) return(PMC_INTEGRATION_INIT_Y);
    abstol_nv = N_VMake_Serial(neq, (realtype*)abstol);
    if (integrate_SUNDIALS_check_flag((void *)abstol, "N_VNew_Serial", 0)) return(PMC_INTEGRATION_INIT_ABSTOL);
    
    reltol_rt = (realtype) reltol;
    
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (integrate_SUNDIALS_check_flag((void *)cvode_mem, "CVodeCreate", 0))
        return(PMC_INTEGRATION_INIT_CVODE_MEM);
    
    /* Set user data */
    flag = CVodeSetUserData(cvode_mem, user_data);
    if (integrate_SUNDIALS_check_flag(&flag, "CVodeSetUserData", 1))
        return(PMC_INTEGRATION_SET_USER_DATA);
    
    /* Call CVodeInit to initialize the integrator memory and specify the
     * right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, f, (realtype) t_initial, y_nv);
    if (integrate_SUNDIALS_check_flag(&flag, "CVodeInit", 1))
        return(PMC_INTEGRATION_INIT_CVODE);
    
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol_rt, abstol_nv);
    if (integrate_SUNDIALS_check_flag(&flag, "CVodeSVtolerances", 1))
        return(PMC_INTEGRATION_SVTOL);

    /* Call CVodeSetMaxNumSteps to specify the maximum nunmber of
     * iterations */
    flag = CVodeSetMaxNumSteps(cvode_mem, max_steps);
    if (integrate_SUNDIALS_check_flag(&flag, "CVodeSetMaxNumSteps", 1))
        return PMC_INTEGRATION_SET_MAX_STEPS;
    
    /* Call CVodeSetMaxConvFails to specify the maximum number of
     * convergence failures */
    flag = CVodeSetMaxConvFails(cvode_mem, PMC_INTEGRATION_MAX_CONV_FAILS);
    if (integrate_SUNDIALS_check_flag(&flag, "CVodeSetMaxConvFails",1))
        return PMC_INTEGRATION_SET_MAX_CONV_FAILS;

#ifdef KRYLOV
    /* Create an SPGMR linear solver */
    LS = SUNSPGMR(y_nv, PREC_LEFT, 0);
    if (integrate_SUNDIALS_check_flag((void *)LS, "SUNSPGMR", 0))
	return(PMC_INTEGRATION_SPGMR_LINEAR_SOLVER);

    /* Set the SPGMR linear solver */
    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    if (integrate_SUNDIALS_check_flag(&flag, "CVSpilsSetLinearSolver", 1))
	return(PMC_INTEGRATION_SET_LINEAR_SOLVER);

    /* Set the preconditioner setup and solve functions */
    flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
    if (integrate_SUNDIALS_check_flag(&flag, "CVSpilsSetPreconditioner", 1))
	return(PMC_INTEGRATION_SET_PRECONDITIONER);

    /* Use stability detection */
    flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
#else
    /* Create sparse SUNMatrix for use in linear solves */
    A = SUNSparseMatrix(neq, neq, n_jac_elem, CSC_MAT);
    if (integrate_SUNDIALS_check_flag((void *)A, "SUNMatrix", 0))
        return(PMC_INTEGRATION_SPARSE_JAC);
    for (int i=0; i<=neq; i++) (SM_INDEXPTRS_S(A))[i] = (sunindextype) jac_col_ptrs[i];
    for (int i=0; i<n_jac_elem; i++) (SM_INDEXVALS_S(A))[i] = (sunindextype) jac_row_ids[i];
    for (int i=0; i<n_jac_elem; i++) (SM_DATA_S(A))[i] = (realtype) 1.0;

    /* Create a clone of the Jacobian to use during calls to the Jacobian
     * function to preserve the original sparse matrix dimensions */
    user_data->jac_init = SUNMatClone(A);
    SUNMatCopy(A, user_data->jac_init);

    /* Create KLU SUNLinearSolver object for use by CVode */
    LS = SUNKLU(y_nv, A);
    if (integrate_SUNDIALS_check_flag((void *)LS, "SUNKLULinearSolver", 0))
	return(PMC_INTEGRATION_KLU_LINEAR_SOLVER);

    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (integrate_SUNDIALS_check_flag(&flag, "CVDlsSetLinearSolver", 1))
	return(PMC_INTEGRATION_SET_LINEAR_SOLVER);

    /* Set the Jacobian routine to Jac */
    flag = CVDlsSetJacFn(cvode_mem, Jac);
    if (integrate_SUNDIALS_check_flag(&flag, "CVDlsSetJacFn", 1))
        return(PMC_INTEGRATION_JAC_FUNC);
#endif

    /* Run the solver */
    flag = CVode(cvode_mem, (realtype) t_final, y_nv, &t_rt, CV_NORMAL);
    if (integrate_SUNDIALS_check_flag(&flag, "CVode", 1))
    {
        // DIAGNOSTIC
        integrate_print_stats(cvode_mem, neq, y_nv);

        return(PMC_INTEGRATION_FAIL);
    }

    /* Free y and abstol vectors */
    /* data pointer is retained, so y and abstol remain allocated */
    N_VDestroy(y_nv);
    N_VDestroy(abstol_nv);
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */
    SUNLinSolFree(LS);

    /* Free the matrix memory */
    SUNMatDestroy(A);
    SUNMatDestroy(user_data->jac_init);

    /* Free the userdata memory */
#ifdef KRYLOV
    DestroyMat(user_data->P);
    DestroyMat(user_data->Jbd);
    destroyArray(user_data->pivot);
#endif
    free(user_data);


    return PMC_INTEGRATION_SUCCESS;
    
#else

    // SUNDIALS is available but the realtype is not set to double precision
    return PMC_INTEGRATION_WRONG_PRECISION;
    
#endif
#endif
    
    // Currently there is no default non-SUNDIALS integration function,
    // so return an error if SUNDIALS is not available
    return PMC_INTEGRATION_NO_AVAIL_SOLVER;
}

#if defined(PMC_USE_SUNDIALS)
/// \brief Check the return value from SUNDIALS call
/**
 * --------------------------------------------
 * Check the return value from a SUNDIALS call.
 *---------------------------------------------
 *  - opt == 0 means SUNDIALS function allocates memory
 *            so check if flagvalue is not a NULL pointer.

 *  - opt == 1 means SUNDIALS function returns a flag so
 *            check if the int pointed to by flagvalue has
 *            flag >= 0.

 *  - opt == 2 means function allocates memory so check
 *            if flagvalue is not a NULL pointer.
 *
 * flagvalue - A pointer to check (either for NULL, or as an
 *             int pointer giving the flag value).
 * funcname  - A string giving the function name returning this
 *             result code.
 * opt       - A flag indicating the type of check to perform.
 * return    - A result code (0 is success).
 */
static int integrate_SUNDIALS_check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
#endif

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/**
 *-----------
 * f routine
 *-----------
 * \param t The current time
 * \param y A pointer to the aq. chemistry state array (atm or M)
 * \param ydot f(t,y) to be calculated
 * \param *user_data A pointer to system data for use in f(t,y) and J(t,y) functions
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    int n_eqn;      // Number of equations
    double *state;  // Pointer to state vector
    double *f;      // Pointer to f(t,y)
    UserData data;  // Integration and model data

    int i;
    
    n_eqn = (int) NV_LENGTH_S(y);
    state = (double*) NV_DATA_S(y);
    f     = (double*) NV_DATA_S(ydot);
    data  = (UserData) user_data;
   
    // Check for negatives
    if (PMC_INTEGRATION_CHECK_NEGATIVE) {
        for(int i=0; i<n_eqn; i++) if(state[i]<0.0) return (1);
    }

    // Call fortran f(t,y) function
    deriv_func(n_eqn, (double) t, state, f, data->sysdata);
    
    return(0);
}

/**
 *------------------
 * Jacobian routine
 *------------------
 * \param t The current time
 * \param y A pointer to the state array 
 * \param fy Calculated f(t,y)
 * \param J Jacobian matrix J(t,y) to be calculated
 * \param *user_data A pointer to system data for use in f(t,y) and J(t,y) functions
 * \param tmp1, tmp2, tmp3 Unused
 */
static int Jac(realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    double *state;  // Pointer to state vector
    double *jac;    // Pointer to J(t,y)
    int N;          // Size of J(t,y)
    UserData data;  // Integration and model data

    int i, j;
    
    state = (double*) NV_DATA_S(y);
    N     = NV_LENGTH_S(y);
    data  = (UserData) user_data;
 
    // Reset the Jacobian dimensions
    if (SM_NNZ_S(J)<SM_NNZ_S(data->jac_init)) {
	SM_INDEXVALS_S(J) = realloc(SM_INDEXVALS_S(J), SM_NNZ_S(data->jac_init)*sizeof(int));
	SM_DATA_S(J) = realloc(SM_DATA_S(J), SM_NNZ_S(data->jac_init)*sizeof(realtype));
    }
    SM_NNZ_S(J) = SM_NNZ_S(data->jac_init);
    for (i=0; i<SM_NNZ_S(J); i++) {
	(SM_DATA_S(J))[i] = (realtype)0.0;
	(SM_INDEXVALS_S(J))[i] = (SM_INDEXVALS_S(data->jac_init))[i];
    }
    for (i=0; i<=SM_NP_S(J); i++) {
	(SM_INDEXPTRS_S(J))[i] = (SM_INDEXPTRS_S(data->jac_init))[i];
    }
    
    // Get a pointer to the data
    jac   = (double*) (SM_DATA_S(J));
    
    // Call fortran jacobian function
    // (jac is square so we can use N for neq)
    jac_func(N, (double) t, state, jac, data->sysdata);
    
    return(0);
}

/**
 *---------------------
 * Preconditioner setup
 *---------------------
 * \param t The current time
 * \param y A pointer to the state array
 * \param fy Calculated f(t,y)
 * \param jok A flag that when SUNTRUE indicates Jacobian data can be reused
 * \param jcurPtr Set to SUNTRUE if Jacobian data was reused, SUNFALSE otherwise
 * \param gamma Scalar gamma in Newton matrix M = I - gamma * J
 * \param user_data A pointer to system data
 */
static int Precond(realtype t, N_Vector y, N_Vector fy, booleantype jok,
		   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  UserData data;       // Preconditioner data and sysdata
  DlsMat P, Jbd;       // Preconditioner and Jacobian matrices
  sunindextype *pivot; // pivot data

  double *state;  // Pointer to state vector
  int n_eqn;      // Number of species
  
  data  = (UserData) user_data;
  P     = data->P;
  Jbd   = data->Jbd;
  pivot = data->pivot;
  n_eqn = (int) NV_LENGTH_S(y);

  if (jok) {
    // Copy Jbd to P
    DenseCopy(Jbd, P);
    *jcurPtr = SUNFALSE;
  } else {
    // Generate new Jacobian data
    state = (double*) NV_DATA_S(y);
    // FIXME Need to adjust to use SUNSparseMatrix Jacobian
    jac_func(n_eqn, (double) t, state, (double*) Jbd->data, data->sysdata);
    DenseCopy(Jbd, P);
    *jcurPtr = SUNTRUE;
  }

  // Scale by -gamma
  DenseScale(-gamma, P);

  // Add identity matrix and do LU decomposition
  AddIdentity(P);
  if (DenseGETRF(P, pivot) != 0) return(1); 

  return(0);

}

/**
 *----------------------
 * Preconditioner solver
 *----------------------
 * \param t The current time
 * \param y A pointer to the state array
 * \param fy Calculated f(t,y)
 * \param r The right-hand side vector of the linear system
 * \param z The vector to compute
 * \param gamma The scaling factor in the Newton matrix M = I - gamma * J
 * \param delta Input tolerance for iterative methods
 * \param lr Flag for left (1) or right (2) preconditioner
 * \param user_data Pointer to preconditioner data and sysdata
 */
static int PSolve(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
		  realtype gamma, realtype delta, int lr, void *user_data)
{
  UserData data;
  DlsMat P;
  sunindextype *pivot;
  int n_eqn;
  realtype *zdata;

  data = (UserData) user_data;
  P = data->P;
  pivot = data->pivot;
  zdata = N_VGetArrayPointer(z);

  N_VScale(1.0, r, z);

  n_eqn = (int) NV_LENGTH_S(y);

  // Solve the system Px = r using LU factors stored in P and
  // pivot data in pivot. Return the solution in z
  DenseGETRS(P, pivot, zdata);

  return(0);
}

#ifndef DOXYGEN_SKIP_DOC
/*
 *------------------------------------
 * DIAGNOSTIC
 * Example f(t,y) and J(t,y) functions
 * and statistic printer function
 * from SUNDIALS cvRoberts_dns.c
 *------------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

static int test_f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, yd1, yd3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
  yd3 = Ith(ydot,3) = RCONST(3.0e7)*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;

  return(0);
}


/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int test_Jac(realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = RCONST(-0.04);
  IJth(J,1,2) = RCONST(1.0e4)*y3;
  IJth(J,1,3) = RCONST(1.0e4)*y2;
  IJth(J,2,1) = RCONST(0.04); 
  IJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
  IJth(J,2,3) = RCONST(-1.0e4)*y2;
  IJth(J,3,2) = RCONST(6.0e7)*y2;

  return(0);
}
#endif

/*
 * Print species concentrations
 */

static void integrate_print_output(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4le      y =%14.6le  %14.6le  %14.6le\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

/*
 * Get and print some final statistics
 */

static void integrate_print_stats(void *cvode_mem, int n_eqn, N_Vector y_nv)
{
    long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
    int flag, i;
    
    realtype sug_tol;
    N_Vector err_weight, est_error;
    
    err_weight = N_VNew_Serial((long int) n_eqn);
    est_error = N_VNew_Serial((long int) n_eqn);
    
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumSteps", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumRhsEvals", 1);
    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumErrTestFails", 1);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
    
    flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
    integrate_SUNDIALS_check_flag(&flag, "CVDlsGetNumJacEvals", 1);
    flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
    integrate_SUNDIALS_check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
    
    flag = CVodeGetNumGEvals(cvode_mem, &nge);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumGEvals", 1);

    flag = CVodeGetTolScaleFactor(cvode_mem, &sug_tol);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetTolScaleFactor", 1);
    flag = CVodeGetErrWeights(cvode_mem, err_weight);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetErrWeights", 1);
    flag = CVodeGetEstLocalErrors(cvode_mem, est_error);
    integrate_SUNDIALS_check_flag(&flag, "CVodeGetEstLocalErrors", 1);

    printf("\nFinal Statistics:\n");
    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
           nst, nfe, nsetups, nfeLS, nje);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
           nni, ncfn, netf, nge);
    printf("suggested tolerance scaling factor = %lg\n\n", sug_tol);
    
    printf("SpecID\tErr Wt\tEst Error\tProduct\tValue\n");
    for (i=0; i<n_eqn; i++) {
        printf("%d\t%lg\t%lg\t%lg\t%lg\n", i+1, Ith(err_weight,i+1), Ith(est_error,i+1),
               Ith(err_weight,i+1)*Ith(est_error,i+1), Ith(y_nv, i+1));
    }
    
}




