/* Copyright (C) 2015 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * This is the c ODE solver for the aqueous chemistry mechanism
 * It is currently set up to use the SUNDIALS BDF method, Newton
 * iteration with the CVDENSE dense linear solver, and a user-supplied
 * Jacobian routine. 
 *
 * It uses a scalar relative tolerance and a vector absolute tolerance.
 *
 * Other integration methods may be used by modifying the 
 * aq_integrate_solver function.
 *
*/
/** \file
 * \brief Interface to c solvers for aqueous-phase chemistry
*/
#include <stdio.h>

/* Header files with a description of contents used */

#if defined(PMC_USE_SUNDIALS)
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#endif

/// Result code indicating successful completion.
#define PMC_AQ_INTEGRATE_SUCCESS        0
/// Result code indicating no available integration routine
#define PMC_AQ_INTEGRATE_NO_AVAIL_SOLVER   12

// SUNDIALS result codes
#if defined(PMC_USE_SUNDIALS)
/// Result code indicating failure to allocate \c y vector.
#define PMC_AQ_INTEGRATE_INIT_Y         1
/// Result code indicating failure to allocate \c abstol vector.
#define PMC_AQ_INTEGRATE_INIT_ABSTOL    2
/// Result code indicating failure to create the solver.
#define PMC_AQ_INTEGRATE_INIT_CVODE_MEM 3
/// Result code indicating failure to initialize the solver.
#define PMC_AQ_INTEGRATE_INIT_CVODE     4
/// Result code indicating failure to set tolerances.
#define PMC_AQ_INTEGRATE_SVTOL          5
/// Result code indicating failure to set maximum steps.
#define PMC_AQ_INTEGRATE_SET_MAX_STEPS  6
/// Result code indicating failure of the solver.
#define PMC_AQ_INTEGRATE_FAIL           7
/// Result code indicating failure to set dense Jacobian solver
#define PMC_AQ_INTEGRATE_DENSE_JAC      8
/// Result code indicating failure to set Jacobian function
#define PMC_AQ_INTEGRATE_JAC_FUNC       9
/// Result code indicating failure to set user data
#define PMC_AQ_INTEGRATE_SET_USER_DATA  10
/// Result code indicating SUNDIALS realtype is not set to double precision
#define PMC_AQ_INTEGRATE_WRONG_PRECISION  11
#endif

// Macros for accessing N_Vector and DlsMat elements with Fortran-type indices
/// Ith numbers components 1..NEQ
#define Ith(v,i)    NV_Ith_S(v,i-1)
/// IJth numbers rows,cols 1..NEQ
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Fortran support subroutines */
#ifndef DOXYGEN_SKIP_DOC
// Fortran subroutine to calcualte f(t,y) for the aqueous-phase chemical mechanism
void aq_integrate_f(int n_eqn_c, double curr_time_c, double *state_c_p,
                    double *f_c_p, void *sysdata_c_p);

// Fortran subroutine to calculate J(t,y) = df/dy for the aqueous-phase chemical mechanism
void aq_integrate_jac(int n_eqn_c, double curr_time_c, double *state_c_p,
                      double *jac_c_p, void *sysdata_c_p);
#endif

/* Functions Called by the Solver */

/// \brief f routine. Compute function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/// \brief Jacobian routine. Compute J(t,y) = df/dy.
static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#ifndef DOXYGEN_SKIP_DOC
static int test_f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int test_Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
/// \brief Private function to print species concentrations
static void aq_integrate_print_output(realtype t, realtype y1, realtype y2, realtype y3);
/// \brief Private function to print final statistics
static void aq_integrate_print_stats(void *cvode_mem, int n_eqn);

#if defined(PMC_USE_SUNDIALS)
/// \brief Private function to check SUNDIALS return values
static int aq_integrate_SUNDIALS_check_flag(void *flagvalue, char *funcname, int opt);
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
 */
int aq_integrate_solver (int neq, double *y, double *abstol, double reltol,
                  int max_steps, double t_initial, double t_final, void *sysdata)
{
    
#if defined(PMC_USE_SUNDIALS)
#if defined(SUNDIALS_DOUBLE_PRECISION)
    /* SUNDIALS VARIABLES */
    realtype reltol_rt; // relative tolerance
    realtype t_rt;      // current time
    N_Vector y_nv;      // state variable (atm or M)
    N_Vector abstol_nv; // absolulte tolerances (atm or M)
    
    void *cvode_mem;    // CVODE solver parameters
    int flag;
    
    // set up N_Vectors and point them to input vectors
    // !!! SUNDIALS realtype must be set to double precision floating point at install !!!
    y_nv = N_VMake_Serial(neq, (realtype*)y);
    if (aq_integrate_SUNDIALS_check_flag((void *)y, "N_VNew_Serial", 0)) return(PMC_AQ_INTEGRATE_INIT_Y);
    abstol_nv = N_VMake_Serial(neq, (realtype*)abstol);
    if (aq_integrate_SUNDIALS_check_flag((void *)abstol, "N_VNew_Serial", 0)) return(PMC_AQ_INTEGRATE_INIT_ABSTOL);
    
    reltol_rt = (realtype) reltol;
    
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (aq_integrate_SUNDIALS_check_flag((void *)cvode_mem, "CVodeCreate", 0))
        return(PMC_AQ_INTEGRATE_INIT_CVODE_MEM);
    
    /* Call CVodeInit to initialize the integrator memory and specify the
     * right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, f, (realtype) t_initial, y_nv);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVodeInit", 1))
        return(PMC_AQ_INTEGRATE_INIT_CVODE);
    
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol_rt, abstol_nv);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVodeSVtolerances", 1))
        return(PMC_AQ_INTEGRATE_SVTOL);

    /* Call CVodeSetMaxNumSteps to specify the maximum nunmber of
     * iterations */
    flag = CVodeSetMaxNumSteps(cvode_mem, max_steps);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVodeSetMaxNumSteps", 1))
        return PMC_AQ_INTEGRATE_SET_MAX_STEPS;
    
    /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, neq);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVDense", 1))
        return(PMC_AQ_INTEGRATE_DENSE_JAC);
    
    /* Set the Jacobian routine to Jac */
    flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVDlsSetDenseJacFn", 1))
        return(PMC_AQ_INTEGRATE_JAC_FUNC);

    /* Set user data */
    flag = CVodeSetUserData(cvode_mem, sysdata);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVodeSetUserData", 1))
        return(PMC_AQ_INTEGRATE_SET_USER_DATA);
    
    /* Run the solver */
    flag = CVode(cvode_mem, (realtype) t_final, y_nv, &t_rt, CV_NORMAL);
    if (aq_integrate_SUNDIALS_check_flag(&flag, "CVode", 1))
    {
        // DIAGNOSTIC
        aq_integrate_print_stats(cvode_mem, neq);

        return(PMC_AQ_INTEGRATE_FAIL);
    }

    /* Free y and abstol vectors */
    /* data pointer is retained, so y and abstol remain allocated */
    N_VDestroy_Serial(y_nv);
    N_VDestroy_Serial(abstol_nv);
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return PMC_AQ_INTEGRATE_SUCCESS;
    
#else

    // SUNDIALS is available but the realtype is not set to double precision
 
	return PMC_AQ_INTEGRATE_WRONG_PRECISION;
    
#endif
#endif
    
    // Currently there is no default non-SUNDIALS integration function,
    // so return an error if SUNDIALS is not available
 
	return PMC_AQ_INTEGRATE_NO_AVAIL_SOLVER;
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
static int aq_integrate_SUNDIALS_check_flag(void *flagvalue, char *funcname, int opt)
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
    
    int i;
    
    n_eqn = (int) NV_LENGTH_S(y);
    state = (double*) NV_DATA_S(y);
    f     = (double*) NV_DATA_S(ydot);
    
    // Call fortran f(t,y) function
    aq_integrate_f(n_eqn, (double) t, state, f, user_data);
    
    return(0);
}

/**
 *------------------
 * Jacobian routine
 *------------------
 * \param N The width of the J matrix 
 * \param t The current time
 * \param y A pointer to the aq. chemistry state array (atm or M)
 * \param fy Calculated f(t,y)
 * \param J Jacobian matrix J(t,y) to be calculated
 * \param *user_data A pointer to system data for use in f(t,y) and J(t,y) functions
 * \param tmp1, tmp2, tmp3 Unused
 */
static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    double *state;  // Pointer to state vector
    double *jac;    // Pointer to J(t,y)
    
    int i, j;
    
    state = (double*) NV_DATA_S(y);
    jac   = (double*) J->data;
    
    // Call fortran jacobian function
    // (jac is square so we can use N for neq)
    aq_integrate_jac(N, (double) t, state, jac, user_data);
    
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

static int test_Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
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

static void aq_integrate_print_output(realtype t, realtype y1, realtype y2, realtype y3)
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

static void aq_integrate_print_stats(void *cvode_mem, int n_eqn)
{
    long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
    int flag, i;
    
    realtype sug_tol;
    N_Vector err_weight, est_error;
    
    err_weight = N_VNew_Serial((long int) n_eqn);
    est_error = N_VNew_Serial((long int) n_eqn);
    
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumSteps", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumRhsEvals", 1);
    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumErrTestFails", 1);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
    
    flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVDlsGetNumJacEvals", 1);
    flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
    
    flag = CVodeGetNumGEvals(cvode_mem, &nge);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetNumGEvals", 1);

    flag = CVodeGetTolScaleFactor(cvode_mem, &sug_tol);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetTolScaleFactor", 1);
    flag = CVodeGetErrWeights(cvode_mem, err_weight);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetErrWeights", 1);
    flag = CVodeGetEstLocalErrors(cvode_mem, est_error);
    aq_integrate_SUNDIALS_check_flag(&flag, "CVodeGetEstLocalErrors", 1);

    printf("\nFinal Statistics:\n");
    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
           nst, nfe, nsetups, nfeLS, nje);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
           nni, ncfn, netf, nge);
    printf("suggested tolerance scaling factor = %lg\n\n", sug_tol);
    
    printf("SpecID\tErr Wt\tEst Error\tProduct\n");
    for (i=0; i<n_eqn; i++) {
        printf("%d\t%lg\t%lg\t%lg\n", i+1, Ith(err_weight,i+1), Ith(est_error,i+1),
               Ith(err_weight,i+1)*Ith(est_error,i+1));
    }
    
}




