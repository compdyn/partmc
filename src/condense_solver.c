/* Copyright (C) 2009-2012 Matthew West
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 */

/** \file
 * \brief Interface to SUNDIALS ODE solver library for condensation.
 */

#include <stdlib.h>
#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spgmr.h>

/** \brief Result code indicating successful completion.
 */
#define PMC_CONDENSE_SOLVER_SUCCESS        0
/** \brief Result code indicating failure to allocate \c y vector.
 */
#define PMC_CONDENSE_SOLVER_INIT_Y         1
/** \brief Result code indicating failure to allocate \c abstol vector.
 */
#define PMC_CONDENSE_SOLVER_INIT_ABSTOL    2
/** \brief Result code indicating failure to create the solver.
 */
#define PMC_CONDENSE_SOLVER_INIT_CVODE_MEM 3
/** \brief Result code indicating failure to initialize the solver.
 */
#define PMC_CONDENSE_SOLVER_INIT_CVODE     4
/** \brief Result code indicating failure to set tolerances.
 */
#define PMC_CONDENSE_SOLVER_SVTOL          5
/** \brief Result code indicating failure to set maximum steps.
 */
#define PMC_CONDENSE_SOLVER_SET_MAX_STEPS  6
/** \brief Result code indicating failure of the solver.
 */
#define PMC_CONDENSE_SOLVER_FAIL           7
/** \brief Result code indicating failure in creation of the linear solver.
 */
#define PMC_CONDENSE_SOLVER_LINSOL_CTOR    8
/** \brief Result code indicating failure in setting the linear solver.
 */
#define PMC_CONDENSE_SOLVER_LINSOL_SET     9
/** \brief Result code indicating failure in setting the preconditioner.
 */
#define PMC_CONDENSE_SOLVER_LINSOL_PREC    10

static int condense_vf(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int condense_check_flag(void *flagvalue, char *funcname, int opt);

/*******************************************************/
// solver block

static int condense_solver_Solve(double t, N_Vector ycur, N_Vector fcur,
				 N_Vector r, N_Vector b, double gamma, double delta, int lr, void *user_data);

/*******************************************************/

void condense_vf_f(int neq, realtype t, double *y_f, double *ydot_f);
void condense_jac_solve_f(int neq, double t, double *ycur_f, double *fcur_f, double *b_f, double gamma);

/** \brief Call the ODE solver.
 *
 * \param neq The number of equations.
 * \param x_f A pointer to a vector of \c neq variables, giving the
 * initial state vector on entry and the final state vector on exit.
 * \param abstol_f A pointer to a vector of \c neq variables, giving
 * the absolute error tolerance for the corresponding state vector
 * component.
 * \param reltol_f The scalar relative tolerance.
 * \param t_initial_f The initial time (s).
 * \param t_final_f The final time (s).
 * \return A result code (0 is success).
 */
int condense_solver(int neq, double *x_f, double *abstol_f, double reltol_f,
		    double t_initial_f, double t_final_f)
{
	realtype reltol, t_initial, t_final, t;
	N_Vector y, abstol;
	void *cvode_mem;
	int flag, i;
	realtype *y_data, *abstol_data;

	y = abstol = NULL;
	cvode_mem = NULL;

	y = N_VNew_Serial(neq);
	if (condense_check_flag((void *)y, "N_VNew_Serial", 0))
                return PMC_CONDENSE_SOLVER_INIT_Y;

	abstol = N_VNew_Serial(neq);
	if (condense_check_flag((void *)abstol, "N_VNew_Serial", 0))
                return PMC_CONDENSE_SOLVER_INIT_ABSTOL;

	y_data = NV_DATA_S(y);
	abstol_data = NV_DATA_S(abstol);
	for (i = 0; i < neq; i++) {
		y_data[i] = x_f[i];
		abstol_data[i] = abstol_f[i];
	}

	reltol = reltol_f;
	t_initial = t_initial_f;
	t_final = t_final_f;

	cvode_mem = CVodeCreate(CV_BDF);
	if (condense_check_flag((void *)cvode_mem, "CVodeCreate", 0))
                return PMC_CONDENSE_SOLVER_INIT_CVODE_MEM;

	flag = CVodeInit(cvode_mem, condense_vf, t_initial, y);
	if (condense_check_flag(&flag, "CVodeInit", 1))
                return PMC_CONDENSE_SOLVER_INIT_CVODE;

	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (condense_check_flag(&flag, "CVodeSVtolerances", 1))
                return PMC_CONDENSE_SOLVER_SVTOL;

	flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
	if (condense_check_flag(&flag, "CVodeSetMaxNumSteps", 1))
                return PMC_CONDENSE_SOLVER_SET_MAX_STEPS;


        SUNLinearSolver LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
        if (condense_check_flag((void *)LS, "SUNLinSol_SPGMR", 0))
                return PMC_CONDENSE_SOLVER_LINSOL_CTOR;
	flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
        if (condense_check_flag(&flag, "CVodeSetLinearSolver", 1))
                return PMC_CONDENSE_SOLVER_LINSOL_SET;
	flag = CVodeSetPreconditioner(cvode_mem, NULL, condense_solver_Solve);
        if (condense_check_flag(&flag, "CVodeSetPreconditioner", 1))
                return PMC_CONDENSE_SOLVER_LINSOL_PREC;

	t = t_initial;
	flag = CVode(cvode_mem, t_final, y, &t, CV_NORMAL);
	if (condense_check_flag(&flag, "CVode", 1))
                return PMC_CONDENSE_SOLVER_FAIL;

	for (i = 0; i < neq; i++) {
		x_f[i] = y_data[i];
	}

	N_VDestroy_Serial(y);
	N_VDestroy_Serial(abstol);
	CVodeFree(&cvode_mem);
	return PMC_CONDENSE_SOLVER_SUCCESS;
}

/** \brief The ODE vector field to integrate.
 *
 * \param t The current time (s).
 * \param y The state vector.
 * \param ydot The rate of change of the state vector.
 * \param user_data A pointer to user-provided data.
 * \return A result code (0 is success).
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
static int condense_vf(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype *y_data, *ydot_data;
	int i, neq;
	double *y_f, *ydot_f;

	neq = NV_LENGTH_S(y);
	y_data = NV_DATA_S(y);
	ydot_data = NV_DATA_S(ydot);

	y_f = malloc(neq * sizeof(double));
	ydot_f = malloc(neq * sizeof(double));

	for (i = 0; i < neq; i++) {
		y_f[i] = y_data[i];
	}
	condense_vf_f(neq, t, y_f, ydot_f);
	for (i = 0; i < neq; i++) {
		ydot_data[i] = ydot_f[i];
	}

	free(y_f);
	free(ydot_f);
	return(0);
}
#pragma GCC diagnostic pop

/** \brief Check the return value from a SUNDIALS call.
 *
 *  - <code>opt == 0</code> means SUNDIALS function allocates memory
 *            so check if \c flagvalue is not a \c NULL pointer.

 *  - <code>opt == 1</code> means SUNDIALS function returns a flag so
 *            check if the \c int pointed to by \c flagvalue has
 *            <code>flag >= 0</code>.

 *  - <code>opt == 2</code> means function allocates memory so check
 *            if \c flagvalue is not a \c NULL pointer.
 *
 * \param flagvalue A pointer to check (either for \c NULL, or as an
 * \c int pointer giving the flag value).
 * \param funcname A string giving the function name returning this
 * result code.
 * \param opt A flag indicating the type of check to perform.
 * \return A result code (0 is success).
 */
static int condense_check_flag(void *flagvalue, char *funcname, int opt)
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

/** \brief Linear solver routine for use by the ODE solver.
 *
 * Should solve the system \f$(I - \gamma J) x = b\f$, where \f$J\f$
 * is the current vector field Jacobian, \f$\gamma\f$ is a given
 * scalar, and \f$b\f$ is a given vector.
 *
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
static int condense_solver_Solve(double t, N_Vector ycur, N_Vector fcur,
				 N_Vector b, N_Vector z, double gamma, double delta, int lr, void *user_data)
{
	realtype *b_data, *ycur_data, *fcur_data, *z_data;
	int i, neq;
	double *b_f, *ycur_f, *fcur_f;

	neq = NV_LENGTH_S(b);
	b_data = NV_DATA_S(b);
	z_data = NV_DATA_S(z);
	ycur_data = NV_DATA_S(ycur);
	fcur_data = NV_DATA_S(fcur);

	b_f = malloc(neq * sizeof(double));
	ycur_f = malloc(neq * sizeof(double));
	fcur_f = malloc(neq * sizeof(double));

	for (i = 0; i < neq; i++) {
		b_f[i] = b_data[i];
		ycur_f[i] = ycur_data[i];
		fcur_f[i] = fcur_data[i];
	}
	condense_jac_solve_f(neq, t, ycur_f, fcur_f, b_f, gamma);
	for (i = 0; i < neq; i++) {
		z_data[i] = b_f[i];
	}

	free(b_f);
	free(ycur_f);
	free(fcur_f);
	return(0);
}
#pragma GCC diagnostic pop
