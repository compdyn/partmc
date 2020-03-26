/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file with some debugging functions for use with camp_solver.c
 *
 */
#ifdef PMC_DEBUG
#define PMC_DEBUG_SPEC_ 0
#define PMC_DEBUG_PRINT(x) \
  pmc_debug_print(sd->cvode_mem, x, false, 0, __LINE__, __func__)
#define PMC_DEBUG_PRINT_INT(x, y) \
  pmc_debug_print(sd->cvode_mem, x, false, y, __LINE__, __func__)
#define PMC_DEBUG_PRINT_FULL(x) \
  pmc_debug_print(sd->cvode_mem, x, true, 0, __LINE__, __func__)
#define PMC_DEBUG_JAC_STRUCT(J, x) pmc_debug_print_jac_struct((void *)sd, J, x)
#define PMC_DEBUG_JAC(J, x) pmc_debug_print_jac((void *)sd, J, x)
void pmc_debug_print(void *cvode_mem, const char *message, bool do_full,
                     const int int_val, const int line, const char *func) {
#ifdef PMC_USE_SUNDIALS
  CVodeMem cv_mem = (CVodeMem)cvode_mem;
  if (!(cv_mem->cv_debug_out)) return;
  printf(
      "\n[DEBUG] line %4d in %-20s(): %-25s %-4.0d t_n = %le h = %le q = %d "
      "hin = %le species %d(zn[0] = %le zn[1] = %le tempv = %le tempv1 = %le "
      "tempv2 = %le acor_init = %le last_yn = %le",
      line, func, message, int_val, cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_q,
      cv_mem->cv_hin, PMC_DEBUG_SPEC_,
      NV_DATA_S(cv_mem->cv_zn[0])[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_zn[1])[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_tempv)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_tempv1)[PMC_DEBUG_SPEC_],
      // NV_DATA_S(cv_mem->cv_tempv2)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_acor_init)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_last_yn)[PMC_DEBUG_SPEC_]);
  if (do_full) {
    for (int i = 0; i < NV_LENGTH_S(cv_mem->cv_y); i++) {
      printf(
          "\n  zn[0][%3d] = % -le zn[1][%3d] = % -le tempv[%3d] = % -le "
          "tempv1[%3d] = % -le tempv2[%3d] = % -le acor_init[%3d] = % -le "
          "last_yn[%3d] = % -le",
          i, NV_DATA_S(cv_mem->cv_zn[0])[i], i, NV_DATA_S(cv_mem->cv_zn[1])[i],
          i, NV_DATA_S(cv_mem->cv_tempv)[i], i, NV_DATA_S(cv_mem->cv_tempv1)[i],
          // i, NV_DATA_S(cv_mem->cv_tempv2)[i], i,
          NV_DATA_S(cv_mem->cv_acor_init)[i], i,
          NV_DATA_S(cv_mem->cv_last_yn)[i]);
    }
  }
#endif
}
void pmc_debug_print_jac_struct(void *solver_data, SUNMatrix J,
                                const char *message) {
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData *)solver_data;

  if (!(sd->debug_out)) return;
  int n_state_var = SM_COLUMNS_S(J);
  int i_elem = 0;
  int next_col = 0;
  printf("\n\n   Jacobian structure (↓ind →dep) - %s\n     ", message);
  for (int i_dep = 0; i_dep < n_state_var; i_dep++) printf("[%3d]", i_dep);
  for (int i_ind = 0; i_ind < n_state_var; i_ind++) {
    printf("\n[%3d]", i_ind);
    next_col = SM_INDEXPTRS_S(J)[i_ind + 1];
    for (int i_dep = 0; i_dep < n_state_var; i_dep++) {
      if (i_dep == SM_INDEXVALS_S(J)[i_elem] && i_elem < next_col) {
        printf(" %3d ", i_elem++);
      } else {
        printf("  -  ");
      }
    }
  }
#endif
}
void pmc_debug_print_jac(void *solver_data, SUNMatrix J, const char *message) {
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData *)solver_data;

  if (!(sd->debug_out)) return;
  int n_state_var = SM_COLUMNS_S(J);
  int i_elem = 0;
  int next_col = 0;
  printf("\n\n   Jacobian (↓ind →dep) - %s\n     ", message);
  for (int i_dep = 0; i_dep < n_state_var; i_dep++)
    printf("      [%3d]", i_dep);
  for (int i_ind = 0; i_ind < n_state_var; i_ind++) {
    printf("\n[%3d]   ", i_ind);
    next_col = SM_INDEXPTRS_S(J)[i_ind + 1];
    for (int i_dep = 0; i_dep < n_state_var; i_dep++) {
      if (i_dep == SM_INDEXVALS_S(J)[i_elem] && i_elem < next_col) {
        printf(" % -1.2le ", SM_DATA_S(J)[i_elem++]);
      } else {
        printf("     -     ");
      }
    }
  }
#endif
}

#else
#define PMC_DEBUG_PRINT(x)
#define PMC_DEBUG_PRINT_INT(x, y)
#define PMC_DEBUG_PRINT_FULL(x)
#define PMC_DEBUG_JAC_STRUCT(J, x)
#define PMC_DEBUG_JAC(J, x)
#endif

/** \brief Print some camp-chem data sizes
 *
 * \param md Pointer to the model data
 */
static void print_data_sizes(ModelData *md) {
  int *ptr = md->rxn_int_data;
  int n_rxn = ptr[0];

  printf("n_rxn: %d ", n_rxn);
  printf("n_state_var: %d", md->n_per_cell_state_var * md->n_cells);
  printf("n_dep_var: %d", md->n_per_cell_dep_var * md->n_cells);
}

/** \brief Print Jacobian matrix in format KLU SPARSE
 *
 * \param M Jacobian matrix
 */
static void print_jacobian(SUNMatrix J) {

  printf("\n NNZ JAC: %lld \n", SM_NNZ_S(J));
  printf("DATA | INDEXVALS:\n");
  for (int i = 0; i < SM_NNZ_S(J); i++) {
    printf("% -le, ", (SM_DATA_S(J))[i]);
    printf("%lld \n", (SM_INDEXVALS_S(J))[i]);
  }
  printf("PTRS:\n");
  for (int i = 0; i < SM_NP_S(J)+1; i++) {
    printf("%lld, ", (SM_INDEXPTRS_S(J))[i]);
  }
}

static void print_jacobian_file(SUNMatrix J, char *filepath) {

  // *A pointer to matrix coefficients
// *jA pointer to matrix colums
// *iA pointer to matrix row indexes
// nrows number of rows
// nnz   number of non zero values
  //void printMatrix(double* A, int* jA, int* iA, int nrows, int nnz)

  FILE *fp;
  fp= fopen("/gpfs/scratch/bsc32/bsc32815/gpupartmc/matrix_basic_1.csr","w");


  fprintf(fp," %d",SM_NNZ_S(J));
  fprintf(fp," %d",SM_NP_S(J));
  fprintf(fp," \n");

  for(int i=0;i<SM_NNZ_S(J);i++)
    fprintf(fp," %lf",SM_DATA_S(J)[i]);
  fprintf(fp," \n");
  for(int i=0;i<SM_NNZ_S(J);i++)
    fprintf(fp," %d",SM_INDEXVALS_S(J)[i]);
  fprintf(fp," \n");
  for(int i=0;i<=SM_NP_S(J);i++)
    fprintf(fp," %d",SM_INDEXPTRS_S(J)[i]);
  fprintf(fp," \n");

  fclose(fp);

}

/** \brief Print derivative array
 *
 * \param deriv Derivative array
 */
static void print_derivative(N_Vector deriv) {
  printf(" deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
    printf("% -le, ", NV_DATA_S(deriv)[i]);
    //printf(" index: %d \n", i);
  }
  printf("\n");

  /*FILE *fp;
  fp= fopen("/gpfs/scratch/bsc32/bsc32815/gpupartmc/rhs_basic2_1_le.csr","w");

  //fprintf(fp," deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
    fprintf(fp," %-le", NV_DATA_S(deriv)[i]);
    //printf(" index: %d \n", i);
  }
  fprintf(fp," \n");

  fclose(fp);
*/
}
