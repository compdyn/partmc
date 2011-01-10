#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng *pmc_rand_gsl_rng = NULL;

#define PMC_RAND_GSL_SUCCESS      0
#define PMC_RAND_GSL_INIT_FAIL    1
#define PMC_RAND_GSL_NOT_INIT     2
#define PMC_RAND_GSL_ALREADY_INIT 3

int pmc_srand_gsl(int seed)
{
        if (pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_ALREADY_INIT;
        }
        gsl_set_error_handler_off(); // turn off automatic error handling
        pmc_rand_gsl_rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (pmc_rand_gsl_rng == NULL) {
                return PMC_RAND_GSL_INIT_FAIL;
        }
        gsl_rng_set(pmc_rand_gsl_rng, seed);
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_finalize_gsl()
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        gsl_rng_free(pmc_rand_gsl_rng);
        pmc_rand_gsl_rng = NULL;
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_gsl(double *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_rng_uniform(pmc_rand_gsl_rng);
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_int_gsl(int n, int *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_rng_uniform_int(pmc_rand_gsl_rng, n) + 1;
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_normal_gsl(double mean, double stddev, double *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_ran_gaussian(pmc_rand_gsl_rng, stddev) + mean;
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_poisson_gsl(double mean, int *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_ran_poisson(pmc_rand_gsl_rng, mean);
        return PMC_RAND_GSL_SUCCESS;
}

int pmc_rand_binomial_gsl(int n, double p, int *harvest)
{
        unsigned int u;

        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        u = n;
        *harvest = gsl_ran_binomial(pmc_rand_gsl_rng, p, u);
        return PMC_RAND_GSL_SUCCESS;
}
