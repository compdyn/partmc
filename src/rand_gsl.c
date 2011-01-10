/* Copyright (C) 2010 Matthew West
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 */

/** \file
 * \brief Wrapper routines for GSL random number functions.
 */

#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/** \brief Private internal-use variable to store the random number
 * generator.
 */
static gsl_rng *pmc_rand_gsl_rng = NULL;

/** \brief Result code indicating successful completion.
 */
#define PMC_RAND_GSL_SUCCESS      0
/** \brief Result code indicating initialization failure.
 */
#define PMC_RAND_GSL_INIT_FAIL    1
/** \brief Result code indicating the generator was not initialized
 * when it should have been.
 */
#define PMC_RAND_GSL_NOT_INIT     2
/** \brief Result code indicating the generator was already
 * initialized when an initialization was attempted.
 */
#define PMC_RAND_GSL_ALREADY_INIT 3

/** \brief Initialize the random number generator with the given seed.
 *
 * This must be called before any other GSL random number functions
 * are called.
 *
 * \param seed The random seed to use.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 * \sa pmc_rand_finalize_gsl() to cleanup the generator.
 */
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

/** \brief Cleanup and deallocate the random number generator.
 *
 * This must be called after pmc_srand_gsl().
 *
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
int pmc_rand_finalize_gsl()
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        gsl_rng_free(pmc_rand_gsl_rng);
        pmc_rand_gsl_rng = NULL;
        return PMC_RAND_GSL_SUCCESS;
}

/** \brief Generate a uniform random number in \f$[0,1)\f$.
 *
 * \param harvest A pointer to the generated random number.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
int pmc_rand_gsl(double *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_rng_uniform(pmc_rand_gsl_rng);
        return PMC_RAND_GSL_SUCCESS;
}

/** \brief Generate a uniform random integer in \f$[1,n]\f$.
 *
 * \param n The upper limit of the random integer.
 * \param harvest A pointer to the generated random number.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
int pmc_rand_int_gsl(int n, int *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_rng_uniform_int(pmc_rand_gsl_rng, n) + 1;
        return PMC_RAND_GSL_SUCCESS;
}

/** \brief Generate a normally-distributed random number.
 *
 * \param mean The mean of the distribution.
 * \param stddev The standard deviation of the distribution.
 * \param harvest A pointer to the generated random number.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
int pmc_rand_normal_gsl(double mean, double stddev, double *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_ran_gaussian(pmc_rand_gsl_rng, stddev) + mean;
        return PMC_RAND_GSL_SUCCESS;
}

/** \brief Generate a Poisson-distributed random integer.
 *
 * \param mean The mean of the distribution.
 * \param harvest A pointer to the generated random number.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
int pmc_rand_poisson_gsl(double mean, int *harvest)
{
        if (!pmc_rand_gsl_rng) {
                return PMC_RAND_GSL_NOT_INIT;
        }
        *harvest = gsl_ran_poisson(pmc_rand_gsl_rng, mean);
        return PMC_RAND_GSL_SUCCESS;
}

/** \brief Generate a Binomial-distributed random integer.
 *
 * \param n The sample size for the distribution.
 * \param p The sample probability for the distribution.
 * \param harvest A pointer to the generated random number.
 * \return PMC_RAND_GSL_SUCCESS on success, otherwise an error code.
 */
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
