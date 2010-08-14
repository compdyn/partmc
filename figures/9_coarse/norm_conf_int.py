#!/usr/bin/env python

import math
import numpy
import scipy.stats
import scipy.optimize

def same_sign(a, b):
    if (a < 0 and b < 0) or (a > 0 and b > 0):
        return True
    else:
        return False

def samp_mean_conf_int(n, X, S, cl):
    """Compute a confidence interval for the sample mean.

    Parameters:
    n : number of samples
    X : sample mean = sum(data) / n
    S : sample standard deviation = sqrt(sum((data - X)**2) / (n - 1))
    cl : confidence level (e.g. 0.95 for 95% confidence interval)

    Returns: (X_l, X_u)
    [X_l, X_u] : confidence interval for X"""
    
    if cl < 0 or cl > 1:
        raise ValueError("confidence level cl must be between 0 and 1")
    if n < 2:
        raise ValueError("number of samples n must be at least 2")
    if S <= 0:
        raise ValueError("sample standard deviation S must be positive")

    r = scipy.stats.t.ppf((1 + cl) / 2, n - 1)
    X_l = X - r * S
    X_u = X + r * S
    return (X_l, X_u)

def samp_var_conf_int(n, S2, cl):
    """Compute a confidence interval for the sample variance.

    Parameters:
    n : number of samples
    S2 : sample variance = sum((data - X)**2) / (n - 1)
    cl : confidence level (e.g. 0.95 for 95% confidence interval)

    Returns: (S2_l, S2_u)
    [S2_l, S2_u] : confidence interval for S2"""

    if cl < 0 or cl > 1:
        raise ValueError("confidence level cl must be between 0 and 1")
    if n < 2:
        raise ValueError("number of samples n must be at least 2")
    if S2 <= 0:
        raise ValueError("sample standard variance S2 must be positive")

    df = n - 1 # degrees-of-freedom
    S2_l = df * S2 / scipy.stats.chi2.ppf(0.5 + cl/2, df)
    S2_u = df * S2 / scipy.stats.chi2.ppf(0.5 - cl/2, df)
    return (S2_l, S2_u)

def samp_std_conf_int(n, S, cl):
    """Compute a confidence interval for the sample standard
    deviation.

    Parameters:
    n : number of samples
    S : sample standard deviation = sqrt(sum((data - X)**2) / (n - 1))
    cl : confidence level (e.g. 0.95 for 95% confidence interval)

    Returns: (S_l, S_u)
    [S_l, S_u] : confidence interval for S"""

    if cl < 0 or cl > 1:
        raise ValueError("confidence level cl must be between 0 and 1")
    if n < 2:
        raise ValueError("number of samples n must be at least 2")
    if S <= 0:
        raise ValueError("sample standard deviation S must be positive")

    (S2_l, S2_u) = samp_var_conf_int(n, S**2, cl)
    return (math.sqrt(S2_l), math.sqrt(S2_u))

def cv_conf_int(n, cv, cl):
    """Compute a confidence interval for the coefficient of variation.

    Parameters:
    n : number of samples
    cv : coefficient of variation = sample_std / sample_mean
    cl : confidence level (e.g. 0.95 for 95% confidence interval)

    Returns: (cv_l, cv_u)
    [cv_l, cv_u] : confidence interval for cv

    References:

    N. L. Johnson and B. L. Welch, "Applications of the non-central
    t-distribution", Biometrika 31(3-4): 362, 1940.

    S. Verrill, "Confidence Bounds for Normal and Lognormal
    Distribution Coefficients of Variation", Research Paper 609, USDA
    Forest Products Laboratory, 2003.

    S. Verrill, "Confidence Bounds and Hypothesis Tests for Normal
    Distribution Coefficients of Variation", Research Paper 638, USDA
    Forest Products Laboratory, 2007.

    S. Verrill and R. A. Johnson, "Confidence Bounds and Hypothesis
    Tests for Normal Distribution Coefficients of Variation",
    Communications in Statistics -- Theory and Methods, 36(12):
    2187-2206, 2007.

    M. G. Vangel, "Confidence Intervals for a Normal Coefficient of
    Variation", The American Statistician 50(1), 1996.
    """

    if cl < 0 or cl > 1:
        raise ValueError("confidence level cl must be between 0 and 1")
    if n < 2:
        raise ValueError("number of samples n must be at least 2")
    if cv == 0:
        raise ValueError("coefficient of variation cv must be non-zero")

    if cv < 0:
        (neg_int_l, neg_int_u) = cv_conf_int(n, -cv, cl)
        return (-neg_int_u, -neg_int_l)

    x = math.sqrt(n) / cv
    df = n - 1 # degrees-of-freedom

    def cv_l_fcn(beta):
        nc = math.sqrt(n) / beta # non-centrality parameter
        return (0.5 - cl/2) - scipy.stats.nct.cdf(x, df, nc)

    def cv_u_fcn(beta):
        nc = math.sqrt(n) / beta # non-centrality parameter
        return (0.5 + cl/2) - scipy.stats.nct.cdf(x, df, nc)

    cv_l_init = cv
    for i in range(10):
        cv_l_init /= 10
        if not same_sign(cv_l_fcn(cv_l_init), cv_l_fcn(cv)):
            break
    else:
        raise ValueError("unable to determine initial conditions for lower-bound search")
    cv_l = scipy.optimize.brentq(cv_l_fcn, cv_l_init, cv)

    limit_cv_u_fcn = (0.5 + cl/2) - scipy.stats.t.cdf(x, df, 0)
    if same_sign(cv_u_fcn(cv), limit_cv_u_fcn):
        cv_u = numpy.inf
    else:
        cv_u_init = cv
        for i in range(10):
            cv_u_init *= 10
            if not same_sign(cv_u_fcn(cv_u_init), cv_u_fcn(cv)):
                break
        else:
            raise ValueError("unable to determine initial conditions for upper-bound search")
        cv_u = scipy.optimize.brentq(cv_u_fcn, cv, cv_u_init)

    return (cv_l, cv_u)

if __name__ == "__main__":
    cl = 0.95
    n = 100

    X = 10
    S = 2
    (X_l, X_u) = samp_mean_conf_int(n, X, S, cl)
    print "number of samples n = {0}".format(n)
    print "sample mean X = {0}".format(X)
    print "sample standard deviation S = {0}".format(S)
    print "X {0}% confidence interval = [{1}, {2}]".format(cl * 100, X_l, X_u)

    (S_l, S_u) = samp_std_conf_int(n, S, cl)
    print
    print "number of samples n = {0}".format(n)
    print "sample standard deviation S = {0}".format(S)
    print "S {0}% confidence interval = [{1}, {2}]".format(cl * 100, S_l, S_u)

    S2 = S**2
    (S2_l, S2_u) = samp_var_conf_int(n, S2, cl)
    print
    print "number of samples n = {0}".format(n)
    print "sample standard variance S2 = {0}".format(S2)
    print "S2 {0}% confidence interval = [{1}, {2}]".format(cl * 100, S2_l, S2_u)

    for cv in [0.2, -0.2, 6.0]:
        (cv_l, cv_u) = cv_conf_int(n, cv, cl)
        print
        print "number of samples n = {0}".format(n)
        print "coefficent of variation CV = {0}".format(cv)
        print "CV {0}% confidence interval = [{1}, {2}]".format(cl * 100, cv_l, cv_u)
