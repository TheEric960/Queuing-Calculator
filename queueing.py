"""
This script computes queueing results for simple queueing models
described in Chapter 6 of Banks, Carson, Nelson and Nicol,
Discrete-Event System Simulation, 5th edition.

Original code courtesy of Professor Barry Nelson
Python conversion by Jack Wilson
"""

import math


def lmda_ge_mu(lmda, mu):
    """Check if lambda is greater than equal to mu."""
    if lmda >= mu:
        raise ValueError('The arrival rate must be less than the service rate')


def lmda_mu_le_0(lmda, mu):
    """Check if lamda or mu is less than equal to 0."""
    if lmda <= 0 or mu <= 0:
        raise ValueError('Arrival rate and service rate must be positive')


def lmda_mu_lt_0(lmda, mu):
    """Check if lambda or mu is less than 0."""
    if lmda < 0 or mu < 0:
        raise ValueError('Arrival rate and service rate must be positive')


def sigma2_lt_0(sigma2):
    """Check if the sigma squared is less than 0."""
    if sigma2 < 0:
        raise ValueError('Variance must be nonnegative')


def lmda_ge_cmu(lmda, c, mu):
    """Check if lambda is greater than equal to c times mu."""
    if lmda >= c * mu:
        raise ValueError('The arrival rate must be less than c times the service rate')


def lmda_et_cmu(lmda, c, mu):
    """Check if lambda is equal to c times mu."""
    if lmda == c * mu:
        raise ValueError('This spreadsheet does not handle the case Lambda equal c*Mu')


def c_lt_1(c):
    """Check if c is less than one."""
    if c < 1:
        raise ValueError('Number of servers must be positive')


def n_lt_c(n, c):
    """Check if n is less than c."""
    if n < c:
        raise ValueError('Capacity must be at least as large as the number of servers')


def k_lt_c(k, c):
    """Check if k is less than c."""
    if k < c:
        raise ValueError('Size of calling population must be at least as large as the number of servers')


def eval_MG1(lmda: float = 0, mu: float = 0, sigma2: float = 0):
    """Evaluate an M/G/1 queue."""
    lmda_ge_mu(lmda, mu)
    lmda_mu_le_0(lmda, mu)
    sigma2_lt_0(sigma2)

    rho = lmda / mu
    l = rho + (rho ** 2 * (1 + sigma2 * mu ** 2)) / 2 / (1 - rho)
    w = l / lmda
    wq = w - 1 / mu
    lq = wq * lmda
    p0 = 1 - rho

    return rho, l, w, wq, lq, p0


def eval_MMc(lmda: float = 0, mu: float = 0, c: int = 0):
    """Evaluate an M/M/c queue."""
    lmda_ge_cmu(lmda, c, mu)
    lmda_mu_le_0(lmda, mu)
    c_lt_1(c)

    rho = lmda / mu / c
    offered_load = lmda / mu
    factor = 1
    p0 = 1

    for i in range(1, c):
        factor *= offered_load / i
        p0 += factor

    cfactorial = math.factorial(c)
    p0 += factor * offered_load / c / (1 - rho)
    p0 = 1 / p0

    l = offered_load + (offered_load ** (c + 1) * p0) / c / cfactorial / (1 - rho) ** 2
    w = l / lmda
    wq = w - 1 / mu
    lq = wq * lmda

    return rho, l, w, wq, lq, p0


def eval_MGc(lmda: float = 0, mu: float = 0, sigma2: float = 0, c: int = 0):
    """Evaluate an M/G/c queue."""
    lmda_ge_cmu(lmda, c, mu)
    lmda_mu_le_0(lmda, mu)
    c_lt_1(c)
    sigma2_lt_0(sigma2)

    cv2 = sigma2 * mu ** 2
    rho = lmda / mu / c
    offered_load = lmda / mu
    factor = 1
    p0 = 1

    for i in range(1, c):
        factor *= offered_load / i
        p0 += factor

    cfactorial = math.factorial(c)
    p0 += factor * offered_load / c / (1 - rho)
    p0 = 1 / p0

    lq = (offered_load ** (c + 1) * p0 / c / cfactorial / (1 - rho) ** 2) * (1 + cv2) / 2
    wq = lq / lmda
    w = wq + 1 / mu
    l = w * lmda

    return rho, l, w, wq, lq


def eval_MMcN(lmda: float = 0, mu: float = 0, c: int = 0, n: int = 0):
    """Evaluate an M/M/c/N queue."""
    lmda_mu_lt_0(lmda, mu)
    c_lt_1(c)
    n_lt_c(n, c)
    lmda_et_cmu(lmda, c, mu)

    rho = lmda / mu / c
    offered_load = lmda / mu
    factor = 1
    p0 = 1

    for i in range(1, c + 1):
        factor *= offered_load / i
        p0 += factor

    if c < n:
        rhosum = rho
        if c < n + 1:
            for i in range(c + 2, n + 1):
                rhosum += rho ** (i - c)
        p0 += factor * rhosum

    cfactorial = math.factorial(c)
    p0 = 1 / p0
    pN = (offered_load ** n / cfactorial / c ** (n - c)) * p0
    lq = p0 * offered_load ** c * rho / cfactorial / (1 - rho) ** 2 * (
            1 - rho ** (n - c) - (n - c) * rho ** (n - c) * (1 - rho)
    )
    lmda_effective = lmda * (1 - pN)
    wq = lq / lmda_effective
    w = wq + 1 / mu
    l = lmda_effective * w
    rho = lmda_effective / c / mu

    return rho, l, w, wq, lq, p0, pN, lmda_effective


def eval_MMcKK(lmda: float = 0, mu: float = 0, c: int = 0, k: int = 0):
    """Evaluate an M/M/c/K/K queue."""
    lmda_mu_lt_0(lmda, mu)
    c_lt_1(c)
    k_lt_c(k, c)

    p = [0] * (k + 1)
    offered_load = lmda / mu
    kfac = math.factorial(k)
    p0 = 1

    if c > 1:
        for i in range(1, c):
            p[i] = (kfac / math.factorial(i) / math.factorial(k - i)) * offered_load ** i
            p0 += p[i]

    cfactorial = math.factorial(c)

    for i in range(c, k + 1):
        p[i] = (kfac / math.factorial(k - i) / cfactorial / c ** (i - c)) * offered_load ** i
        p0 += p[i]

    p0 = 1 / p0
    l = 0
    lq = 0
    lmda_effective = k * lmda * p0

    for i in range(1, k + 1):
        p[i] *= p0
        l += i * p[i]
        lq += max(0, i - c) * p[i]
        lmda_effective += lmda * (k - i) * p[i]

    w = l / lmda_effective
    wq = lq / lmda_effective
    rho = lmda_effective / c / mu

    return rho, l, w, wq, lq, p0, lmda_effective
