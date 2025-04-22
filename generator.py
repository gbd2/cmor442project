import numpy as np
import pandas as pd
from itertools import product

def lognormal_params(mean, std):
    """
    Converts mean and standard deviation to log-normal parameters.
    The log-normal distribution is defined by its mean and standard deviation
    in the log space.
    
    Parameters:
        mean (float): Desired mean of the log-normal distribution.
        std (float): Desired standard deviation of the log-normal distribution.
    Returns:
        mu_n (float): Mean of the underlying normal distribution.
        sigma_n (float): Standard deviation of the underlying normal distribution.
    """
    sigma_n_sq = np.log(1 + (std**2) / (mean**2))
    mu_n = np.log(mean) - sigma_n_sq / 2
    sigma_n = np.sqrt(sigma_n_sq)
    return mu_n, sigma_n


def generate_joint_log_normal_samples(mu1, sigma1, mu2, sigma2, N):
    """
    Generates N joint log-normal samples for DEMAND1 and DEMAND2
    with specified target means and standard deviations.

    Parameters:
            mu1 (float): Mean of DEMAND1.
            sigma1 (float): Standard deviation of DEMAND1.
            mu2 (float): Mean of DEMAND2.
            sigma2 (float): Standard deviation of DEMAND2.
            N (int): Number of samples to generate. 
    Returns:
        samples: list of dicts with keys 'DEMAND1', 'DEMAND2', 'prob'
    """
    
    mu1_n, sigma1_n = lognormal_params(mu1, sigma1)
    mu2_n, sigma2_n = lognormal_params(mu2, sigma2)

    # Sample from independent normal distributions
    demand1_vals = np.random.normal(mu1_n, sigma1_n, size=N)
    demand2_vals = np.random.normal(mu2_n, sigma2_n, size=N)

    demand1 = np.exp(demand1_vals)
    demand2 = np.exp(demand2_vals)

    # Combine into structured format with equal probabilities
    prob = 1.0 / N
    samples = [{"DEMAND1": d1, "DEMAND2": d2, "prob": prob} for d1, d2 in zip(demand1, demand2)]

    return samples

def generate_indep_log_normal_samples(mu1, sigma1, mu2, sigma2, N1, N2):
    """
    Generates N1 independent samples for DEMAND1 and N2 for DEMAND2
    using log-normal distributions with specified means and stds.

    Parameters:
        mu1 (float): Mean of DEMAND1.
        sigma1 (float): Standard deviation of DEMAND1.
        mu2 (float): Mean of DEMAND2.
        sigma2 (float): Standard deviation of DEMAND2.
        N1 (int): Number of samples for DEMAND1.
        N2 (int): Number of samples for DEMAND2.
    Returns:
        scenarios: list of dicts with keys 'DEMAND1', 'DEMAND2', 'prob'
    """

    mu1_n, sigma1_n = lognormal_params(mu1, sigma1)
    mu2_n, sigma2_n = lognormal_params(mu2, sigma2)

    # Sample independently
    demand1_vals = np.random.normal(mu1_n, sigma1_n, size=N1)
    demand2_vals = np.random.normal(mu2_n, sigma2_n, size=N2)
    
    demand1 = np.exp(demand1_vals)
    demand2 = np.exp(demand2_vals)    

    prob = 1.0 / (N1 * N2)
    scenarios = [{"DEMAND1": d1, "DEMAND2": d2, "prob": prob} for d1, d2 in product(demand1, demand2)]

    return scenarios
