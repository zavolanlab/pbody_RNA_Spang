import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.stats import entropy
from .simulation import simulate_log2_nCPM_from_one_condition
import logging

def kl_divergence_hist(samples_p, samples_q, bins=50):
    # the higher the number of bins -> the higher will be KL divergence value
    # Determine common range for both histograms
    # samples_p represents reference, while samples_q represents query
    
    combined = np.concatenate([samples_p, samples_q])
    bin_edges = np.linspace(np.min(combined), np.max(combined), bins + 1)
    
    # Calculate histograms (densities)
    p_hist, _ = np.histogram(samples_p, bins=bin_edges, density=True)
    q_hist, _ = np.histogram(samples_q, bins=bin_edges, density=True)
    
    # Avoid zeros
    p_hist += 1e-10
    q_hist += 1e-10
    
    # Normalize to ensure sum is 1
    p_hist /= p_hist.sum()
    q_hist /= q_hist.sum()
    kl_div = entropy(p_hist, q_hist)
    
    return kl_div

def fit_size_factors(sf_values):
    """Fits LogNorm to provided size factors."""
    s, loc, scale = stats.lognorm.fit(sf_values, floc=0)
    return [s, loc, scale]

class Fitter:
    def __init__(self, obs_total, obs_frac, sf_params, n_replicates=3):
        self.obs_total = obs_total
        self.obs_frac = obs_frac
        self.n_g = len(obs_total)
        self.n_r = n_replicates
        self.sf_params = sf_params # s, loc, scale
        
        # Placeholders for fitted params
        self.tot_params = None
        self.frac_params = None

    # KL divergence for Total Expression
    def _KLdiv_total_total(self, TotalExpr_params, Nsimul=20,return_list_of_vals=False):
        
        # Simulation
        kl_vals = []
        for _ in range(Nsimul): # how many times to simulate
            sim_log2_norm_expr_df = simulate_log2_nCPM_from_one_condition(
                TotalExpr_params,
                self.sf_params,
                N_g=self.n_g, N_r=self.n_r
            )
            kl = kl_divergence_hist(self.obs_total, sim_log2_norm_expr_df['log2_nCPM_Total'])
            kl_vals.append(kl)
        if return_list_of_vals:
            return kl_vals
        else:
            return np.mean(kl_vals)

    def optimize_total(self, n_starts=10):
        # Bounds: df, loc, scale
        bounds = np.array([[1, 10], [1, 12], [0, 12]])
        best_res = None
        
        for i in range(n_starts): # Multi-start local optimization

            if sys.version_info > (3, 0):
                # Check for interruption between heavy optimization starts
                import signal
                if signal.getsignal(signal.SIGINT) == signal.default_int_handler:
                    # Only checks if the handler is default
                    pass
            
            x0 = np.random.uniform(bounds[:, 0], bounds[:, 1])
            res = minimize(self._KLdiv_total_total, x0, bounds=bounds, 
                           method='Nelder-Mead', tol=0.1, 
                          options={'adaptive': True,'maxiter':20})
            if best_res is None or res.fun < best_res.fun:
                best_res = res
            logging.info(f"{str(i+1)} out of {n_starts} starts of NM optimization done, curr best KLdiv is {str(np.round(best_res.fun,4))}")
            
        self.tot_params = best_res.x
        return best_res.x
