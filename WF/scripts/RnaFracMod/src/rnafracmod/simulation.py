import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import logging
from . import processing

def quantile_normalize_to_reference(target, reference):
    """
    Normalizes 'target' so it matches the distribution of 'reference'.
    """
    # 1. Sort the reference to get the desired distribution curve
    sorted_reference = np.sort(reference)
    
    # 2. Get the rank of the target values
    # argsort twice gives the rank (0, 1, ... N-1) of each element
    rank_indices = np.argsort(np.argsort(target))
    
    # 3. Calculate percentiles (0.0 to 1.0) for the target ranks
    # Use (N - 1) to ensure the range covers 0.0 to 1.0 inclusive
    target_percentiles = rank_indices / (len(target) - 1)
    
    # 4. Calculate percentiles for the reference data
    reference_percentiles = np.linspace(0, 1, len(reference))
    
    # 5. Interpolate
    # Find what value in the sorted reference corresponds to the target's percentile
    normalized_target = np.interp(target_percentiles, reference_percentiles, sorted_reference)
    
    return normalized_target

def hybrid_poisson_sampler(log_lambdas, cutoff=5):
    """
    Samples from a vector of log2 logarithms of lambdas using Poisson for small values
    and a Normal approximation for large values.
    """
    log_lambdas = np.asanyarray(log_lambdas)
    samples = np.zeros_like(log_lambdas)

    # --- Clamp values to prevent Overflow ---
    # 2^700 is safely within float64 limits (~1.8e308 is 2^1024)
    # This prevents the "RuntimeWarning: overflow encountered in power"
    log_lambdas = np.clip(log_lambdas, a_min=None, a_max=700)
    
    # Create masks for small and large lambdas
    is_large = log_lambdas > cutoff
    is_small = ~is_large
    
    # 1. Poisson sampling for small lambdas
    if np.any(is_small):
        samples[is_small] = np.random.poisson(2**log_lambdas[is_small])
        
    # 2. Normal approximation for large lambdas: N(mu=lam, sigma=sqrt(lam))
    if np.any(is_large):
        mu = 2**log_lambdas[is_large]
        sigma = np.sqrt(mu)
        # Use round to keep results as integers (Poisson properties)
        # and clip to 0 to avoid rare negative values from Normal tails
        norm_samples = np.random.normal(loc=mu, scale=sigma)
        samples[is_large] = np.maximum(0, np.round(norm_samples))
        
    return samples

def simulate_log2_nCPM_from_one_condition(
                        TotalExpr_params=[2,5,3],
                       SF_params = [1,0,1],
                       N_r = 3, N_g = 7000, log2_f_alpha=0,log2_f_beta=0,
                    debug=True):
    TotExpr_df, log_TotExpr_loc, log_TotExpr_scale = TotalExpr_params
    SF_s, SF_loc, SF_scale = SF_params
    # N_r is a number of replicate samples per condition
    # N_g is a number of genes that are truly expressed
    # f_alpha and f_beta define the parameters of Beta distribution specifying the proportion of each 
    # molecule in a sub-fraction
    
    log2_r_1 = stats.t.rvs(df = TotExpr_df,loc=log_TotExpr_loc,scale=log_TotExpr_scale, size=N_g) # distribution of molecule counts per gene
    # r_1 = np.round(2**log2_r_1,0)+1 # need at least one molecule per gene, make discrete

    f_alpha,f_beta = 1+2**log2_f_alpha,1+2**log2_f_beta # during optimization we don't use constraints, but we want values below 1 here
    p_1 = stats.beta.rvs(a=f_alpha, b=f_beta, size=N_g) # proportion of molecules in a sub-fraction, condition 1, aka recruitment efficiency (RE)!
    min_positive_p = np.min(p_1[p_1>0])
    p_1 = p_1+min_positive_p # like pseudocount
    
    # f_1 = np.round(p_1*r_1,0) # True number of molecules in a fraction, condition 1
    log2_f_1 = np.log2(p_1)+log2_r_1
    
    # here, proportion in a sub-fraction, and total counts are assumed to be independent
    # so, we don't expect e.g. that higher-expressed genes would have higher RE
    
    # now, we do Poisson sampling for defined number of bioreplicates
    expr_data,cols = [],[]
    total_c1_samples,frac_c1_samples = [],[]
    for repl in range(N_r):
        # true log2 size factor is a random variable following Normal Distrib
        sf_1,f_sf_1 = stats.lognorm.rvs(s = SF_s, loc=SF_loc, scale=SF_scale, size=2) # _1 specifies condition 1
        log_sf_1,log_f_sf_1 = np.log(sf_1),np.log(f_sf_1)
        if debug:
            pass
            # print('log size factors: '+str([log_sf_1,log_f_sf_1]))
            
        log_lam_1,log_f_lam_1 = log2_r_1+log_sf_1,log2_f_1+log_f_sf_1 # lambda parameters for count sampling
        if debug:
            pass
            # print('log lambda vals Total: '+str([min(log_lam_1),max(log_lam_1)]))
            # print('log lambda vals SubFrac: '+str([min(log_f_lam_1),max(log_f_lam_1)]))
        
        expr_data.append(hybrid_poisson_sampler(log_lam_1))
        expr_data.append(hybrid_poisson_sampler(log_f_lam_1))
    
        cols = cols + [str(repl)+'_c1_Total',str(repl)+'_c1_Frac']
        total_c1_samples.append(str(repl)+'_c1_Total')
        frac_c1_samples.append(str(repl)+'_c1_Frac')
        
    expr_df = pd.DataFrame(expr_data).transpose()
    expr_df.columns = cols
    
    # do normalization like Deseq2
    norm_expr_df, sfs_df = processing.deseq2_normalize(expr_df, cols)
    log2_norm_expr_df = np.log2(norm_expr_df) # should work because of pseudocount

    log2_norm_expr_df['log2_nCPM_Total'] = log2_norm_expr_df[total_c1_samples].mean(1)
    log2_norm_expr_df['log2_nCPM_SubFraction'] = log2_norm_expr_df[frac_c1_samples].mean(1)
    
    return log2_norm_expr_df