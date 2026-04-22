import argparse
import pandas as pd
import numpy as np
import sys
import os
import logging

from scipy import stats

from . import reporting
from . import fitting
from . import simulation
from . import processing

def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
    """
    Catch all crashes and log them to file before exiting.
    """
    if issubclass(exc_type, KeyboardInterrupt):
        # Allow Ctrl+C to stop the program without logging a crash
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    # Log the full traceback to the file
    logging.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))

# Register the hook
sys.excepthook = handle_uncaught_exception

def setup_logging(outdir, log_file_arg=None):
    """
    Configures logging to output to both a file and the console.
    """
    # Determine log path
    if log_file_arg:
        log_path = log_file_arg
        # Ensure directory for explicit log file exists
        log_dir = os.path.dirname(log_path)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
    else:
        log_path = os.path.join(outdir, "rnafracmod.log")

    # Create handlers
    file_handler = logging.FileHandler(log_path, mode='w')
    console_handler = logging.StreamHandler(sys.stdout)

    # set formatting
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Apply configuration
    logging.basicConfig(
        level=logging.INFO,
        handlers=[file_handler, console_handler],
        force=True # Reset any existing handlers
    )
    return log_path

def validate_metadata(meta):
    """
    Validates the metadata dataframe structure and content.
    """
    # 1. Check for required columns
    required_cols = {'sample', 'EXP_CTL', 'fraction'}
    if not required_cols.issubset(meta.columns):
        missing = required_cols - set(meta.columns)
        raise ValueError(f"Metadata is missing required columns: {missing}")

    # 2. Check for duplicate sample IDs
    if meta['sample'].duplicated().any():
        dups = meta.loc[meta['sample'].duplicated(), 'sample'].unique()
        raise ValueError(f"Metadata contains duplicate sample IDs: {dups}")

    # 3. Validate 'EXP_CTL' (Conditions)
    valid_conditions = {'EXP', 'CTL'}
    unique_conditions = set(meta['EXP_CTL'].unique())
    if not unique_conditions.issubset(valid_conditions):
        raise ValueError(f"Column 'EXP_CTL' contains invalid values. Allowed: {valid_conditions}. Found: {unique_conditions}")
    
    # 4. Validate 'fraction'
    valid_fractions = {'T', 'F'}
    unique_fractions = set(meta['fraction'].unique())
    if not unique_fractions.issubset(valid_fractions):
        raise ValueError(f"Column 'fraction' contains invalid values. Allowed: {valid_fractions}. Found: {unique_fractions}")

    # 5. Ensure every condition has both 'T' and 'F' fractions
    for condition in unique_conditions:
        # Get all fractions associated with this specific condition
        fractions_in_condition = set(meta.loc[meta['EXP_CTL'] == condition, 'fraction'])
        
        # Check if both 'T' and 'F' are present
        if not {'T', 'F'}.issubset(fractions_in_condition):
            raise ValueError(f"Condition '{condition}' is incomplete. It must contain at least one Total ('T') "
                             f"and at least one Fraction ('F') sample. Found: {fractions_in_condition}")
    return True

def main():
    parser = argparse.ArgumentParser(description="RNA Fractionation Modelling Tool")
    parser.add_argument('--count_matrix', required=True, help="Path to TSV count matrix (genes x samples)")
    parser.add_argument('--metadata', required=True, help="Path to sample metadata TSV")
    parser.add_argument('--outdir', required=True, help="Output directory")
    parser.add_argument('--log_file', help="Explicit path to log file. Defaults to <outdir>/rnafracmod.log")
    parser.add_argument('--plot_diagnostics', action='store_true', 
                        help="If set, generates diagnostic plots in <outdir>/plots/")
    
    args = parser.parse_args()
    
    # --- 0.1 Setup Output Directory & Logging ---
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        
    log_path = setup_logging(args.outdir, args.log_file)
    logging.info(f"Starting analysis. Log file: {log_path}")

    # 0.2 Setup Plots Directory
    plot_diagnostics = args.plot_diagnostics
    plots_dir = None
    if plot_diagnostics:
        plots_dir = os.path.join(args.outdir, "plots")
        if not os.path.exists(plots_dir):
            os.makedirs(plots_dir)
        logging.info(f"Diagnostic plotting enabled. Plots will be saved to: {plots_dir}")
    
    # --- 1. Load Data ---
    logging.info(f"Loading metadata from {args.metadata}...")
    try:
        # forcing "sample" column to be string to avoid int/string mismatch with count matrix columns
        meta_df = pd.read_csv(args.metadata, sep='\t', dtype={'sample': str})
        validate_metadata(meta_df)
    except Exception as e:
        logging.error(f"Metadata validation failed.\n{e}")
        sys.exit(1)

    logging.info(f"Loading counts from {args.count_matrix}...")
    try:
        counts_df = pd.read_csv(args.count_matrix, sep='\t', index_col=0)
    except Exception as e:
        logging.error(f"Could not read count matrix file. Details: {e}")
        sys.exit(1)

    # --- NEW: Check for Duplicate Gene IDs in the Index ---
    if counts_df.index.duplicated().any():
        num_dups = counts_df.index.duplicated().sum()
        # Get up to 5 examples to show the user
        dup_examples = counts_df.index[counts_df.index.duplicated()].unique().tolist()[:5]
        logging.error(f"ERROR: The first column of the count matrix (Gene IDs) contains {num_dups} duplicate entries. "
                 f"Each row must have a unique identifier. "
                 f"Examples of duplicates found: {dup_examples}...")
        sys.exit(1)
    
    # --- 2. Data Integrity Check ---
    # Ensure all samples in metadata exist in count matrix
    meta_samples = set(meta_df['sample'])
    count_samples = set(counts_df.columns)
    
    missing_samples = meta_samples - count_samples
    if missing_samples:
        logging.error(f"ERROR: The following samples are in metadata but missing from count matrix:\n{missing_samples}")
        sys.exit(1)
    
    # --- 3. library size normalization ---
    logging.info("Normalizing counts (DESeq2 method)...")
    
    try:
        # convert to list for further processing
        meta_samples = list(meta_samples)
        count_samples = list(count_samples)
        norm_counts_df, sfs_df = processing.deseq2_normalize(counts_df, meta_samples)
        log2_norm_counts_df = np.log2(norm_counts_df + 1)
        
        # Save the size factor stats for the user to see
        sf_path = os.path.join(args.outdir, "size_factors.tsv")
        sfs_df.to_csv(sf_path, sep='\t', index=False)
        spearman_corr = stats.spearmanr(a=sfs_df['sf'],b=sfs_df['read_sum'])[0]
        logging.info(f"Size factors saved to {sf_path}, sf correlation with read sum is {str(np.round(spearman_corr,2))}")
        if plot_diagnostics:
            reporting.plot_size_factors(sfs_df, meta_df, plots_dir)
    except Exception as e:
        logging.error(f"Normalization failed: {e}")
        sys.exit(1)

    # --- 4.1 Fit size factor distribution ---
    logging.info("fitting size factor distribution...")
    try:
        sf_params = fitting.fit_size_factors(sfs_df['sf']) # list, containing s, loc, scale params

        sf_params_txt = [str(np.round(elem,2)) for elem in sf_params]
        logging.info("estimated size factor params, s, loc, scale: "+','.join(sf_params_txt))
        
        if plot_diagnostics:
            reporting.plot_simul_vs_observed_size_factors(sfs_df, sf_params, plots_dir)
    except Exception as e:
        logging.error(f"fitting size factor distribution failed: {e}")
        sys.exit(1)          
    
    # --- 4.2 Fit expression parameters ---    
    logging.info("fitting Total-vs-Fraction expression parameteres independently in conditions...")
    try:
        conditions = list(meta_df['EXP_CTL'].unique())
        for condition in conditions:
            condition_meta_df = meta_df.loc[meta_df['EXP_CTL'] == condition]
            
            samples_total = condition_meta_df[condition_meta_df['fraction'] == 'T']['sample'].values
            samples_frac  = condition_meta_df[condition_meta_df['fraction'] == 'F']['sample'].values
    
            if len(samples_total) == 0 or len(samples_frac) == 0:
                logging.error(f"ERROR: Could not find both 'T' and 'F' fractions for the {condition} condition.")
                sys.exit(1)
            
            # Calculate mean observed expression for the condition
            obs_total = log2_norm_counts_df[samples_total].mean(axis=1).values
            obs_frac  = log2_norm_counts_df[samples_frac].mean(axis=1).values
            
            logging.info(f"Log2 Total expression observed mean and std for {condition} condition are {str(np.round(np.mean(obs_total),2))} and {str(np.round(np.std(obs_total),2))}")
            logging.info(f"fitting Total expression parameteres for {condition} condition using Student's T distibution")
            fitter = fitting.Fitter(obs_total=obs_total, obs_frac=obs_frac, sf_params = sf_params, n_replicates=len(samples_total))
            best_params = fitter.optimize_total()
            best_params_txt = [str(np.round(elem,2)) for elem in best_params]
            logging.info(f"Expression parameters for Total fraction of {condition} successfully fit, df, loc, scale for t distr of log2 molecule abundance: "+','.join(best_params_txt))
            
    except Exception as e:
        logging.error(f"fitting Total-vs-Fraction expression parameteres failed: {e}")
        sys.exit(1)    
    # # --- 4. Fit ---
    # fitter = Fitter(obs_total=obs_total, obs_frac=obs_frac, n_replicates=len(samples_total))
    
    # # Pass the specific size factors corresponding to these samples
    # relevant_sfs = size_factors[np.concatenate([samples_total, samples_frac])]
    # fitter.fit_size_factors(relevant_sfs.values)
    
    # best_params = fitter.optimize_total()
    
    # # --- 5. Final Simulation & Report ---
    # sim_df = simulate_log2_nCPM_from_one_condition(
    #     best_params, fitter.sf_params, N_g=len(obs_total)
    # )
    
    # generate_summary(args.outdir, obs_total, sim_df['log2_nCPM_Total'], best_params)
    logging.info(f"Success! Output written to {args.outdir}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        # Standard exit code for Ctrl+C is 130
        sys.stderr.write("\nProcess interrupted by user. Exiting.\n")
        sys.exit(130)