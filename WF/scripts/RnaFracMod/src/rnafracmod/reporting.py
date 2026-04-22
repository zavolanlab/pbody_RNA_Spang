import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import logging
import numpy as np
from adjustText import adjust_text
from scipy import stats

def generate_summary(outdir, obs_total, sim_total, params, sf_data=None):
    os.makedirs(outdir, exist_ok=True)
    
    # 1. Plot Total Distribution
    plt.figure(figsize=(6, 4))
    sns.histplot(obs_total, color='green', alpha=0.5, label='Observed', stat='density')
    sns.histplot(sim_total, color='orange', alpha=0.5, label='Simulated', stat='density')
    plt.legend()
    plt.title(f"Total Expression Fit\nParams: {np.round(params, 2)}")
    plot_path = os.path.join(outdir, 'total_fit.png')
    plt.savefig(plot_path)
    plt.close()

    # 2. Generate HTML
    html_content = f"""
    <html>
    <head><title>RnaFracMod Summary</title></head>
    <body style="font-family: sans-serif; padding: 20px;">
        <h1>Analysis Summary</h1>
        <h2>Fitted Parameters</h2>
        <p>Total Expression (df, loc, scale): {params}</p>
        
        <h2>Plots</h2>
        <div style="display: flex;">
            <img src="total_fit.png" width="500">
        </div>
    </body>
    </html>
    """
    
    with open(os.path.join(outdir, 'summary.html'), 'w') as f:
        f.write(html_content)

###
# Diagnostic plots
###

def plot_size_factors(sfs_df, meta_df, outdir):
    """
    Plots the distribution of Size Factors and Read Sums.
    """
    data = pd.merge(sfs_df.copy().reset_index(drop=True), 
                    meta_df.copy().reset_index(drop=True), 
                    how="left", on="sample")
    
    alpha_param, s_param = 0.7, 40
    sns.set(font_scale=1.2)
    sns.set_style("white")
    fig, axes = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(4, 4))
    
    x_feature, y_feature = "read_sum_mln", "sf"
    
    hue = "fraction"
    hue_order = ["F", "T"]
    palette = ["orange", "teal"]
    
    ax = sns.scatterplot(
        ax=axes,
        data=data,
        x=x_feature,
        y=y_feature,
        s=s_param,
        alpha=alpha_param,
        edgecolor="black",
        linewidth=0.5,
        hue=hue,
        hue_order=hue_order,
        palette=palette,
    )

    spearman_corr = stats.spearmanr(a=data[x_feature],b=data[y_feature])[0]
    
    ax.set(xlabel="# reads, mln",
           ylabel="size factor,\nDeseq2 normalization",
          title="spearman corr = "+str(np.round(spearman_corr,2)))
    ax.tick_params(left=True, bottom=True)
    
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc=2, borderaxespad=0.0, title="fraction", ncols=1)
    fig.savefig(
        os.path.join(outdir, "library_size_vs_SF.png"),
        bbox_inches="tight",
        dpi=600
    )
    fig.savefig(
        os.path.join(outdir, "library_size_vs_SF.pdf"),
        bbox_inches="tight",
        dpi=600,
    )

def plot_simul_vs_observed_size_factors(sfs_df, sf_params, outdir):
    s_mle, loc_mle, scale_mle = sf_params
    
    sf_simul = stats.lognorm.rvs(s=s_mle,loc=loc_mle, scale=scale_mle, size=1000)
    log2_sf_simul = np.log2(sf_simul)

    plot_feature = 'log2_sf'
    
    bins = 10
    bin_edges = np.linspace(np.min(sfs_df[plot_feature]), np.max(sfs_df[plot_feature]), bins + 1)
    
    sns.set(font_scale=0.7)
    sns.set_style("white")
    fig, axes = plt.subplots(1, 1, sharey=True, sharex=True, figsize=(2.5, 1.2))
    
    ax = sns.histplot(sfs_df[plot_feature],stat='proportion',bins=bin_edges,color='blue',label='observed')
    ax = sns.histplot(log2_sf_simul,stat='proportion',bins=bin_edges,color='orange',label='simulated')
    
    sf_params_txt = [str(np.round(elem,2)) for elem in sf_params]
    ax.set(title="$log_2$ size factors\nfitted params, s, loc, scale: "+','.join(sf_params_txt))
    ax.tick_params(left=True, bottom=True)
    
    ax.legend(
        bbox_to_anchor=(1.05, 1.00),
        loc="upper left",
        borderaxespad=0.0,
        title="",
        ncols=1,
    )
    fig.savefig(
        os.path.join(outdir, "simul_vs_observed_size_factors.png"),
        bbox_inches="tight",
        dpi=600
    )
    fig.savefig(
        os.path.join(outdir, "simul_vs_observed_size_factors.pdf"),
        bbox_inches="tight",
        dpi=600,
    )    