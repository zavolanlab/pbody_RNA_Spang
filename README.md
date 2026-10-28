# Analysis of P-body associated RNA in yeast (In collaboration with Anne Spang lab) - Analysis and Pipelines

This repository contains the computational workflows and downstream analysis notebooks related to the analysis of P-body associated RNA, based on RNA-seq and cCLAP experiments conduct in Saccharomyces cerevisiae.

The repository is optimized for running BOTH the workflows and analysis in jupyter notebook on HPC cluster.

On sciCORE HPC, running jupyter notebook on a computational node is nicely enabled by [OnDemand service](https://docs.scicore.unibas.ch/HPC%20Cluster/interactivecomputing/#open-ondemand-ood).

We utilize a hybrid approach: **Snakemake** for robust, scalable data processing on HPC clusters (sciCORE), and **Jupyter Notebooks** for interactive downstream analysis and visualization.

# Current state
Currently, we're analyzing the bulk RNA-seq and cCLAP data from yeast experiments.

## Repository Structure

```text
.
├── pbody_RNA_Spang.current.ipynb       # a Jupyter notebook dedicated to the project, includes analysis and workflow configuration
├── pbody_RNA_Spang.fracSanity.ipynb    # P-body pull-down vs Total with Sanity_utils: global fraction F, per-gene fractions, imaging ground truth
├── pbody_RNA_Spang.TotalDE.ipynb       # differential expression in the Total fraction, MA plots, gene-set analysis (Sanity_utils)
├── environment.yml                     # conda environment of the project (Sanity_utils, zavolab_pyutils and notebook packages)
├── colab/                              # Google Colab notebook for collaborators and the SWITCHdrive sync set-up
├── pbody_RNA_Spang.template.env        # Template for required environment variables/paths
├── upload_to_hf.no_token.py            # python script which used HF API to upload selected data to Hugging Face
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-basic                 # Snakemake pipeline for RNA-seq and cCLAP data processing optimized for SMART-Seq Total RNA Library Prep with ZapR Depletion (with UMIs)
    ├── Snakefile-fracSanity            # Sanity_utils pipeline on the featureCounts tables: Sanity, fracSanity, gene sets, sequence/structure features, motifs
    ├── config.fracSanity.template.yaml # Template configuration of Snakefile-fracSanity
    ├── config.template.yaml            # Template configuration for Snakemake parameters
    ├── envs/                           # Conda environments isolated for specific Snakemake rules
    ├── profile/                        # SLURM execution profile for the HPC
    └── scripts/                        # Python and R scripts utilized by both Snakemake and Jupyter
```

## Data from external databases/resources
!Attention: clicking on links below will automatically start downloading of big files.

**yeast genome and annotation**

[genome assembly R64-1-1 from Ensembl release 109 .fasta](https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz)

[corresponding genome annotation, .gtf](https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz)

## Quick Start & Setup

To ensure strict reproducibility and security, this project uses `.env` files to manage all absolute paths (data directories, genome annotations, etc.). **Do not hardcode paths into the Python or Snakemake files.**

### 1. Clone the Repository
Clone this repository into your local user space (`$HOME`):
```bash
git clone https://github.com/zavolanlab/pbody_RNA_Spang.git
cd pbody_RNA_Spang
```

### 2. Configure Environment Paths
You must map the project to your local HPC paths. 
**First**, copy the template, rename it, and fill in your absolute paths, for example like that:
  ```bash
  cp pbody_RNA_Spang.template.env pbody_RNA_Spang.scicore.env
  # Open .env and edit the "Base Directories" section to match your system
  ```
* **Recommended if you are a Zavolan group member on sciCORE:** move the `pbody_RNA_Spang.scicore.env` to Project GROUP folder and symlink into your local repository directory from step 1 (`pbody_RNA_Spang`):
  ```bash
  ln -s <a file with specified sciCORE paths> pbody_RNA_Spang.scicore.env
  ```
This way `pbody_RNA_Spang.scicore.env` will be automatically accessible by group members but will not be tracked by git.
*(Note: `*.env` files are ignored by git to protect private cluster paths, except the `pbody_RNA_Spang.template.env` file).*
**(`pbody_RNA_Spang.scicore.env` does exist in the GROUP folder of the Project on Scicore. Look for README there.)

### 3. Install the conda environment of the project
The analysis notebooks and the Sanity_utils workflow (`WF/Snakefile-fracSanity`) use one conda environment, defined in [environment.yml](environment.yml). It installs [Sanity_utils](https://github.com/nimwegenLab/Sanity_utils) (Sanity normalization, condition averages, differential expression, fracSanity, gene-set analysis; private repository, pip needs an SSH key with access to `nimwegenLab`) and [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/separated_functionalities) (general-purpose plots such as `pca_plot`) from GitHub, plus the notebook packages (e.g. `adjustText`, `matplotlib-venn`) and the tools used by the workflow scripts (Biopython, ViennaRNA's `RNAfold`).

```bash
cd pbody_RNA_Spang
PYTHONNOUSERSITE=1 conda env create -f environment.yml
conda activate pbody_rna_spang
# build the C++ Sanity binary into the env; a portable CPU flag, since SLURM jobs run on nodes with different CPUs
sanity-utils-install-sanity --arch-flag=-march=x86-64-v2
# Jupyter kernel (e.g. for sciCORE OnDemand); it sets PYTHONNOUSERSITE so that packages in ~/.local do not shadow the env's own
python -m ipykernel install --user --name pbody_rna_spang --display-name "Python (pbody_rna_spang)" --env PYTHONNOUSERSITE 1
```

Select the kernel "Python (pbody_rna_spang)" in Jupyter. For the workflow, set `sanity_utils_python` and `sanity_utils_bin` in the config (done by the notebooks) to this env's `bin/python` and `bin/`. Snakemake itself runs from a separate env with the SLURM executor plugin (`snakemake`, `snakemake-executor-plugin-slurm`). FIMO (motif scan) is taken from the sciCORE module `MEME/5.5.7-gompi-2023b` (config key `fimo_bin`). To update Sanity_utils or zavolab_pyutils later, re-run the corresponding `pip install` line of `environment.yml` with `--upgrade --no-deps`, and restart the Jupyter kernel.

### 4. Essential for developpers! Install nbstripout
When in the `pbody_RNA_Spang` directory, run (nbstripout is part of the conda environment):
```bash
nbstripout --install
```
This will automatically hide the output of cells in juputer notebooks when pushed to github! Otherwise there is a risk of exposing your HPC cluster paths to public.

### 5. Use the juputer notebook to configure the workflow and input table preparation
Configuration of the workflows (i.e. creation of input .tsv with sample specification and .yaml config is done **inside** the [jupyter notebook](pbody_RNA_Spang.current.ipynb))

### 6. Executing the Workflows

The heavy lifting is divided into (currently, three) separate Snakemake/Nextflow workflows.

`Bash` commands are also prepared inside the jupyter notebook. They should be further copied into command line and executed.

**On an HPC cluster like sciCORE**, workflows should in general be executed on a **login** node. Snakemake/Nextflow further automatically submits jobs to computational nodes.

### 7. Downstream Analysis
Once the workflows are complete, all results are routed to the shared group directories defined in your `.env` file.

Use respective sections of the Jupyter Notebook to analyze the outputs. 

The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate all workflow results, figures, and metadata regardless of where you cloned this repository.
### 8. HTML reports for collaborators
The analysis notebooks (`pbody_RNA_Spang.TotalDE.ipynb`, `pbody_RNA_Spang.fracSanity.ipynb`, `pbody_RNA_Spang.translation.ipynb`) double as reports. Each saved figure is shown with its file name and links to the full-size PNG/PDF. Cells tagged `report-skip` (setup, configuration, raw table dumps) are left out of the report.

After running the notebooks (Kernel > Restart & Run All), export them without code. Without notebook arguments, all three are exported; list notebooks to export only those:

```bash
conda activate sanity_utils
python WF/scripts/make_html_reports.py --out <report_folder>        # all three notebooks
python WF/scripts/make_html_reports.py --out <report_folder> pbody_RNA_Spang.TotalDE.ipynb   # one notebook only
```

The report folder contains one HTML file per exported notebook, an `index.html`, and `figures/` with the full-size files that the reports link to (only the figures linked from the exported notebooks are copied). Share it as a whole (zip, SWITCHdrive). Exports so far, in `$SHARED_PROJECT_DIR/reports/`:

| Folder | Notebooks | Content |
|---|---|---|
| `2026-10-28_TotalDE` | `pbody_RNA_Spang.TotalDE.ipynb` | Total fraction: normalization, PCA, differential expression, gene-set activities (literature sets and top-$z_m$ sets per database, biotype covariates), TIN, selected genes; 36 figure files (11 MB) |
