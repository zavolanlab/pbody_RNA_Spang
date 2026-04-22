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
├── pbody_RNA_Spang.template.env        # Template for required environment variables/paths
├── upload_to_hf.no_token.py            # python script which used HF API to upload selected data to Hugging Face
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-basic                 # Snakemake pipeline for RNA-seq and cCLAP data processing optimized for SMART-Seq Total RNA Library Prep with ZapR Depletion (with UMIs)
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

### 3. Install the conda environment with zavolab_pyutils
Analysis in the notebook is largely based on the functions from [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/dev) repository.
Follow the instruction from that repo "Developer Setup from source, with conda environment".
Use the created conda environment "zavolab_pyutils" to execute the Jupyter Notebook.

### 4. Essential for developpers! Install nbstripout
When in the `pbody_RNA_Spang` directory, run:
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