# RnaFracMod

[![license][badge-license]][badge-url-license]

RnaFracMod is a command-line tool for modelling count data from bulk RNA-seq-like experiments involving total fraction and a sub-fraction, for example RNAseq-RiboSeq or RNAseq from pull-down or subcellular fraction experiments

## Installation

### Using conda (recommended)

Create the conda environment with necessary dependencies:
```bash
conda env create -f install/RnaFracMod.yaml
conda activate rnafracmod
```

Inside created conda environment, run:
```bash
pip install -e .
```
## Usage

```bash
rnafracmod [-h] --count_matrix COUNT_MATRIX --metadata METADATA --outdir OUTDIR

RNA Fractionation Modelling Tool

options:
  -h, --help            show this help message and exit
  --count_matrix COUNT_MATRIX
                        Path to TSV count matrix
  --metadata METADATA   Path to sample metadata TSV
  --outdir OUTDIR       Output directory
  --plot_diagnostics    If set, generates diagnostic plots in <outdir>/plots/
```

### Input Data Format

To run `rnafracmod`, you need two input files: a count matrix and a metadata file.

#### 1. Count Matrix (`--count_matrix`)
A standard tab-separated (`.tsv`) file where rows are genes and columns are samples.
* **Rows:** Gene IDs / Transcript IDs or similar, should not have duplicated entries
* **Columns:** Sample IDs (must match metadata)
* **Values:** Raw read counts (integers)

#### 2. Metadata File (`--metadata`)
A tab-separated (`.tsv`) file describing your experiment. **Strict formatting is required.** It must contain at least these three columns:

| Column | Description | Allowed Values |
| :--- | :--- | :--- |
| `sample` | Unique identifier for the sample | Must match columns in Count Matrix |
| `EXP_CTL` | Experimental Condition | `CTL` (Control) or `EXP` (Experiment) |
| `fraction` | The type of lysate | `T` (Total) or `F` (Sub-fraction/Pull-down/Ribo-Seq etc) |

**Example `metadata.tsv`:**
```tsv
sample	EXP_CTL	fraction
Sample_01	CTL	T
Sample_02	CTL	F
Sample_03	EXP	T
Sample_04	EXP	F
```
**Note** current fitting algorithm does not account for **pairing** (e.g., that Sample_01 Total comes from the same cell culture as Sample_02 Fraction).

## License

GNU General Public License v2.0

[badge-license]: <https://img.shields.io/badge/License-GPL_v2-blue.svg>
[badge-url-license]: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>