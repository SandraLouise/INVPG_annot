# INVPG_annot
A tool to annotate inversions from pangenome graph bubbles.

## Installation

> [!WARNING]\
> **Prerequisites:**
> - minimap2
> - python $\geq$ 3.10

```bash
git clone https://github.com/SandraLouise/INVPG_annot.git
cd INVPG_annot
pip install -r requirements.txt --upgrade
python -m pip install . --quiet
```

## Usage

> [!WARNING]\
> **Prerequisites:**
> - `GFA` - A pangenome graph in GFA format.
> - `VCF` - The bubbles extracted from the `.gfa` in VCF format. Note that INVPG-annot has only been tested on and developed based on the formats of VCFs produced by `vg deconstruct` or the [`minigraph-call` pipeline](https://github.com/lh3/minigraph?tab=readme-ov-file#sv-calling-showcase-human-mhc).

You can use a single command to execute the whole pipeline or do a step-by-step analysis (see below).

```bash
usage: invpg [-h] -v INPUT_VCF_FILE -g INPUT_GFA_FILE [-d DIV_PERCENTAGE] [-k] [-t THREADS] [-m MINCOV]

A tool to annotate inversions from pangenome graph bubbles.
  -h, --help            show this help message and exit
  -v INPUT_VCF_FILE, --input_vcf_file INPUT_VCF_FILE
                        Path to a VCF file.
  -g INPUT_GFA_FILE, --input_gfa_file INPUT_GFA_FILE
                        Path to a GFA-like file.
  -d DIV_PERCENTAGE, --div_percentage DIV_PERCENTAGE
                        Originally intended as the estimated percentage of genome divergence. 
                        This parameter controls the leniency of the algorithm towards allele 
                        size difference (in nt) in the first step of variant/bubble filtering. 
                        Now advised to be set as `-d 10` regardless of genome divergence level.  
  -k, --keep_files      Keep temporary files after pipeline completion (mostly for debugging purposes).
  -r REFERENCE_PATH, --reference_path REFERENCE_PATH
                        ID for reference to use in output.
  -t THREADS, --threads THREADS
                        Number of threads used for parallelization (minimap2).
  -m MINCOV, --mincov MINCOV
                        Minimum coverage of inversion signal. Advised to be set at 0.5.
```

Example command line:

```bash
invpg -v bubbles.vcf -g graph.gfa -d 10 -t 8 -m 0.5
```

### Impact of the `-d` and `-m` parameters

The `-d` parameter controls the leniency of the algorithm towards allele size difference (in nucleotides) in the first step of variant/bubble filtering. Only the bubbles that pass this filtering step will be processed in the annotation step. The main motive behind the filtering step is to speed up the annotation step by ignoring the bubbles that are considered unlikely to represent inversions due to their allele sizes (e.g. SNPs, indels, deletions, insertions), given that inversion bubbles are expected to represent only a small portion of the total bubbles in a pangenome graph. As the sizes of two inversion alleles may not be strictly identical in a pangenome graph due to biological factors (inner genomic variation) and/or artificial factors (alignment artefacts generated during pangenome graph inference), the algorithm uses the `-d` parameter to define the level of allele size difference acceptable in a potential inversion bubble as `len_Am * d / 100` (`len_Am` being the size of the largest allele of the bubble). With `-d 0`, only the variants/bubbles that have at least two alleles with identical size will go through the annotation step. The higher the `-d` value used, the less likely inversion bubbles will be wrongly discarded due to high size difference between alleles, but the longer the annotation step will take. Based on several tests, we advise to use `-d 10` (allowing for an allele size difference of 10% of the largest allele) even with low divergence between genomes.


The `-m` parameter sets the minimum coverage of inversion signal (as a fraction of the bubble nucleotidic length) that must be found on a bubble for it to be reported in the output. We advise to set is as `-m 0.5`. For more details on how this coverage is calculated, please see the Method section of our [paper](https://doi.org/10.1101/2025.03.14.643331).

## Running INVPG-annot step by step

### 1. Selecting the bubbles to process

Selects bubbles corresponding to putative balanced SVs.

```bash
usage: invpg filtvcf [-h] [-d DIV_PERCENTAGE] input_vcf_file

positional arguments:
  input_vcf_file        Path to a VCF file.

options:
  -h, --help            show this help message and exit
  -d DIV_PERCENTAGE, --div_percentage DIV_PERCENTAGE
                        Estimated percentage of genome divergence (for variants filtering).
```

- `input_vcf_file`  Unfiltered VCF.
- `divPct`  Estimated percentage of divergence of the genomes in the pangenome graph. Defines the leniency to consider a variant as balanced.

Output: 
- `input_vcf_file`.balancedSV.vcf  A VCF file with selected bubbles.

### 2. Annotating the selected bubbles

> [!WARNING]\
> **Requires minimap2.**

Annotates the bubbles as "INV:path" or "INV:aln".

```bash
usage: invpg annot [-h] [-t THREADS] [-m MINCOV] input_vcf_file input_gfa_file

positional arguments:
  input_vcf_file        Path to a VCF file.
  input_gfa_file        Path to a GFA-like file.

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads used for parallelization (minimap2).
  -m MINCOV, --mincov MINCOV
                        Minimum coverage of inversion signal.
```

- `input_vcf_file`  Filtered VCF file.
- `threads`  Number of threads to use for the sequence alignment (minimap2).
- `mincov` Minimum coverage of inversion signal.

Output:
- `input_vcf_file`.annot.tsv  A TSV (tabular separated) file with INV annotated bubbles, one bubble per line.

### 3. Detecting one-node inversions from the graph [DEPRECATED]

> [!WARNING]\
> This step is no longer part of the INVPG-annot default pipeline and is deprecated.

Detects one-node inversions that may be missing from `vg deconstruct` VCF.

```bash
usage: invpg rescue [-h] [-g INPUT_GFA_FILE] [-b INPUT_BED_FILE] [-r REFERENCE_PATH]

options:
  -h, --help            show this help message and exit
  -g INPUT_GFA_FILE, --input_gfa_file INPUT_GFA_FILE
                        Path to a GFA-like file. Should be provided solely when not using minigraph graphs.
  -b INPUT_BED_FILE, --input_bed_file INPUT_BED_FILE
                        Path to a BED file. Should be provided solely when working with minigraph graphs.
  -r REFERENCE_PATH, --reference_path REFERENCE_PATH
                        ID for reference to use in output.
```

### 4. Filtering annotations [DEPRECATED]

> [!WARNING]\
> This step is no longer part of the INVPG-annot default pipeline and is deprecated.

```bash
usage: invpg filtannot [-h] [-b INPUT_BED_FILE] [-r REFERENCE_PATH] [-m MINCOV]

options:
  -h, --help            show this help message and exit
  -b INPUT_BED_FILE, --input_bed_file INPUT_BED_FILE
                        Path to a BED file. Should be provided solely when working with minigraph graphs.
  -r REFERENCE_PATH, --reference_path REFERENCE_PATH
                        ID for reference to use in output.
  -m MINCOV, --mincov MINCOV
                        Minimum coverage of inversion signal.
```

## Citation

Romain, S., Dubois, S., Legeai, F., & Lemaitre, C. (2025). Investigating the topological motifs of inversions in pangenome graphs. *bioRxiv*, 2025-03, https://doi.org/10.1101/2025.03.14.643331.

## Contact

INVPG-annot is a [Genscale](https://team.inria.fr/genscale/) tool developed by Sandra Romain, Siegfried Dubois, Fabrice Legeai and Claire Lemaitre. For any bug report or feedback, please use the Github Issues form.