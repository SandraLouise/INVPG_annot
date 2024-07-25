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
> - `GFA` - A pangenome graph in GFA format. Can also need a `.bed` file in the case you're using minigraph.
> - `VCF` - The bubbles extracted from the `.gfa` in VCF format using `vg deconstruct`.

You can use a single command to execute the whole pipeline or do a step-by-step analysis (see below).

```bash
usage: invpg [-h] -v INPUT_VCF_FILE [-g INPUT_GFA_FILE] [-b INPUT_BED_FILE] [-d DIV_PERCENTAGE] [-k] [-r REFERENCE_PATH] [-t THREADS] [-m MINCOV]

A tool to annotate inversions from pangenome graph bubbles.
  -h, --help            show this help message and exit
  -v INPUT_VCF_FILE, --input_vcf_file INPUT_VCF_FILE
                        Path to a VCF file.
  -g INPUT_GFA_FILE, --input_gfa_file INPUT_GFA_FILE
                        Path to a GFA-like file. Should be provided solely when not using minigraph graphs.
  -b INPUT_BED_FILE, --input_bed_file INPUT_BED_FILE
                        Path to a BED file. Should be provided solely when working with minigraph graphs.
  -d DIV_PERCENTAGE, --div_percentage DIV_PERCENTAGE
                        Estimated percentage of genome divergence (for variants filtering).
  -k, --keep_files      Keep temporary files after pipeline completion (mostly for debugging purposes).
  -r REFERENCE_PATH, --reference_path REFERENCE_PATH
                        ID for reference to use in output.
  -t THREADS, --threads THREADS
                        Number of threads used for parallelization (minimap2).
  -m MINCOV, --mincov MINCOV
                        Minimum coverage of inversion signal.
```

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
  input_gfa_file        Path to a GFA-like file. Should be provided solely when not using minigraph graphs.

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

### 3. Detecting one-node inversions from the graph

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

### 4. Filtering annotations

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
