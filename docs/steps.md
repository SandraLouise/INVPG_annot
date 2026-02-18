# Running INVPG-annot step by step

## 1. Selecting the bubbles to process

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

## 2. Annotating the selected bubbles

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