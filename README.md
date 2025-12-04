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
> Please see [Commands](docs/pipeline.md) to learn more about how to prepare your files for invpg-annot.

You can use a single command to execute the whole pipeline or do a step-by-step analysis (see [Steps](docs/steps.md)).

```bash
usage: invpg [-h] [-v INPUT_VCF_FILE] [-g INPUT_GFA_FILE] [-o OUTPUT_PREFIX] [-d DIV_PERCENTAGE] [-m MINCOV] [-k] [-t THREADS] [-O OUTPUT_VCF_FILE]

A tool to annotate inversions from pangenome graph bubbles.
  -h, --help            show this help message and exit
  -v, --input_vcf_file INPUT_VCF_FILE
                        Path to a VCF file.
  -g, --input_gfa_file INPUT_GFA_FILE
                        Path to a GFA-like file. Should be provided solely when not using minigraph graphs.
  -o, --output_prefix OUTPUT_PREFIX
                        Name/path of output VCF file. If parent folder of output VCF file doesn't already exist, it will be created.
  -d, --div_percentage DIV_PERCENTAGE
                        This parameter controls the leniency of the algorithm towards allele size difference (in nt) in the first step of variant/bubble
                        filtering. Only the non-reference alleles that have a size difference <= (d * max allele size / 100) will go through the annotation
                        step. (default: 10)
  -m, --mincov MINCOV   Minimum coverage of inversion signal as fraction of bubble length. (default: 0.5)
  -k, --keep_files      Keep temporary files after pipeline completion (mostly for debugging purposes).
  -t, --threads THREADS
                        Number of threads used for parallelization (minimap2).
```

### Test with a small dataset

To check that INVPG-annot behaves as expected on your device, you can run:

```bash
cd test-dir/
invpg -v test_bubbles.vcf -g test_graph.gfa -o test_annotation.bed -m 0.5 -d 10
diff expected_annotation.bed test_annotation.bed
```

To explore the intermediate output files (described [here](https://github.com/SandraLouise/INVPG_annot?tab=readme-ov-file#intermediate-files-when-using--k-parameter)) on a small dataset, run:

```bash
mkdir outputfiles
cd outputfiles
invpg -v ../test_bubbles.vcf -g ../test_graph.gfa -o test_annotation.bed -k -m 0.5 -d 10
cd res_*
```

### Impact of the `-d` and `-m` parameters

The `-d` parameter controls the leniency of the algorithm towards allele size difference (in nucleotides) in the first step of variant/bubble filtering. Only the bubbles that pass this filtering step will be processed in the annotation step. The main motive behind the filtering step is to speed up the annotation step by ignoring the bubbles that are considered unlikely to represent inversions due to their allele sizes (e.g. SNPs, indels, deletions, insertions), given that inversion bubbles are expected to represent only a small portion of the total bubbles in a pangenome graph. As the sizes of two inversion alleles may not be strictly identical in a pangenome graph due to biological factors (inner genomic variation) and/or artificial factors (alignment artefacts generated during pangenome graph inference), the algorithm uses the `-d` parameter to define the level of allele size difference acceptable in a potential inversion bubble as `len_Am * d / 100` (`len_Am` being the size of the largest allele of the bubble). With `-d 0`, only the variants/bubbles that have at least two alleles with identical size will go through the annotation step. The higher the `-d` value used, the less likely inversion bubbles will be wrongly discarded due to high size difference between alleles, but the more time the annotation step will take. Based on several tests, we advise to use `-d 10` (allowing for an allele size difference of 10% of the largest allele) even with low divergence between genomes.


The `-m` parameter sets the minimum coverage of inversion signal (as a fraction of the bubble nucleotidic length) that must be found on a bubble for it to be reported in the output. We advise to set is as `-m 0.5`. For more details on how this coverage is calculated, please see the Method section of our [paper](https://doi.org/10.1101/2025.03.14.643331).

### Output files

#### Final BED annotation file

By default, the `invpg` command outputs a single BED file in which each line describes a bubble annotated as inversion. The three first fields contain classical BED information (*i.e.* reference chromosome ID, start position, end position). The fourth field contains additionnal information about the annotation for each non-reference allele that passed the first step filter (separated by `;`) in the form "INV:`signaltype`:`cov`,`x`" or "DIV" (when the inversion signal was < `--mincov`). 

- `signaltype` can be either "path" or "aln" depending on the source of the signal used for the annotation (path-based or alignment-based).
- `cov` is the fraction of inversion signal coverage.
- `x` corresponds to additionnal statistics that can be ignored (used for development that will soon be removed from the output).

#### Intermediate files (when using `-k` parameter)

There are two types of intermediate files: one VCF file and multiple PAF files.

- The VCF file (name ending with `.balancedSV.vcf`) contains all input VCF lines that pass the first step filter. It is the file used as input for the second step of annotation.
- The PAF files contain the minimap2 alignment results for each allele going through the second step of annotation. The number of PAF files can be very large depending on the input VCF contents. We plan to optimize the number of PAF files generated in the future.

## Citation

Romain, S., Dubois, S., Legeai, F., & Lemaitre, C. (2025). Investigating the topological motifs of inversions in pangenome graphs. *bioRxiv*, 2025-03, https://doi.org/10.1101/2025.03.14.643331.

## Contact

INVPG-annot is a [Genscale](https://team.inria.fr/genscale/) tool developed by Sandra Romain, Siegfried Dubois, Fabrice Legeai and Claire Lemaitre. For any bug report or feedback, please use the Github Issues form.
