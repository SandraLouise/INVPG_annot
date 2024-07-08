# INVPG_annot
A tool to annotate inversions from pangenome graph bubbles.

## Usage

### 0. Prerequisites

- `GFA`  A pangenome graph in GFA format.
- `VCF`  The bubbles extracted from the `GFA` in VCF format using `vg deconstruct`.
- minimap2

### 1. Selecting the bubbles to process

Selects bubbles corresponding to putative balanced SVs.

```bash
python variant_filter.py <VCF> <outPrefix> <divPct>
```

- `outPrefix`  Prefix of the output file.
- `divPct`  Estimated percentage of divergence of the genomes in the pangenome graph. Defines the leniency to consider a variant as balanced.

Output: 
- `outPrefix`.balancedSV.vcf  A VCF file with selected bubbles.

### 2. Annotating the selected bubbles

> **Requires minimap2.**

Annotates the bubbles as "INV:path" or "INV:aln".

```bash
python inv_annot.py <balancedVCF> <GFA> <outPrefix> <threads>
```

- `outPrefix`  Prefix of the output files.
- `threads`  Number of threads to use for the sequence alignment (minimap2).

Output:
- `outPrefix`.annot.tsv  A TSV (tabular separated) file with INV annotated bubbles, one bubble per line.

### 3. Detecting one-node inversions from the graph

Detects one-node inversions that may be missing from `vg deconstruct` VCF.

```bash
python rescue_1node_inv.py <GFA> gfa <refPath>
```

- `refPath`  ID of the reference path in the graph (same one as used by `vg deconstruct`).
