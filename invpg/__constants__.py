# Software description
SOFT_DESCRIPTION: str = "A tool to annotate inversions from pangenome graph bubbles."
# Global commands help
HELP_COMMAND_INVANNOT: str = "Annotates the bubbles as \"INV:path\" or \"INV:aln\"."
HELP_COMMAND_FILTERVCF: str = "Selects bubbles corresponding to putative balanced SVs."
HELP_COMMAND_RESCUEINV: str = "[DEPRECATED] Detects one-node inversions that may be missing from vg deconstruct VCF."
HELP_COMMAND_FILTANNOT: str = "[DEPRECATED] Filters annotations of inversions obtained from the annotation step."
# Input help strings
HELP_INPUT_FILE_GFA: str = "Path to a GFA-like file.  Should be provided solely when not using minigraph graphs."
HELP_INPUT_FILE_VCF: str = "Path to a VCF file."
HELP_INPUT_FILE_BED: str = "Path to a BED file, output of the main INVPG_annot pipeline."
# Parameters help strings
HELP_PARAM_THREADS: str = "Number of threads used for parallelization (minimap2)."
HELP_PARAM_PERCENTAGE: str = "Originally intended as the estimated percentage of genome divergence. This parameter controls the leniency of the algorithm towards allele size difference (in nt) in the first step of variant/bubble filtering. Now advised to be set as `-d 10` regardless of genome divergence level. "
HELP_PARAM_MINCOV: str = "Minimum coverage of inversion signal. Advised to be set at 0.5."
HELP_PARAM_REFID: str = "ID for reference to use in output."
HELP_PARAM_KEEP: str = "Keep temporary files after pipeline completion (mostly for debugging purposes)."
HELP_PARAM_OUTPUT_MAIN: str = "Output name, must be a path to a .bed file."
HELP_PARAM_OUTPUT: str = "Output path."
# Default value for parameters
DEFAULT_PERCENTAGE: int = 10
DEFAULT_MINCOV: float = .5
DEFAULT_THREADS: int = 1
