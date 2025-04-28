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
HELP_PARAM_PERCENTAGE: str = "This parameter controls the leniency of the algorithm towards allele size difference (in nt) in the first step of variant/bubble filtering. Only the non-reference alleles that have a size difference <= (d * max allele size / 100) will go through the annotation step. (default: 10)"
HELP_PARAM_MINCOV: str = "Minimum coverage of inversion signal as fraction of bubble length. (default: 0.5)"
HELP_PARAM_REFID: str = "ID for reference to use in output."
HELP_PARAM_KEEP: str = "Keep temporary files after pipeline completion (mostly for debugging purposes)."
HELP_PARAM_OUTPUT_MAIN: str = "Name of output BED file (without file extension). Can be a path to control output directory."
HELP_PARAM_OUTPUT: str = "Output path."
# Default value for parameters
DEFAULT_PERCENTAGE: int = 10
DEFAULT_MINCOV: float = .5
DEFAULT_THREADS: int = 1
