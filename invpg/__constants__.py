# Software description
SOFT_DESCRIPTION: str = "A tool to annotate inversions from pangenome graph bubbles."
# Global commands help
HELP_COMMAND_INVANNOT: str = "Annotates the bubbles as \"INV:path\" or \"INV:aln\"."
HELP_COMMAND_FILTERVCF: str = "Selects bubbles corresponding to putative balanced SVs."
HELP_COMMAND_RESCUEINV: str = "Detects one-node inversions that may be missing from vg deconstruct VCF."
HELP_COMMAND_FILTANNOT: str = "<help string for filter_annot.py"
# Input help strings
HELP_INPUT_FILE_GFA: str = "Path to a GFA-like file.  Should be provided solely when not using minigraph graphs."
HELP_INPUT_FILE_VCF: str = "Path to a VCF file."
HELP_INPUT_FILE_BED: str = "Path to a BED file. Should be provided solely when working with minigraph graphs."
# Parameters help strings
HELP_PARAM_THREADS: str = "Number of threads used for parallelization (minimap2)."
HELP_PARAM_PERCENTAGE: str = "Estimated percentage of genome divergence (for variants filtering)."
HELP_PARAM_MINCOV: str = "Minimum coverage of inversion signal."
HELP_PARAM_REFID: str = "ID for reference to use in output."
HELP_PARAM_KEEP: str = "Keep temporary files after pipeline completion (mostly for debugging purposes)."
# Default value for parameters
DEFAULT_PERCENTAGE: int = 1
DEFAULT_MINCOV: float = .5
DEFAULT_THREADS: int = 1
