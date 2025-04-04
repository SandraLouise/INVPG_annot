#!/usr/bin/env python3
from invpg.__constants__ import *
from pathlib import Path
from argparse import ArgumentParser
from sys import argv
from os import remove, listdir
from shutil import rmtree
from invpg.inv_annot import invannot
from invpg.variant_filter import filter_main
from invpg.rescue_1node_inv import search_bed
from invpg.filter_annot import filterannot_main

parser: ArgumentParser = ArgumentParser(
    description=SOFT_DESCRIPTION,
    add_help=True
)
subparsers = parser.add_subparsers(
    help='Available subcommands',
    dest="subcommands"
)
parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'


#####################
## GENERAL PARSER ##
####################

parser.add_argument(
    "-v",
    "--input_vcf_file",
    type=str,
    help=HELP_INPUT_FILE_VCF,
)
parser.add_argument(
    "-g",
    "--input_gfa_file",
    type=str,
    help=HELP_INPUT_FILE_GFA,
)
parser.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    help=HELP_PARAM_OUTPUT,
)
parser.add_argument(
    "-d",
    "--div_percentage",
    type=int,
    help=HELP_PARAM_PERCENTAGE,
    default=DEFAULT_PERCENTAGE,
)
parser.add_argument(
    "-k",
    "--keep_files",
    help=HELP_PARAM_KEEP,
    action='store_true',
    default=False,
)
parser.add_argument(
    "-t",
    "--threads",
    type=int,
    help=HELP_PARAM_THREADS,
    default=DEFAULT_THREADS,
)
parser.add_argument(
    "-m",
    "--mincov",
    type=float,
    help=HELP_PARAM_MINCOV,
    default=DEFAULT_MINCOV,
)

########################
## INDIVIDUAL PARSERS ##
########################

## Subparser for invannot ##

parser_invannot: ArgumentParser = subparsers.add_parser(
    'annot',
    help=HELP_COMMAND_INVANNOT,
)
parser_invannot.add_argument(
    "input_vcf_file",
    type=str,
    help=HELP_INPUT_FILE_VCF,
)
parser_invannot.add_argument(
    "input_gfa_file",
    type=str,
    help=HELP_INPUT_FILE_GFA,
)
parser_invannot.add_argument(
    "-t",
    "--threads",
    type=int,
    help=HELP_PARAM_THREADS,
    default=DEFAULT_THREADS,
)
parser_invannot.add_argument(
    "-m",
    "--mincov",
    type=float,
    help=HELP_PARAM_MINCOV,
    default=DEFAULT_MINCOV,
)

## Subparser for filtervcf ##

parser_filtervcf: ArgumentParser = subparsers.add_parser(
    'filtvcf',
    help=HELP_COMMAND_FILTERVCF,
)
parser_filtervcf.add_argument(
    "input_vcf_file",
    type=str,
    help=HELP_INPUT_FILE_VCF,
)
parser_filtervcf.add_argument(
    "-d",
    "--div_percentage",
    type=int,
    help=HELP_PARAM_PERCENTAGE,
    default=DEFAULT_PERCENTAGE,
)

parser_rescueinv: ArgumentParser = subparsers.add_parser(
    'rescue',
    help=HELP_COMMAND_RESCUEINV,
)
parser_rescueinv.add_argument(
    "-g",
    "--input_gfa_file",
    type=str,
    default=None,
    help=HELP_INPUT_FILE_GFA,
)


parser_filtannot: ArgumentParser = subparsers.add_parser(
    'filtannot',
    help=HELP_COMMAND_FILTANNOT,
)
parser_filtannot.add_argument(
    "-m",
    "--mincov",
    type=float,
    help=HELP_PARAM_MINCOV,
    default=DEFAULT_MINCOV,
)

#######################################
args = parser.parse_args()
#######################################


def main() -> None:
    "Main call for subprograms"
    if len(argv) == 1:
        print(
            "You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit(1)
    match args.subcommands:
        case 'annot':
            invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=args.input_vcf_file,
                temp_folder=f"tmp_{Path(args.gfa_file).stem}/",
                mincov=args.mincov,
                threads=args.threads,
            )
        case 'filtvcf':
            filter_main(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
            )
        case 'rescue':
            search_bed(
                in_file=args.input_gfa_file,
            )
        case 'filtannot':
            filterannot_main(
                input_annotation_file=args.input_bed_file,
                reference_name=args.reference_path,
                minimum_coverage=args.mincov,
            )
        case _:
            # First we filter the VCF file
            temp_output_vcf: str = filter_main(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
            )
            # Then we rescue nodes in inversions that weren't described in the VCF
            bed_file_raw = invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=temp_output_vcf,
                temp_folder=(
                    temp_folder := f"tmp_{Path(args.input_gfa_file).stem}/"
                ),
                mincov=args.mincov,
                threads=args.threads,
            )
            if not args.keep_files:
                remove(temp_output_vcf)
                for file in listdir(temp_folder):
                    remove(f"{temp_folder}{file}")
                rmtree(temp_folder)

    exit(0)
