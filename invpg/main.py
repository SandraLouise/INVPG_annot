#!/usr/bin/env python3
from invpg.__constants__ import *
from pathlib import Path
from argparse import ArgumentParser
from sys import argv
from os import listdir, remove
from datetime import datetime
from invpg.inv_annot import invannot
from invpg.variant_filter import variant_filter
from invpg.rescue_1node_inv import search_bed
from invpg.filter_annot import filterannot

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
    help=HELP_PARAM_OUTPUT_MAIN,
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
    "-o",
    "--output_prefix",
    type=str,
    help=HELP_PARAM_OUTPUT,
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
    "-o",
    "--output_prefix",
    type=str,
    help=HELP_PARAM_OUTPUT,
)
parser_filtervcf.add_argument(
    "-d",
    "--div_percentage",
    type=int,
    help=HELP_PARAM_PERCENTAGE,
    default=DEFAULT_PERCENTAGE,
)

## Subparser for rescue (DEPRECATED) ##

parser_rescueinv: ArgumentParser = subparsers.add_parser(
    'rescue',
    description=HELP_COMMAND_RESCUEINV,
)
parser_rescueinv.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    help=HELP_PARAM_OUTPUT,
)
parser_rescueinv.add_argument(
    "-b",
    "--input_bed_file",
    type=str,
    help=HELP_INPUT_FILE_BED,
)

## Subparser for filtannot (DEPRECATED) ##

parser_filtannot: ArgumentParser = subparsers.add_parser(
    'filtannot',
    description=HELP_COMMAND_FILTANNOT,
)
parser_filtannot.add_argument(
    "-b",
    "--input_bed_file",
    type=str,
    help=HELP_INPUT_FILE_BED,
)
parser_filtannot.add_argument(
    "-r",
    "--reference_name",
    type=str,
    help=HELP_PARAM_REFID,
)
parser_filtannot.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    help=HELP_PARAM_OUTPUT,
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

    # This timestamp helps identify temporary files of this run
    ts = datetime.now()
    print(f"Starting job {str(ts)}")
    timestamp: str = str(ts).replace(' ', '_')
    match args.subcommands:
        case 'annot':
            invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=args.input_vcf_file,
                temp_folder=f"tmp_{Path(args.gfa_file).stem}/",
                mincov=args.mincov,
                threads=args.threads,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
        case 'filtvcf':
            filter(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
        case 'rescue':
            # DEPRECATED
            search_bed(
                in_file=args.input_gfa_file,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
        case 'filtannot':
            # DEPRECATED
            filterannot(
                input_annotation_file=args.input_bed_file,
                reference_name=args.reference_name,
                minimum_coverage=args.mincov,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
        case _:
            # First we filter the VCF file
            print("[" + str(datetime.now()) + "] STEP 1: filtering VCF file")
            temp_output_vcf: str = variant_filter(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
            print(f"Results output in {temp_output_vcf}")
            # Then we rescue nodes in inversions that weren't described in the VCF
            print("[" + str(datetime.now()) +
                  "] STEP 2: rescuing nodes in inversions")
            invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=temp_output_vcf,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
                mincov=args.mincov,
                threads=args.threads,
            )
            print("[" + str(datetime.now()) + "] DONE!")
            if not args.keep_files:
                if '/' not in args.output_prefix:
                    temp_folder = './'
                else:
                    temp_folder = '/'.join(
                        [x for x in args.output_prefix.split('/')][:-1]
                    ) + '/'
                for file_name in listdir(temp_folder):
                    if file_name.startswith(timestamp):
                        remove(f"{temp_folder}{file_name}")
    exit(0)
