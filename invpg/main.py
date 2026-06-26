#!/usr/bin/env python3

"""*******************************************************************************
    Name: INVPG-annot
    Description: INVPG-annot aims to annotate inversions among bubbles extracted from a pangenome graph.
    Authors: Sandra Romain, Siegfried Dubois
    Contact: claire.lemaitre@inria.fr, Inria/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
    
    Copyright (C) 2024 Inria
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

from invpg.__constants__ import *
from argparse import ArgumentParser
from sys import argv
from os import path
from shutil import rmtree, which
from datetime import datetime
from invpg.inv_annot import invannot
from invpg.variant_filter import variant_filter

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
    default=DEFAULT_OUTPUT_MAIN,
)
parser.add_argument(
    "-d",
    "--div_percentage",
    type=int,
    help=HELP_PARAM_PERCENTAGE,
    default=DEFAULT_PERCENTAGE,
)
parser.add_argument(
    "-m",
    "--mincov",
    type=float,
    help=HELP_PARAM_MINCOV,
    default=DEFAULT_MINCOV,
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
    default=DEFAULT_OUTPUT,
)
parser_invannot.add_argument(
    "-m",
    "--mincov",
    type=float,
    help=HELP_PARAM_MINCOV,
    default=DEFAULT_MINCOV,
)
parser_invannot.add_argument(
    "-t",
    "--threads",
    type=int,
    help=HELP_PARAM_THREADS,
    default=DEFAULT_THREADS,
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
    default=DEFAULT_OUTPUT,
)
parser_filtervcf.add_argument(
    "-d",
    "--div_percentage",
    type=int,
    help=HELP_PARAM_PERCENTAGE,
    default=DEFAULT_PERCENTAGE,
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

    # This timestamp helps identify temporary files of this run (for deletion later on)
    ts = datetime.now()
    print(f"Starting job @{str(ts)}")
    timestamp: str = str(ts).replace(' ', '_').replace(':', '')

    if not which('minimap2'):
        raise RuntimeError('Minimap2 is not installed or not in path.')

    match args.subcommands:
        case 'annot':
            invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=args.input_vcf_file,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
                mincov=args.mincov,
                threads=args.threads,
            )
        case 'filtvcf':
            filter(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
        case _:
            # We check if both -v and -g are given
            if not args.input_vcf_file or not args.input_gfa_file:
                raise RuntimeError('Both -v (input vcf file) and -g (input gfa file) must be given when using global command.')
            # First we filter the VCF file
            print("[" + str(datetime.now()) + "] STEP 1: filtering VCF file")
            temp_output_vcf: str = variant_filter(
                in_vcf=args.input_vcf_file,
                div_pct=args.div_percentage,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
            )
            if args.keep_files:
                print(f"Selected bubbles output in {temp_output_vcf}")
            # Then we rescue nodes in inversions that weren't described in the VCF
            print("[" + str(datetime.now()) +
                  "] STEP 2: rescuing nodes in inversions"
                  )
            invannot(
                gfa_file=args.input_gfa_file,
                vcf_file=temp_output_vcf,
                output_prefix=args.output_prefix,
                timestamp=timestamp,
                mincov=args.mincov,
                threads=args.threads,
            )
            print("[" + str(datetime.now()) + "] DONE!")
            print(f"Results output in files {path.splitext(args.output_prefix)[0]}.vcf and {path.splitext(args.output_prefix)[0]}.stats")
            if not args.keep_files:
                if '/' not in args.output_prefix:
                    temp_folder = f'./res_{timestamp}/'
                else:
                    temp_folder = '/'.join(
                        [x for x in args.output_prefix.split('/')][:-1]
                    ) + f'/res_{timestamp}/'
                rmtree(temp_folder)

    exit(0)
