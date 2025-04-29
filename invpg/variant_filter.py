#! /bin/python3

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

from typing import Any
from os import path
from pathlib import Path


def allele_length_mg(info: str):
    """Extracts reference and alternative allele lengths from VCF produced by gfatools-minigraph pipeline.

    Parameters
    ----------
    info : str
        INFO field of VCF line

    Returns
    -------
    int
        Reference allele length
    list[int]
        List of alternative allele(s) length
    """

    alen: str = info.split("ALEN=")[1]

    if ";" in alen:
        alen: str = alen.split(";")[0]

    if "," in alen:
        list_alen: list[str] = alen.split(",")

        ref_len: int = int(list_alen[0])
        alt_len: list[int] = list(map(int, list_alen[1:]))
    else:
        ref_len: int = int(alen)
        alt_len: list[int] = [0]

    return ref_len, alt_len


def allele_length_vg(
    ref: str,
    alt: str
):
    """Extracts reference and alternative allele lengths from VCF produced by vg deconstruct.

    Parameters
    ----------
    ref : str
        REF field of VCF line
    alt : str
        ALT field of VCF line

    Returns
    -------
    int
        Reference allele length
    list[int]
        List of alternative allele(s) length
    """

    ref_len: int = len(ref)

    if "," in alt:
        list_alt: list[str] = alt.split(",")
    else:
        list_alt: list[str] = [alt]

    alt_len: list[int] = [len(a) for a in list_alt]

    return ref_len, alt_len


def allele_length(
    ref: str,
    alt: str,
    info: str
):
    """Extracts reference and alternative allele lengths from VCF line.

    Parameters
    ----------
    ref : str
        REF field of VCF line
    alt : str
        ALT field of VCF line
    info : str
        INFO field of VCF line

    Returns
    -------
    list[int, [int]]
        List of reference allele length [0] and of alternative allele(s) length list [1]
    """

    if "ALEN=" in info:
        return allele_length_mg(info)

    else:
        return allele_length_vg(ref, alt)


def parse_vcf_line(line: str) -> dict[str, Any]:
    """Parse a line from a VCF file

    Parameters
    ----------
    line : str
        the raw line from the input VCF

    Returns
    -------
    dict[str, Any]
        formatted information
    """

    tab_parsed: list[str] = line.rstrip().split("\t")
    ref: str = tab_parsed[3]
    alt: str = tab_parsed[4]
    info: str = tab_parsed[7]

    ref_len, alt_len = allele_length(ref, alt, info)

    return {
        "ref_id": tab_parsed[0],
        "ref_start": int(tab_parsed[1]),
        "ref_end": int(tab_parsed[1]) + len(ref) - 1,
        "ref_seq": ref,
        "alt_seq": alt,
        "ref_len": ref_len,
        "alt_len": alt_len,
        "lvl": int(tab_parsed[7].split(";LV=")[1].split(";")[0]) if ";LV=" in tab_parsed[7] else "NA"
    }


def is_balanced(
    ref_len: int,
    alt_len: int,
    div: int,
) -> bool:
    """Analyses if the inversion is balanced or not

    Parameters
    ----------
    ref_len : int
        length on the reference path
    alt_len : int
        length on the alternate path

    Returns
    -------
    bool
        if the inversion is balanced
    """
    max_a: int = max(ref_len, alt_len)
    min_a: int = min(ref_len, alt_len)

    return max_a - min_a <= (div * max_a)


def variant_filter(
    in_vcf: str,
    div_pct: int,
    output_prefix: str,
    timestamp: str,
) -> str:
    """Selects bubbles corresponding to putative balanced SVs.

    Parameters
    ----------
    in_vcf : str
        input VCF path to file
    div_pct : int
        Originally intended as the estimated percentage of genome divergence. This parameter controls the leniency of the algorithm towards allele size difference (in nt) in the first step of variant/bubble filtering.

    Returns
    -------
    str
        Path to the balanced SV VCF file.
    """
    count_BL_entries: int = 0
    div: float = float(div_pct) / 100
    min_len: int = 50
    f_INFO: int = 7

    # Path for temporary files
    if '/' not in output_prefix:
        temp_folder = './'
    else:
        temp_folder = '/'.join(
            [x for x in output_prefix.split('/')][:-1]
        ) + '/'

    Path(temp_folder).mkdir(parents=True, exist_ok=True)

    with open(outVCF := f"{temp_folder}{timestamp}_balancedSV.vcf", 'w', encoding='utf-8') as out_vcf_balanced:
        with open(in_vcf, 'r', encoding='utf-8') as file:
            for line in file:

                if line.startswith("#"):
                    out_vcf_balanced.write(line)
                    continue

                parsed_line: dict[str, Any] = parse_vcf_line(line)

                if parsed_line["ref_len"] < min_len and all([alt < min_len for alt in parsed_line["alt_len"]]):
                    continue

                i_bal: list[str] = list()

                for i_alt in range(len(parsed_line["alt_len"])):

                    alt_len: int = parsed_line["alt_len"][i_alt]

                    if is_balanced(
                        ref_len=parsed_line["ref_len"],
                        alt_len=alt_len,
                        div=div,
                    ):
                        i_bal.append(str(i_alt + 1))

                # ---------------------------------------------------
                # Output balanced SV
                # ---------------------------------------------------
                if len(i_bal) > 0:

                    vcf_array: list[str] = line.rstrip().split("\t")

                    bal_INFO: str = f"BL={','.join(i_bal)}"
                    new_INFO: str = ";".join([vcf_array[7], bal_INFO])

                    vcf_array[f_INFO] = new_INFO
                    out_vcf_balanced.write("\t".join(vcf_array)+"\n")

                    count_BL_entries += 1

    print(f"Balanced variants entries: {str(count_BL_entries)}")
    return outVCF
