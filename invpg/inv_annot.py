#! /bin/python3
from re import split as rexpsplit
from pathlib import Path
from typing import Any
from subprocess import run
from os import path
from shutil import rmtree


def parse_path(
    allele_path: str
) -> list[int]:
    """From VCF allele path, parse to get list of signed nodes 
    (sign = node traversal strand). 

    Parameters
    ----------
    allele_path : _type_
        _description_

    Returns
    -------
    list[int]
    Examples: 
        - ">43>45>46" --> [43,45,46]
        - ">43<44>46" --> [43,-44,46]
    """
    parsed_path: list[str] = rexpsplit(r'(\D+)', allele_path)[1:]

    int_path: list[int] = list()

    for i in range(0, len(parsed_path), 2):
        int_path.append(
            int(parsed_path[i+1]) if parsed_path[i] == ">" else -int(parsed_path[i+1]))

    return int_path


def is_INV_fromPath(
    a0_path: list[int],
    a1_path: list[int],
) -> tuple[bool, list[int]]:
    """ Checks whether paths from a bubble indicate the presence of an INVersion.

    Input: parsed path from parse_path() function
    Returns: boolean
    """
    inverted_nodes: list[int] = list()

    is_patternFound: bool = False

    a1_set: set[int] = set(a1_path)
    for i in a0_path:

        if -i in a1_set:
            is_patternFound = True
            inverted_nodes.append(i)

    return (is_patternFound, inverted_nodes)


def is_INV_fromAln(
    aln_paf: str
) -> tuple[bool, float, int, float, int]:
    """Checks whether the nodes sequences are reverse complement

    Parameters
    ----------
    aln_paf : str
        Path to an alignment file

    Returns
    -------
    tuple[bool, float, int, float, int]
        _description_
    """
    is_compRev: bool = False
    cum_for_len: int = 0
    cum_rev_len: int = 0
    n_for_aln: int = 0
    n_rev_aln: int = 0
    len_0: int = -1

    with open(aln_paf, 'r', encoding='utf-8') as file:
        for line in file:
            len_0, a_start, a_end, strand = line.split("\t")[1:5]

            if strand == "-" and int(a_end) - int(a_start) >= 50:
                cum_rev_len += int(a_end) - int(a_start)
                n_rev_aln += 1

            elif strand == "+" and int(a_end) - int(a_start) >= 50:
                cum_for_len += int(a_end) - int(a_start)
                n_for_aln += 1

    is_compRev = cum_rev_len/int(len_0) > cum_for_len/int(len_0)

    return (is_compRev, cum_rev_len/int(len_0), n_rev_aln, cum_for_len/int(len_0), n_for_aln)


def write_fasta(
    fasta_name: str,
    seq_id: str,
    sequence: str
) -> None:
    """_summary_

    Parameters
    ----------
    fasta_name : str
        _description_
    seq_id : str
        _description_
    sequence : str
        _description_
    """
    with open(fasta_name, 'w', encoding='utf-8') as fasta:
        fasta.write(f">{seq_id}\n{sequence}")


def get_len_node(
    d_nodes: dict[str, int],
    nodeID: int
) -> int:
    """_summary_

    Parameters
    ----------
    d_nodes : dict[str, int]
        _description_
    nodeID : int
        _description_

    Returns
    -------
    int
        _description_
    """
    try:
        return d_nodes[str_nodeID := str(abs(nodeID))]
    except KeyError:
        print(f"Error: node {str_nodeID} not found in GFA")

# ===========================================================
# MAIN
# ===========================================================


def invannot_main(
    gfa_file: str,
    vcf_file: str,
    temp_folder: str,
    threads: int,
) -> str:
    """_summary_

    Parameters
    ----------
    gfa_file : str
        _description_
    vcf_file : str
        _description_
    temp_folder : str
        _description_
    threads : int
        _description_
    output_file : str
        _description_
    """
    d_nodes: dict[str, int] = dict()
    with open(gfa_file, 'r', encoding='utf-8') as input_gfa_file:
        for line in input_gfa_file:
            if line.startswith("S"):
                __, nID, nSeq = line.rstrip().split("\t")[:3]
                d_nodes[nID] = len(nSeq)

    Path(temp_folder).mkdir(parents=True, exist_ok=True)

    with open(outBED := f"{path.splitext(vcf_file)[0]}.raw.bed", 'w', encoding='utf-8') as output_bed_file:
        with open(vcf_file, 'r', encoding='utf-8') as input_vcf_file:
            for line in input_vcf_file:

                if line[0] == "#":
                    continue

                # ---------------------------------------------------
                # Retrieve coordinates of bubble
                # ---------------------------------------------------
                chrom, pos, bubble = line.split("\t")[0:3]

                # ---------------------------------------------------
                # Retrieve allele paths and sequences from line
                # ---------------------------------------------------
                if ";" in line.split("\t")[7]:
                    aPaths: str = line.split("\t")[7].split("AT=")[
                        1].split(";")[0]
                else:
                    aPaths: str = line.rstrip().split("\t")[7].split("AT=")[1]

                a0Seq: str = line.split("\t")[3]
                a1Seqs: str = line.split("\t")[4]

                aPaths: str = aPaths.split(",")
                a1Seqs: str = a1Seqs.split(",")

                n_a1: int = len(a1Seqs)
                are_INV: list = [None] * n_a1

                # For potential alignment
                a0Fasta: str = f"{temp_folder}/{chrom}.{pos}.a0.fa"
                write_fasta(a0Fasta, "a0", a0Seq)

                # ---------------------------------------------------
                # Check if alleles present INV pattern
                # ---------------------------------------------------
                a0Path: list[int] = parse_path(aPaths[0])

                # Get balanced a1
                f_INFO: str = line.split("\t")[7]
                i_bal: list[int] = f_INFO.split(";BL=")[1]

                if "," in i_bal:
                    i_bal: list[int] = list(map(int, i_bal.split(",")))

                else:
                    i_bal: list[int] = [int(i_bal)]

                for i in i_bal:

                    # -----------------------------------------------
                    # Check for pattern in path
                    # -----------------------------------------------
                    a1Path = parse_path(aPaths[i])

                    is_inv_fromPath, rev_nodes = is_INV_fromPath(
                        a0Path, a1Path)

                    if is_inv_fromPath:

                        len_rev = 0
                        for n in rev_nodes:
                            len_rev += get_len_node(d_nodes, n)

                        are_INV[i-1] = (True, "path", ",".join([str(round(len_rev/len(a0Seq), 2)), str(
                            len(rev_nodes))]), str(len(a1Path)-2-len(rev_nodes)))

                        # print(chrom, pos, i)

                    # -----------------------------------------------
                    # Check for pattern in alignment
                    # -----------------------------------------------
                    else:

                        a1Seq = a1Seqs[i-1]

                        a1Fasta = f"{temp_folder}/{chrom}.{pos}.a{str(i+1)}.fa"
                        write_fasta(a1Fasta, "a1", a1Seq)
                        # run(f"echo '>a1' > {a1Fasta}", shell=True)
                        # run(f"echo {a1Seq} >> {a1Fasta}", shell=True)

                        # Run minimap2
                        alnPAF: str = f"{temp_folder}/{chrom}.{pos}.a{str(i+1)}.paf"
                        run(
                            f"minimap2 -cx asm20 --cs -r2k -t {threads} {a0Fasta} {a1Fasta} > {alnPAF} ",
                            shell=True,
                        )

                        is_inv_fromAln, frac_rev, n_rev_aln, frac_for, n_for_aln = is_INV_fromAln(
                            alnPAF)

                        if is_inv_fromAln:
                            are_INV[i-1] = (True, "aln", ",".join([str(round(frac_rev, 2)), str(
                                n_rev_aln)]), ",".join([str(round(frac_for, 2)), str(n_for_aln)]))
                        else:
                            are_INV[i-1] = (False, ".")

                for i in range(len(are_INV)):

                    if are_INV[i] == None:
                        are_INV[i] = (False, ".")

                # ---------------------------------------------------
                # Output results
                # ---------------------------------------------------
                output_bed_file.write("\t".join([
                    chrom, pos, bubble,
                    str(len(a0Seq)),
                    ",".join([str(len(a1)) for a1 in a1Seqs]),
                    ";".join(
                        ["INV:" + ":".join(b[1:]) if b[0]
                            else "DIV" for b in are_INV]
                    )
                ]) + "\n")

    return outBED
