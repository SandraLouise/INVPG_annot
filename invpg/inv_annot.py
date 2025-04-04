#! /bin/python3
import re
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
            int(parsed_path[i+1]) if parsed_path[i][0] == ">" else -int(parsed_path[i+1]))

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


def get_node_len(
    d_nodes: dict[str, str],
    nodeID: str
) -> int:
    """_summary_

    Parameters
    ----------
    d_nodes : dict[str, str]
        Dictionnary containing node ID as keys and node sequence as values.
    nodeID : str
        ID of the target node to retrieve the length of.

    Returns
    -------
    int
        Length of the target node.
    """
    try:
        return len(d_nodes[str_nodeID := str(abs(nodeID))])

    except KeyError:
        print(f"Error: node {str_nodeID} not found in GFA")


def index_node_seq(
    d_nodes: dict[str, str],
    Sline: str
) -> dict[str, str]:
    """Index the node ID with their sequence in a dictionnary.

    Parameters
    ----------
    d_nodes : dict[str, str]
        Dictionnary containing node ID as keys and node sequence as values.
    Sline : str
        S line of a GFA.

    Returns
    -------
    dict[str, str]
        Updated index dictionnary.
    """

    __, nID, nSeq = Sline.rstrip().split("\t")[:3]
    nID = str(re.sub("[^0-9]", "", nID))

    d_nodes[nID] = nSeq

    return d_nodes


def allele_walk(info: str):
    """Retrieve the alleles walk in the bubble.

    Parameters
    ----------
    info : str
        INFO field of a VCF line.

    Returns
    -------
    list[str]
        List of allele walks.
    """
    # vg deconstruct
    if "AT=" in info:
        aWalks: str = info.split("AT=")[1]

        if ";" in aWalks:
            aWalks: str = aWalks.split(";")[0]

        aWalks: list[str] = aWalks.split(",")

        # remove source and sink nodes from vg allele walks
        corrected_aWalks: list[str] = list()
        for w in aWalks:
            int_walk: list[int] = parse_path(w)
            int_walk = int_walk[1:-1]

            str_walk: str = ""
            for i in range(0, len(int_walk)):
                if int_walk[i] > 0:
                    str_walk = str_walk + ">{}".format(str(int_walk[i]))
                else:
                    str_walk = str_walk + "<{}".format(str(abs(int_walk[i])))

            corrected_aWalks.append(str_walk)
        return corrected_aWalks

    # gfatools-minigraph pipeline
    elif "AWALK=" in info:
        aWalks: str = info.split("AWALK=")[1]

        if ";" in aWalks:
            aWalks: str = aWalks.split(";")[0]

        aWalks: list[str] = aWalks.split(",")
        return aWalks


def get_node_seq(
    d_nodes: dict[str, str],
    node_ID: str
) -> str:
    """Retrieve sequence of a node from index dictionnary.

    Parameters
    ----------
    d_nodes : dict[str, str]
        Dictionnary containing node ID as keys and node sequence as values.
    node_ID : str
        ID of the target node.

    Returns
    -------
    str
        Sequence of the target node.
    """

    return d_nodes[node_ID]


def get_allele_seq(
    allele_walk: list[int],
    d_nodes: dict[str, str]
) -> str:
    """Recontruct the sequence of an allele from its walk through the bubble.

    Parameters
    ----------
    allele_walk : list[int]
        Parsed walk of an allele through its bubble (as an int list).
    d_nodes : dict[str, str]
        Dictionnary containing node ID as keys and node sequence as values.

    Returns
    -------
    str
        Sequence of the target node.
    """

    allele_seq: str = ""

    for node_int in allele_walk:

        # Forward traversal
        if node_int > 0:
            allele_seq = allele_seq + d_nodes[str(node_int)]

        # Reverse traversal
        else:
            allele_seq = allele_seq + \
                reverse_complement(d_nodes[str(abs(node_int))])

    return allele_seq


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement sequence of a sequence.

    Parameters
    ----------
    seq : str
        Original sequence

    Returns
    -------
    str
        Reverse complement sequence
    """

    d = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }

    revcomp: str = ""

    for pos in range(len(seq)-1, -1, -1):
        revcomp = revcomp + d[seq[pos]]

    return revcomp


def invannot(
    gfa_file: str,
    vcf_file: str,
    output_prefix: str,
    timestamp: str,
    mincov: float,
    threads: int,
) -> None:
    """_summary_

    Parameters
    ----------
    gfa_file : str
        Path to a valid GFA file
    vcf_file : str
        Path to a valid VCF file
    output_prefix : str
        Path to a .bed output file
    mincov : float
        Minimum coverage
    threads : int
        Number of threads for minimap2
    """
    d_nodes: dict[str, str] = dict()
    with open(gfa_file, 'r', encoding='utf-8') as input_gfa_file:
        for line in input_gfa_file:
            if line.startswith("S"):
                d_nodes = index_node_seq(d_nodes, line)

    # Defining output path
    if not output_prefix.endswith('.bed'):
        output_bed_file = output_prefix + '.bed'
    else:
        output_bed_file = output_prefix

    # Path for temporary files
    if '/' not in output_prefix:
        temp_folder = './'
    else:
        temp_folder = '/'.join(
            [x for x in output_prefix.split('/')][:-1]
        ) + '/'

    Path(temp_folder).mkdir(parents=True, exist_ok=True)

    bubble_count: int = 0

    with open(output_bed_file, 'w', encoding='utf-8') as output_bed_file:
        with open(vcf_file, 'r', encoding='utf-8') as input_vcf_file:
            for line in input_vcf_file:

                if line[0] == "#":
                    continue

                # ---------------------------------------------------
                # Retrieve coordinates of bubble
                # ---------------------------------------------------
                chrom, pos = line.split("\t")[0:2]
                bubble_count += 1
                # ---------------------------------------------------
                # Retrieve allele paths and sequences from line
                # ---------------------------------------------------
                info: str = line.split("\t")[7]
                aWalks: list[str] = allele_walk(info)

                a0Walk: list[int] = parse_path(aWalks[0])
                a0Seq: str = get_allele_seq(a0Walk, d_nodes)

                n_a1: int = len(aWalks) - 1
                are_INV: list = [None] * n_a1

                # For potential alignment
                a0Fasta: str = f"{temp_folder}{timestamp}{chrom}.{pos}.a0.fa"
                write_fasta(a0Fasta, "a0", a0Seq)

                # ---------------------------------------------------
                # Check if alleles present INV pattern
                # ---------------------------------------------------

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
                    a1Walk = parse_path(aWalks[i])

                    is_inv_fromPath, rev_nodes = is_INV_fromPath(
                        a0Walk, a1Walk)
                    path_coverage: float = .0

                    if is_inv_fromPath:

                        len_rev: int = 0
                        for n in rev_nodes:
                            len_rev += get_node_len(d_nodes, n)

                        path_coverage: float = float(len_rev)/float(len(a0Seq))

                    if path_coverage >= mincov:
                        are_INV[i-1] = (True, "path", ",".join([str(path_coverage),
                                        str(len(rev_nodes))]), str(len(a1Walk)-len(rev_nodes)))

                    # -----------------------------------------------
                    # Check for pattern in alignment
                    # -----------------------------------------------
                    else:

                        # a1Seq = a1Seqs[i-1]
                        a1Seq = get_allele_seq(a1Walk, d_nodes)

                        a1Fasta = f"{temp_folder}{timestamp}{chrom}.{pos}.a{str(i+1)}.fa"
                        write_fasta(a1Fasta, "a1", a1Seq)

                        # Run minimap2
                        alnPAF: str = f"{temp_folder}{timestamp}{chrom}.{pos}.a{str(i+1)}.paf"
                        run(
                            f"minimap2 -cx asm20 --cs -r2k -t {threads} {a0Fasta} {a1Fasta} 1> {alnPAF} 2> /dev/null ",
                            shell=True,
                        )

                        is_inv_fromAln, frac_rev, n_rev_aln, frac_for, n_for_aln = is_INV_fromAln(
                            alnPAF)
                        aln_coverage: float = .0

                        if is_inv_fromAln:
                            aln_coverage = float(frac_rev)

                        if aln_coverage >= mincov:
                            are_INV[i-1] = (True, "aln", ",".join([str(aln_coverage), str(
                                n_rev_aln)]), ",".join([str(round(frac_for, 2)), str(n_for_aln)]))
                        else:
                            are_INV[i-1] = (False, ".")

                for i in range(len(are_INV)):

                    if are_INV[i] == None:
                        are_INV[i] = (False, ".")

                # ---------------------------------------------------
                # Output results
                # ---------------------------------------------------
                if any([b[0] for b in are_INV]):
                    output_bed_file.write(
                        "\t".join([
                            chrom,
                            pos,
                            str(int(pos) + len(a0Seq) - 1),
                            ";".join(
                                [
                                    "INV:" + ":".join(b[1:]) if b[0] else "DIV" for b in are_INV
                                ]
                            )
                        ]) + "\n"
                    )
            print("Total number of bubbles: " + str(bubble_count))
