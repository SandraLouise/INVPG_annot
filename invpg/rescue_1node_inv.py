#! /bin/python3
from os import path
from pathlib import Path


def str_path_to_int(
    str_path: str
) -> list[int]:
    """Conerts a path of str to ints

    Parameters
    ----------
    str_path : str
        The P-line path field

    Returns
    -------
    list[int]
        a representation as a series of positive and negative ints of the string
    """
    return [
        int(str_node[:-1]) if str_node[-1] == "+" else -int(str_node[:-1]) for str_node in str_path.split(',')
    ]


def int_path_to_str(
    int_path: list[int]
) -> str:
    """Converts a path of ints to str format

    Parameters
    ----------
    int_path : list[int]
        the P-line path as a series of ints

    Returns
    -------
    str
        the character chain describing the same path
    """
    return ",".join(
        [
            f'{str(abs(int_node))}+' if int_node > 0 else f'{str(abs(int_node))}-' for int_node in int_path
        ]
    )


def parse_P_line(
    P_line: str
) -> tuple[str, str, list[int]]:
    """Parses a P-line from GFAspec

    Parameters
    ----------
    P_line : str
        Full line from a GFA file

    Returns
    -------
    tuple[str, str, list[int]]
        formatted information from the line
    """

    path_ID, str_path = P_line.rstrip().split("\t")[1:3]

    if str_path[0] == "s":
        str_nodes = [n[1:] for n in str_path.split(",")]
        str_path = ",".join(str_nodes)

    int_path: list[int] = str_path_to_int(str_path)

    return (path_ID, str_path, int_path)


def find_rev_pattern(
    d_int_paths: dict[str, list[int]]
):
    """Find initial pattern ('+x,-y,+z' in any path p)

    Parameters
    ----------
    d_int_paths : dict[str, list[int]]
        Dictionary with path_ID as key and list of int nodes as value.
    """

    list_IDs = list(d_int_paths.keys())

    d_rev_patterns = {}
    list_rev_duplicates = []

    for path_ID in list_IDs:

        int_path = d_int_paths[path_ID]

        for i_node in range(len(int_path) - 3):

            if all([int_path[i_node] > 0,
                    int_path[i_node+1] < 0,
                    int_path[i_node+2] > 0]):

                rev_pat_str = int_path_to_str(int_path[i_node: i_node+3])

                if rev_pat_str not in d_rev_patterns.keys():
                    d_rev_patterns[rev_pat_str] = []

                # If rev_pat_str duplicated in same path, remove
                elif path_ID in d_rev_patterns[rev_pat_str]:
                    list_rev_duplicates.append(rev_pat_str)
                    del d_rev_patterns[rev_pat_str]

                if rev_pat_str not in list_rev_duplicates:
                    d_rev_patterns[rev_pat_str].append(path_ID)

    rev_pat_str = list(d_rev_patterns.keys())

    for pat_str in rev_pat_str:
        if len(d_rev_patterns[pat_str]) == len(list_IDs):
            del d_rev_patterns[pat_str]

    return d_rev_patterns


def find_whole_pattern(
    d_str_paths: dict,
    d_rev_patterns: dict,
):
    """Find whole pattern ('+x,+y,+z' in path != p)

    Parameters
    ----------
    d_str_paths : dict
        dictionnary of paths
    d_rev_patterns : dict
        reverse patterns

    Returns
    -------
    dict
        a subselection of reverse patterns filtered using d_str_paths
    """

    d_inv_pattern = {}
    # key = rev_pat_str
    # value = [[path_ID_rev], [path_ID_for]]

    for rev_pat_str, list_path_ID_rev in d_rev_patterns.items():

        rev_pat_int = str_path_to_int(rev_pat_str)
        for_pat_int = [rev_pat_int[0], abs(rev_pat_int[1]), rev_pat_int[2]]
        for_pat_str = int_path_to_str(for_pat_int)

        # Look for for_pat_str in paths with ID != rev_list_IDs
        # -----------------------------------------------------------------------
        paths_to_search = list(d_str_paths.keys())

        for path_ID in list_path_ID_rev:

            if path_ID not in paths_to_search:
                print(rev_pat_str, list_path_ID_rev)

            paths_to_search.remove(path_ID)

        for p in paths_to_search:

            if for_pat_str in d_str_paths[p]:

                if rev_pat_str not in d_inv_pattern.keys():
                    d_inv_pattern[rev_pat_str] = [list_path_ID_rev, []]

                d_inv_pattern[rev_pat_str][1].append(p)

    return d_inv_pattern


def write_fasta(
    fasta_name: str,
    seq_id: str,
    sequence: str,
):
    """Writes to disk a sequence as a .fasta file

    Parameters
    ----------
    fasta_name : str
        path + name of the fasta file to be written
    seq_id : str
        ID of the fasta sequence
    sequence : str
        Nucleotidic sequence
    """

    with open(fasta_name, "w") as fasta:

        fasta.write(f">{seq_id}\n")
        fasta.write(sequence)


def is_INV_fromAln(
    aln_paf: str,
) -> bool:
    """Checks whether the nodes sequences are reverse complement

    Parameters
    ----------
    aln_paf : str
        Path to a alignment file

    Returns
    -------
    bool
        nodes are reverse complement
    """

    is_compRev = False
    cum_for_len = 0
    cum_rev_len = 0
    n_for_aln = 0
    n_rev_aln = 0
    len_0 = -1

    with open(aln_paf, "r") as file:
        for line in file:

            len_0, a_start, a_end, strand = line.split("\t")[1:5]

            if strand == "-" and int(a_end) - int(a_start) >= 50:
                cum_rev_len += int(a_end) - int(a_start)
                n_rev_aln += 1

            elif strand == "+" and int(a_end) - int(a_start) >= 50:
                cum_for_len += int(a_end) - int(a_start)
                n_for_aln += 1

    if cum_rev_len/int(len_0) > cum_for_len/int(len_0):
        is_compRev = True

    return is_compRev, cum_rev_len/int(len_0), n_rev_aln, cum_for_len/int(len_0), n_for_aln


def search_gfa(
    in_gfa: str,
    ref_path: str,
) -> None:
    """Search for 1-sized in GFA file

    Parameters
    ----------
    in_gfa : str
        Ipnut GFA file path
    ref_path : str
        Name (ID) of the reference path in the graph
    """

    d_node_len: dict = {}
    d_str_paths: dict = {}
    d_int_paths: dict = {}

    with open(f"{path.splitext(in_gfa)[0]}.complementary.bed"):
        with open(in_gfa, "r") as file:

            for line in file:
                if line.startswith("S"):
                    __, node_id, node_seq = line.rstrip().split("\t")[:3]
                    d_node_len[int(node_id)] = len(node_seq)

                elif line.startswith("P"):

                    path_ID, str_path, int_path = parse_P_line(line)

                    d_str_paths[path_ID] = str_path
                    d_int_paths[path_ID] = int_path

        d_rev_patterns = find_rev_pattern(d_int_paths)

        d_inv_patterns = find_whole_pattern(d_str_paths, d_rev_patterns)

        # print(f"#INV patterns found: {len(d_inv_patterns.keys())}")

        for inv_pattern, path_IDs in d_inv_patterns.items():

            # Get inv node id
            inv_node_id = abs(str_path_to_int(inv_pattern)[1])

            # Get start & end of inversion
            i = 0
            start = 0
            node_id = abs(d_int_paths[ref_path][i])
            while node_id != inv_node_id:
                start += d_node_len[node_id]
                i += 1
                node_id = abs(d_int_paths[ref_path][i])

            end = start + d_node_len[inv_node_id] - 1

            # Output inv
            print("\t".join([ref_path, str(start), str(end), "INV:path"]))


def search_bed(
    in_bed: str,
    output_prefix: str,
    timestamp: str,
):
    """Detects one-node inversions that may be missing from vg deconstruct VCF.

    Parameters
    ----------
    in_bed : str
        path to a input bed file
    output_prefix : str
        path to a folder to store temporary files and results
    timestamp : str
        timecode for temporary files
    """
    # Path for temporary files
    if '/' not in output_prefix:
        temp_folder = './'
    else:
        temp_folder = '/'.join(
            [x for x in output_prefix.split('/')][:-1]
        ) + '/'

    Path(temp_folder).mkdir(parents=True, exist_ok=True)

    with open(f"{temp_folder}{timestamp}.rescue1node.bed") as out_bed:
        with open(in_bed, "r", encoding='utf-8') as file:
            for line in file:

                parsed_line = line.rstrip().split("\t")
                chrom, pos, end = parsed_line[:3]
                size_bubble = int(parsed_line[3])
                common = int(parsed_line[5])
                a0Len, a1Len = parsed_line[6:8]
                bubble = parsed_line[11].split(",")

                if common == 1 and bubble[1] == bubble[2]:

                    print("\t".join([
                        chrom, pos, end,
                        a0Len, a1Len,
                        ",".join(bubble),
                        "INV:path",
                    ]), file=out_bed)
