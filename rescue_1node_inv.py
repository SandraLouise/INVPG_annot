#! /bin/python3

import sys
import subprocess

def str_path_to_int(str_path):

    int_path = []

    for str_node in str_path.split(','):

        if str_node[-1] == "+":
            int_path.append(int(str_node[:-1]))
        elif str_node[-1] == "-":
            int_path.append(- int(str_node[:-1]))

    return int_path

def int_path_to_str(int_path):

    str_nodes = []

    for int_node in int_path:

        if int_node > 0:
            str_node = f"{str(abs(int_node))}+"
        elif int_node < 0:
            str_node = f"{str(abs(int_node))}-"
        
        str_nodes.append(str_node)
    
    str_path = ",".join(str_nodes)

    return str_path

def parse_P_line(P_line):

    __, path_ID, str_path = P_line.rstrip().split("\t")[:3]

    if str_path[0] == "s":
        str_nodes = [n[1:] for n in str_path.split(",")]
        str_path = ",".join(str_nodes)
    
    int_path = str_path_to_int(str_path)
    
    return path_ID, str_path, int_path

def parse_W_line(W_line):

    split_line = W_line.rstrip().split("\t")
    pass

def find_rev_pattern(d_int_paths):
    """ Find initial pattern ('+x,-y,+z' in any path p) """

    list_IDs = list(d_int_paths.keys())

    d_rev_patterns = {}
    list_rev_duplicates = []

    for path_ID in list_IDs:

        int_path = d_int_paths[path_ID]

        for i_node in range(len(int_path) - 3):

            if all([int_path[i_node] > 0,
                    int_path[i_node+1] < 0,
                    int_path[i_node+2] > 0]):
                
                rev_pat_str = int_path_to_str(int_path[i_node : i_node+3])

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

def find_whole_pattern(d_str_paths, d_rev_patterns):
    """ Find whole pattern ('+x,+y,+z' in path != p) """

    d_inv_pattern = {}
    #key = rev_pat_str
    #value = [[path_ID_rev], [path_ID_for]]

    for rev_pat_str, list_path_ID_rev in d_rev_patterns.items():
        
        rev_pat_int = str_path_to_int(rev_pat_str)
        for_pat_int = [rev_pat_int[0], abs(rev_pat_int[1]), rev_pat_int[2]]
        for_pat_str = int_path_to_str(for_pat_int)

        # Look for for_pat_str in paths with ID != rev_list_IDs
        #-----------------------------------------------------------------------
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

def write_fasta(fasta_name, seq_id, sequence):

    with open(fasta_name, "w") as fasta:

        fasta.write(f">{seq_id}\n")
        fasta.write(sequence)

def is_INV_fromAln(aln_paf):
    """ Checks whether the nodes sequences are reverse complement

    Input: name of the aln file
    Returns: boolean
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

def search_gfa(in_gfa):

    d_node_len = {}
    d_str_paths = {}
    d_int_paths = {}

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
        node_id = abs(d_int_paths[REF_PATH][i])
        while node_id != inv_node_id:
            start += d_node_len[node_id]
            i += 1
            node_id = abs(d_int_paths[REF_PATH][i])

        end = start + d_node_len[inv_node_id] - 1

        # Output inv
        print("\t".join([REF_PATH, str(start), str(end), "INV:path"]))

def search_bed(in_bed):

    # tmpDir = "tmp_aln"
    # subprocess.run(f"mkdir {tmpDir}", shell=True)

    # outBed = "inv.bed"
    # out = open(outBed, "w")
    
    with open(in_bed, "r") as file:
        for line in file:

            parsed_line = line.rstrip().split("\t")
            chrom, pos, end = parsed_line[:3]
            size_bubble = int(parsed_line[3])
            common = int(parsed_line[5])
            a0Len, a1Len = parsed_line[6:8]
            bubble = parsed_line[11].split(",")

            if common == 1 and bubble[1] == bubble[2]:

                # out.write("\t".join([
                #     chrom, pos, end,
                #     a0Len, a1Len,
                #     ",".join(bubble),
                #     "INV:path"
                #     ]) + "\n")
                
                print("\t".join([
                    chrom, pos, end,
                    a0Len, a1Len,
                    ",".join(bubble),
                    "INV:path"
                    ]))

            # elif size_bubble > 3:
            #     a0Seq, a1Seq = parsed_line[12:14]

            #     a0Fasta = f"{tmpDir}/{chrom}.{pos}.a0.fa"
            #     write_fasta(a0Fasta, "a0", a0Seq)
            #     a1Fasta = f"{tmpDir}/{chrom}.{pos}.a1.fa"
            #     write_fasta(a1Fasta, "a1", a1Seq)

            #     alnPAF = f"{tmpDir}/{chrom}.{pos}.paf"
            #     c_minimap2 = f"minimap2 -cx asm20 --cs -r2k -t 1 {a0Fasta} {a1Fasta} > {alnPAF} "
            #     subprocess.run(c_minimap2, shell=True)
            #     is_inv_fromAln, frac_rev, __, frac_for, __ = is_INV_fromAln(alnPAF)

            #     if is_inv_fromAln:
            #         out.write("\t".join([
            #             chrom, pos, end,
            #             a0Len, a1Len,
            #             ",".join(bubble),
            #             f"INV:aln:{str(round(frac_rev, 2))}:{str(round(frac_for, 2))}"
            #             ]) + "\n")
    
    # out.close()
    # subprocess.run(f"rm -r {tmpDir}", shell=True)

def main(in_file, file_format):

    if file_format == "gfa":
        # print("Searching GFA file...")
        search_gfa(in_file)
    
    # elif file_format == "bed":
    #     print("Searching BED file...")
    #     search_bed(in_file)

#===============================================================================               
# MAIN
#=============================================================================== 

in_file, file_format, REF_PATH = sys.argv[1:]
main(in_file, file_format)