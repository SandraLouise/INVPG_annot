#! /bin/python3

import sys
import subprocess
import re

inVCF, inGFA, inPREF, threads, MINCOV = sys.argv[1:]

MINCOV = float(MINCOV)

tmpDir = f"tmp_{inPREF}"
outName = f"{inPREF}.annot.bed"

#===========================================================
# FUNCTIONS
#===========================================================

def parse_path(allele_path):
    """ From VCF allele path, parse to get list of signed nodes 
    (sign = node traversal strand). 
    
    Returns: list of ints

    Examples: 
        - ">43>45>46" --> [43,45,46]
        - ">43<44>46" --> [43,-44,46]
    """
    parsed_path = re.split(r'(\D+)', allele_path)[1:]

    int_path = []

    for i in range(0, len(parsed_path), 2):

        if parsed_path[i] == ">":
            int_path.append(int(parsed_path[i+1]))
        
        else:
            int_path.append(-int(parsed_path[i+1]))

    return int_path

def is_INV_fromPath(a0_path, a1_path):
    """ Checks whether paths from a bubble indicate the presence of an INVersion.
    
    Input: parsed path from parse_path() function
    Returns: boolean
    """
    inverted_nodes = []

    is_patternFound = False

    for i in a0_path:

        if -i in a1_path:
            is_patternFound = True
            inverted_nodes.append(i)
    
    #-------------------------------------------------------
    # IDEA: better check of INV pattern
    #
    # 1) Save positions of i and -i in a0 and a1 paths
    # 2) Check that:
    #    if pos(i1) < pos(i2) in a0:
    #       pos(-i1) > pos(-i2) in a1
    #-------------------------------------------------------

    return is_patternFound, inverted_nodes

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

def write_fasta(fasta_name, seq_id, sequence):

    with open(fasta_name, "w") as fasta:

        fasta.write(f">{seq_id}\n")
        fasta.write(sequence)

def get_len_node(d_nodes, nodeID):

    str_nodeID = str(abs(nodeID))
    
    if str_nodeID in d_nodes.keys():
        return d_nodes[str_nodeID]
    else:
        print(f"Error: node {str_nodeID} not found in GFA")
        exit


#===========================================================
# MAIN
#===========================================================
        
d_nodes = {}
with open(inGFA, "r") as file:
    for line in file:
        if line[0] == "S":
            __, nID, nSeq = line.rstrip().split("\t")[:3]
            d_nodes[nID] = len(nSeq)

subprocess.run(f"mkdir {tmpDir}", shell=True)

out = open(outName, "w")

with open(inVCF, "r") as file:

    for line in file:

        if line[0] == "#":
            continue

        #---------------------------------------------------
        # Retrieve coordinates of bubble
        #---------------------------------------------------
        chrom, pos, bubble = line.split("\t")[0:3]

        #---------------------------------------------------
        # Retrieve allele paths and sequences from line
        #---------------------------------------------------
        if ";" in line.split("\t")[7]:
            aPaths = line.split("\t")[7].split("AT=")[1].split(";")[0]
        else:
            aPaths = line.rstrip().split("\t")[7].split("AT=")[1]

        a0Seq = line.split("\t")[3]
        a1Seqs = line.split("\t")[4]

        aPaths = aPaths.split(",")
        a1Seqs = a1Seqs.split(",")

        n_a1 = len(a1Seqs)
        are_INV = [None] * n_a1

        # For potential alignment
        a0Fasta = f"{tmpDir}/{chrom}.{pos}.a0.fa"
        write_fasta(a0Fasta, "a0", a0Seq)

        #---------------------------------------------------
        # Check if alleles present INV pattern
        #---------------------------------------------------
        a0Path = parse_path(aPaths[0])

        # Get balanced a1
        f_INFO = line.split("\t")[7]
        i_bal = f_INFO.split(";BL=")[1]

        if "," in i_bal:
            i_bal = list(map(int, i_bal.split(",")))

        else:
            i_bal = [int(i_bal)]

        for i in i_bal:

            #-----------------------------------------------
            # Check for pattern in path
            #-----------------------------------------------
            a1Path = parse_path(aPaths[i])

            is_inv_fromPath, rev_nodes = is_INV_fromPath(a0Path, a1Path)
            path_coverage = 0

            if is_inv_fromPath:
                
                len_rev = 0
                for n in rev_nodes:
                    len_rev += get_len_node(d_nodes, n)
                
                path_coverage = round(len_rev/len(a0Seq), 2)
            
            if path_coverage >= MINCOV:
                are_INV[i-1] = (True, "path", ",".join([str(path_coverage), str(len(rev_nodes))]), str(len(a1Path)-2-len(rev_nodes)))
            
            else:

                #-----------------------------------------------
                # Check for pattern in alignment
                #-----------------------------------------------
                a1Seq = a1Seqs[i-1]

                a1Fasta = f"{tmpDir}/{chrom}.{pos}.a{str(i+1)}.fa"
                write_fasta(a1Fasta, "a1", a1Seq)

                # Run minimap2
                alnPAF = f"{tmpDir}/{chrom}.{pos}.a{str(i+1)}.paf"
                c_minimap2 = f"minimap2 -cx asm20 --cs -r2k -t {threads} {a0Fasta} {a1Fasta} > {alnPAF} "
                subprocess.run(c_minimap2, shell=True)

                is_inv_fromAln, frac_rev, n_rev_aln, frac_for, n_for_aln = is_INV_fromAln(alnPAF)
                aln_coverage = 0

                if is_inv_fromAln:
                    aln_coverage = round(frac_rev, 2)
                
                if aln_coverage >= MINCOV:
                    are_INV[i-1] = (True, "aln", ",".join([str(aln_coverage), str(n_rev_aln)]), ",".join([str(round(frac_for, 2)), str(n_for_aln)]))
                
                else:
                    are_INV[i-1] = (False, ".")

            # Delete paf file
            subprocess.run(f"rm {alnPAF}", shell=True)   

        for i in range(len(are_INV)):

            if are_INV[i] == None:
                are_INV[i] = (False, ".")

        if any([b[0] for b in are_INV]):

            out.write("\t".join([
                chrom, pos,
                str(int(pos) + len(a0Seq) - 1),
                ";".join(["INV:" + ":".join(b[1:]) if b[0] else "DIV" for b in are_INV])
            ]) + "\n")

        #---------------------------------------------------
        # Clean tmp fasta files
        #---------------------------------------------------
        subprocess.run(f"rm {tmpDir}/*.fa", shell=True)

out.close()

#-----------------------------------------------------------
# Clean tmp files
#-----------------------------------------------------------
# subprocess.run(f"rm -r {tmpDir}", shell=True)         
                         