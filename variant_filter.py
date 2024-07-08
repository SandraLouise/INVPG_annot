#! /bin/python3

""" Outputs result in stdout """

import sys
# import statistics as stat

in_vcf, out_pref, div_pct = sys.argv[1:]

#===========================================================
# Functions
#===========================================================
def parse_vcf_line(line):

    d = {
        "ref_id" : "",
        "ref_start" : -1,
        "ref_end" : -1,
        "ref_seq" : -1,
        "alt_seq" : -1,
        "ref_len" : -1,
        "alt_len" : [],
        "lvl" : -1
    }

    tab_parsed = line.rstrip().split("\t")

    d["ref_id"] = tab_parsed[0]
    d["ref_start"] = int(tab_parsed[1])

    #-------------------------------------------------------
    # Get alleles lengths and end on ref
    #-------------------------------------------------------
    ref = tab_parsed[3]
    alt = tab_parsed[4]

    # Parse alt alleles if more than one
    if "," in alt:
        alt = alt.split(",")
    else:
        alt = [alt]

    d["ref_seq"] = ref
    d["alt_seq"] = alt

    # Calculate lengths
    ref_len = len(ref)
    alt_len = [len(a) for a in alt]

    # Save lengths
    d["ref_len"] = ref_len
    d["alt_len"] = alt_len

    # Calculate end pos on ref
    d["ref_end"] = d["ref_start"] + ref_len -1

    #-------------------------------------------------------
    # Get bubble level
    #-------------------------------------------------------
    if ";LV=" in tab_parsed[7]:
        d["lvl"] = int(tab_parsed[7].split(";LV=")[1].split(";")[0])
    
    else:
        d["lvl"] = "NA"

    return d

def is_balanced(ref_len, alt_len):

    balanced = False

    max_a = max(ref_len, alt_len)
    min_a = min(ref_len, alt_len)

    if max_a - min_a <= (div * max_a):
        balanced = True
    
    return balanced

#===========================================================
# Main
#===========================================================

#-----------------------------------------------------------
# Statistics
count_BL_entries = 0
#-----------------------------------------------------------

# outBED = f"{out_pref}.sv.typed.bed"
# out_4intersect = f"{out_pref}.forIntersect.bed"
outVCF = f"{out_pref}.balancedSV.vcf"
# outVCF = f"{out_pref}.sv.vcf"

# out_bed = open(outBED, "w")
# out_intersect = open(out_4intersect, "w")
out_vcf_balanced = open(outVCF, "w")
# out_vcf = open(outVCF, "w")

with open(in_vcf, "r") as file:

    for line in file:

        #---------------------------------------------------
        # Parse vcf line
        #---------------------------------------------------
        if line.startswith("#"):
            out_vcf_balanced.write(line)
            continue

        # Ignore variants with conflict 
        # (if more traversals were found than expected from given ploidy)
        # potential tandem repeats ?
        # if ";CONFLICT=" in line:
        #     continue

        parsed_line = parse_vcf_line(line)

        # Same for variants with more than 3 alt alleles
        # if len(parsed_line["a1_len"]) > 3:
        #     continue

        #---------------------------------------------------
        # Filter SVs (at least 51 bp for either a0 or a1)
        #---------------------------------------------------
        min_len = 50
        if parsed_line["ref_len"] <= min_len and all([alt <= min_len for alt in parsed_line["alt_len"]]):
            continue

        #---------------------------------------------------
        # Filter balanced SV
        #---------------------------------------------------
        div = float(div_pct) / 100

        i_bal = []

        for i_alt in range(len(parsed_line["alt_len"])):

            alt_len = parsed_line["alt_len"][i_alt]

            if is_balanced(parsed_line["ref_len"], alt_len):
                
                # Incrementing alt index by 1 to start counting at 1
                # + skip ref walk at index 0 in bubble walks during annot
                i_bal.append(str(i_alt + 1))
        
        #---------------------------------------------------
        # Output balanced SV
        #---------------------------------------------------
        if len(i_bal) > 0:

            vcf_array = line.rstrip().split("\t")
            f_INFO = 7

            bal_INFO = "BL=" + ",".join(i_bal)
            new_INFO = ";".join([vcf_array[7], bal_INFO])

            vcf_array[f_INFO] = new_INFO

            new_line = "\t".join(vcf_array)
            out_vcf_balanced.write(new_line + "\n")

            count_BL_entries += 1


        # Group a1 lengths
        # a1_groups = group_a1(parsed_line["a1_len"], div)

        # sv_type, sv_len = define_type(parsed_line["a0_len"], a1_groups, div)

        #---------------------------------------------------
        # Save VCF lines of characterized DIVs for INV rescue
        #--------------------------------------------------- 
        # if "DIV" in sv_type:
        #     out_vcf_balanced.write(line)

        #---------------------------------------------------
        # Update stats
        #---------------------------------------------------
        # for i in range(len(sv_type)):
        #     count[sv_type[i]].append(sv_len[i])

        #---------------------------------------------------
        # Format output
        #---------------------------------------------------
        # out_bed.write(format_normal_output(parsed_line, sv_type) + "\n")
        # out_vcf.write(line)
        
# out_bed.close()
# out_intersect.close()
out_vcf_balanced.close()
# out_vcf.close()

#-----------------------------------------------------------
# Print statistics
#-----------------------------------------------------------
print(f"Balanced variants entries: {str(count_BL_entries)}")