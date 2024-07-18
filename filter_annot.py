#! /bin/python3

import sys

INPUT_ANNOT = sys.argv[1]
REF_ID = sys.argv[2]
MIN_COV = float(sys.argv[3])

def main():

    entries = []

    # Read input entries and save INV entries
    with open(INPUT_ANNOT, "r") as file:
        for line in file:

            if is_inv(line):

                entry = read_input(line)

                entry[3] = best_inv_annot(entry[3])

                if "na" in entry[3]:
                    continue

                # Filter on signal coverage
                if "INV" in entry[3] and signal_cov(entry[3]) < float(MIN_COV):
                    continue
                
                entries.append(entry)

    #==========================================#
    # FILTER ON NESTED BUBBLES WAS DEACTIVATED #
    #==========================================#

    # # Filter INV entries
    # i = 1
    # while i < len(entries):

    #     prev = entries[i-1]
    #     current = entries[i]

    #     if i == 0:
    #         i += 1

    #     elif is_nested(prev, current):

    #         to_remove = i - lowest_cov_entry(prev, current)

    #         # Remove entry with lowest signal coverage
    #         del entries[to_remove]

    #         # Update i
    #         i -= 1
        
    #     else:

    #         i += 1
    #==========================================#

    # Output filtered entries
    for e in entries:
        print(format_entry(e))


def read_input(line):

    f_ref = 0
    f_start = 1
    f_len = 3
    f_annot = 5

    line = line.rstrip().split("\t")

    entry = [line[f_ref], int(line[f_start]), int(line[f_start])+int(line[f_len])-1, line[f_annot]]

    return entry

def is_inv(line):

    if "INV" in line:
        return True
    
    return False

def best_inv_annot(annot):

    if ";" in annot:
        best_annot = "na:na:0.0"
        for a in list(annot.split(";")):

            if a == "DIV":
                continue
            
            if signal_cov(a) > signal_cov(best_annot):
                best_annot = a
        
        return best_annot

    else:
        return annot

def is_nested(prev_entry, current_entry):
    
    start = 1
    end = 2

    if current_entry[start] < prev_entry[end]:
        return True
    
    return False

def signal_cov(annot):
    return float(annot.split(":")[2].split(",")[0])

def lowest_cov_entry(prev_entry, current_entry):

    i_prev = 1
    i_current = 0

    start = 1
    end = 2
    annot = 3

    prev_cov_len = (prev_entry[end]-prev_entry[start]+1)*signal_cov(prev_entry[annot])
    current_cov_len = (current_entry[end]-current_entry[start]+1)*signal_cov(current_entry[annot])

    if prev_cov_len > current_cov_len:
        lowest = i_current
    
    elif current_cov_len > prev_cov_len:
        lowest = i_prev
    
    else:
        prev_entry_len = prev_entry[end] - prev_entry[start] + 1
        current_entry_len = current_entry[end] - current_entry[start] + 1

        if prev_entry_len > current_entry_len:
            lowest = i_current
        else:
            lowest = i_prev
    
    return lowest

def format_entry(entry):

    start = 1
    end = 2
    annot = 3

    output_line = "\t".join([REF_ID, str(entry[start]), str(entry[end]), entry[annot]])

    return output_line

main()