#! /bin/python3
from typing import Any
from os import path


def read_input(
    line: str
) -> list[str | int]:
    """_summary_

    Parameters
    ----------
    line : str
        _description_

    Returns
    -------
    list[str | int]
        _description_
    """

    f_ref: int = 0
    f_start: int = 1
    f_len: int = 3
    f_annot: int = 5

    line = line.rstrip().split("\t")

    entry = [line[f_ref], int(line[f_start]), int(
        line[f_start])+int(line[f_len])-1, line[f_annot]]

    return entry


def is_inv(
    line
) -> bool:
    """_summary_

    Parameters
    ----------
    line : _type_
        _description_

    Returns
    -------
    bool
        _description_
    """
    return "INV" in line


def best_inv_annot(
    annot
):

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


def is_nested(
    prev_entry,
    current_entry,
) -> bool:

    start: int = 1
    end: int = 2

    return current_entry[start] < prev_entry[end]


def signal_cov(
    annot: str
) -> float:
    """_summary_

    Parameters
    ----------
    annot : str
        _description_

    Returns
    -------
    float
        _description_
    """
    return float(annot.split(":")[2].split(",")[0])


def lowest_cov_entry(
    prev_entry: list[Any],
    current_entry: list[Any],
) -> int:
    """_summary_

    Parameters
    ----------
    prev_entry : list[Any]
        _description_
    current_entry : list[Any]
        _description_

    Returns
    -------
    int
        _description_
    """
    i_prev: int = 1
    i_current: int = 0

    start: int = 1
    end: int = 2
    annot: int = 3

    prev_cov_len = (prev_entry[end]-prev_entry[start]+1) * \
        signal_cov(prev_entry[annot])
    current_cov_len = (
        current_entry[end]-current_entry[start]+1)*signal_cov(current_entry[annot])

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


def format_entry(
    entry: list[str | int],
    reference_name: str,
) -> str:
    """_summary_

    Parameters
    ----------
    entry : list[str  |  int]
        _description_
    reference_name : str
        _description_

    Returns
    -------
    str
        _description_
    """
    start: int = 1
    end: int = 2
    annot: int = 3

    return "\t".join(
        [reference_name, str(entry[start]), str(entry[end]), entry[annot]]
    )


def filterannot_main(
    input_annotation_file: str,
    reference_name: str,
    minimum_coverage: float,
) -> None:
    """_summary_

    Parameters
    ----------
    input_annotation_file : str
        _description_
    reference_name : str
        _description_
    minimum_coverage : float
        _description_
    """
    # Read input entries and save INV entries
    with open(f"{path.splitext(input_annotation_file)[0]}.filtered.bed", "w", encoding='utf-8') as output_file:
        with open(input_annotation_file, "r", encoding='utf-8') as input_file:
            for line in input_file:
                if is_inv(line):

                    entry: list[str | int] = read_input(line)

                    entry[3] = best_inv_annot(entry[3])

                    if "na" in entry[3]:
                        continue

                    # Filter on signal coverage
                    if "INV" in entry[3] and signal_cov(entry[3]) < float(minimum_coverage):
                        continue

                    output_file.write(
                        format_entry(
                            entry=entry,
                            reference_name=reference_name
                        )
                    )

        # ==========================================#
        # FILTER ON NESTED BUBBLES WAS DEACTIVATED #
        # ==========================================#

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
        # ==========================================#
