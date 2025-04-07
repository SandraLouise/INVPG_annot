#! /bin/python3
from typing import Any
from pathlib import Path


def read_input(
    line: str
) -> list[str | int]:
    """Transfors ma bed line in a list, extracting meaningful information

    Parameters
    ----------
    line : str
        the bed line

    Returns
    -------
    list[str | int]
        formats the entry as a list
    """

    f_ref: int = 0
    f_start: int = 1
    f_len: int = 3
    f_annot: int = 5

    line: list[str] = line.rstrip().split("\t")

    return [line[f_ref], int(line[f_start]), int(
        line[f_start])+int(line[f_len])-1, line[f_annot]]


def is_inv(
    line: str,
) -> bool:
    """Analyses the annotation tag from a line of the bed file

    Parameters
    ----------
    line : str
        text describing the variation

    Returns
    -------
    bool
        if the inversion is flagged as an inversion
    """
    return "INV" in line


def best_inv_annot(
    annot: str,
) -> str:
    """Give a list of annotations, returns the best annotation (maximizing coverage)

    Parameters
    ----------
    annot : str
        annotations separated by semicolumns

    Returns
    -------
    str
        the descriptor of the best annotation
    """
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
    prev_entry: list[Any],
    current_entry: list[Any],
) -> bool:
    """Searches if the bubble is within a bubble or not

    Parameters
    ----------
    prev_entry : list[Any]
        descriptor of the previous entry
    current_entry : list[Any]
        descriptor of the current entry

    Returns
    -------
    bool
        information about the bubble being nested
    """
    start: int = 1
    end: int = 2

    return current_entry[start] < prev_entry[end]


def signal_cov(
    annot: str
) -> float:
    """Extracts coverage inforamtion

    Parameters
    ----------
    annot : str
        string describing the annotation

    Returns
    -------
    float
        coverage value
    """
    return float(annot.split(":")[2].split(",")[0])


def lowest_cov_entry(
    prev_entry: list[Any],
    current_entry: list[Any],
) -> int:
    """Searches for the input with the lowest coverage

    Parameters
    ----------
    prev_entry : list[Any]
        descriptor of the previous entry
    current_entry : list[Any]
        descriptor of the current entry

    Returns
    -------
    int
        id of the lowest input
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
    """Concatenates the input in a specific format

    Parameters
    ----------
    entry : list[str  |  int]
        File input
    reference_name : str
        Name of the path

    Returns
    -------
    str
        A formated string
    """
    start: int = 1
    end: int = 2
    annot: int = 3

    return "\t".join(
        [reference_name, str(entry[start]), str(entry[end]), entry[annot]]
    )


def filterannot(
    input_annotation_file: str,
    reference_name: str,
    minimum_coverage: float,
    output_prefix: str,
    timestamp: str,
) -> None:
    """Filters annotations of inversions obtained from the annotation step.

    Parameters
    ----------
    input_annotation_file : str
        Path to input bed file
    reference_name : str
        Name of the reference path
    minimum_coverage : float
        Minimum coverage of inversion signal. 
    """
    # Path for temporary files
    if '/' not in output_prefix:
        temp_folder = './'
    else:
        temp_folder = '/'.join(
            [x for x in output_prefix.split('/')][:-1]
        ) + '/'

    Path(temp_folder).mkdir(parents=True, exist_ok=True)

    # Read input entries and save INV entries
    with open(f"{temp_folder}{timestamp}.filtered.bed", "w", encoding='utf-8') as output_file:
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
