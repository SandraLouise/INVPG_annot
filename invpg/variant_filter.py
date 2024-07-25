#! /bin/python3
from typing import Any
from os import path

# ===========================================================
# Functions
# ===========================================================


def parse_vcf_line(line: str) -> dict[str, Any]:
    """_summary_

    Parameters
    ----------
    line : str
        _description_

    Returns
    -------
    dict[str, Any]
        _description_
    """

    tab_parsed: list[str] = line.rstrip().split("\t")

    if "," in tab_parsed[4]:
        alt: list[str] = alt.split(",")
    else:
        alt: list[str] = [alt]
    ref: str = tab_parsed[3]

    return {
        "ref_id": tab_parsed[0],
        "ref_start": int(tab_parsed[1]),
        "ref_end": int(tab_parsed[1]) + len(ref) - 1,
        "ref_seq": ref,
        "alt_seq": alt,
        "ref_len": len(ref),
        "alt_len": [len(a) for a in alt],
        "lvl": int(tab_parsed[7].split(";LV=")[1].split(";")[0]) if ";LV=" in tab_parsed[7] else "NA"
    }


def is_balanced(
    ref_len: int,
    alt_len: int,
    div: int,
) -> bool:
    """_summary_

    Parameters
    ----------
    ref_len : int
        _description_
    alt_len : int
        _description_

    Returns
    -------
    bool
        _description_
    """
    max_a: int = max(ref_len, alt_len)
    min_a: int = min(ref_len, alt_len)

    return max_a - min_a <= (div * max_a)

# ===========================================================
# Main
# ===========================================================


def filter_main(
    in_vcf: str,
    div_pct: int
) -> str:
    """_summary_

    Parameters
    ----------
    in_vcf : str
        _description_
    div_pct : int
        _description_

    Returns
    -------
    str
        Path to the balanced SV VCF file.
    """
    count_BL_entries: int = 0
    div: float = float(div_pct) / 100
    min_len: int = 50
    f_INFO: int = 7

    with open(outVCF := f"{path.splitext(in_vcf)[0]}.balancedSV.vcf", 'w', encoding='utf-8') as out_vcf_balanced:
        with open(in_vcf, 'r', encoding='utf-8') as file:
            for line in file:

                if line.startswith("#"):
                    out_vcf_balanced.write(line)
                    continue

                parsed_line: dict[str, Any] = parse_vcf_line(line)

                if parsed_line["ref_len"] <= min_len and all([alt <= min_len for alt in parsed_line["alt_len"]]):
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
