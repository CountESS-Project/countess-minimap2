""" CountESS Minimap2 Plugin"""

import re
from typing import Mapping, Optional

import mappy  # type: ignore
import pandas as pd
from countess.core.logger import Logger
from countess.core.parameters import (
    BooleanParam,
    ChoiceParam,
    ColumnChoiceParam,
    FileParam,
    IntegerParam,
    StringCharacterSetParam,
    StringParam,
)
from countess.core.plugins import PandasTransformSingleToDictPlugin

VERSION = "0.0.9"

CS_STRING_RE = r"(=[ACTGTN]+|:[0-9]+|(?:\*[ACGTN][ACGTN])+|\+[ACGTN]+|-[ACGTN]+)"
MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


def cs_to_hgvs(cs_string: str, offset: int = 1) -> str:
    """Turn the Minimap2 "difference string" into a HGVS string"""

    hgvs_ops = []
    for op in re.findall(CS_STRING_RE, cs_string.upper()):
        if op[0] == ":":
            offset += int(op[1:])
        elif op[0] == "=":
            offset += len(op) - 1
        elif op[0] == "*":
            # regex can match multiple operations like "*AT*AT*GC"
            if len(op) > 3:
                hgvs_ops.append(f"{offset}delins{op[2::3]}")
            else:
                hgvs_ops.append(f"{offset}{op[1]}>{op[2]}")
            offset += len(op) // 3
        elif op[0] == "+":
            hgvs_ops.append(f"{offset}_{offset+1}ins{op[1]}")
            offset += 1
        elif op[0] == "-":
            hgvs_ops.append(f"{offset}del")
            offset += 1
    if len(hgvs_ops) == 0:
        return "g.="
    elif len(hgvs_ops) == 1:
        return "g." + hgvs_ops[0]
    else:
        return "g.[" + ";".join(hgvs_ops) + "]"


class MiniMap2Plugin(PandasTransformSingleToDictPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    # XXX what is up with the CIGAR string not showing all variants?

    name = "MiniMap2 Plugin"
    description = "Finds variants using Minimap2"
    additional = "Note that the CIGAR string doesn't always show all variants."
    version = VERSION
    link = "https://github.com/CountESS-Project/countess-minimap2#readme"

    FILE_TYPES = [("MMI", "*.mmi"), ("FASTA", "*.fa *.fasta *.fa.gz *.fasta.gz")]
    CHARACTER_SET = set(["A", "C", "G", "T"])

    parameters = {
        "column": ColumnChoiceParam("Input Column", "sequence"),
        "prefix": StringParam("Output Column Prefix", "mm"),
        "ref": FileParam("Ref FA / Ref MMI", file_types=FILE_TYPES),
        "seq": StringCharacterSetParam("*OR* Ref Sequence", character_set=CHARACTER_SET),
        "preset": ChoiceParam("Preset", "sr", choices=MM2_PRESET_CHOICES),
        "min_length": IntegerParam("Minimum Match Length", 0),
        "drop": BooleanParam("Drop Unmatched", False),
    }

    aligner = None

    def prepare(self, sources: list[str], row_limit: Optional[int]):
        if self.parameters["seq"].value:
            self.aligner = mappy.Aligner(seq=self.parameters["seq"].value, preset=self.parameters["preset"].value)
        elif self.parameters["ref"].value:
            self.aligner = mappy.Aligner(self.parameters["ref"].value, preset=self.parameters["preset"].value)
        else:
            self.aligner = None

    def process_value(self, value: str, logger: Logger):
        assert isinstance(self.parameters["column"], ColumnChoiceParam)
        if not self.aligner:
            return None

        prefix = self.parameters["prefix"].value
        min_length = self.parameters["min_length"].value

        # XXX only returns first match
        x = self.aligner.map(value, cs=True)
        for z in x:
            if z.r_en - z.r_st >= min_length:
                return {
                    prefix + "_ctg": z.ctg,
                    prefix + "_r_st": z.r_st,
                    prefix + "_r_en": z.r_en,
                    prefix + "_strand": z.strand,
                    prefix + "_cigar": z.cigar_str,
                    prefix + "_cs": z.cs,
                    prefix + "_hgvs": cs_to_hgvs(z.cs, z.r_st + 1),
                }

        if self.parameters["drop"].value:
            return None
        else:
            return {
                prefix + "_ctg": None,
                prefix + "_r_st": 0,
                prefix + "_r_en": 0,
                prefix + "_strand": 0,
                prefix + "_cigar": "",
                prefix + "_cs": "",
                prefix + "_hgvs": "",
            }
