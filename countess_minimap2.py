""" CountESS Minimap2 Plugin"""

import re
from typing import Any, Mapping, Optional

import pandas as pd

import mappy  # type: ignore

from countess.core.logger import Logger
from countess.core.parameters import (
    BooleanParam,
    ChoiceParam,
    ColumnOrIndexChoiceParam,
    IntegerParam,
    StringParam,
    FileParam,
)
from countess.core.plugins import PandasTransformPlugin

VERSION = '0.0.5'

CS_STRING_RE = r"(=[ACTGTN]+|:[0-9]+|(?:\*[ACGTN][ACGTN])+|\+[ACGTN]+|-[ACGTN]+)"
MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


def cs_to_hgvs(cs_string: str, offset: int=1) -> str:
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
        return "g.{offset}="
    elif len(hgvs_ops) == 1:
        return "g." + hgvs_ops[0]
    else:
        return "g.[" + ";".join(hgvs_ops) + "]"

class MiniMap2Plugin(PandasTransformPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    # XXX what is up with the CIGAR string not showing all variants?

    name = "MiniMap2 Plugin"
    description = """
        Finds variants using Minimap2.  Note that the CIGAR string doesn't always
        show all variants.
    """
    version = VERSION
    link = "https://github.com/nickzoic/countess-minimap2#readme"

    FILE_TYPES = [("MMI", "*.mmi"), ("FASTA", "*.fa(sta)?")]

    parameters = {
        "column": ColumnOrIndexChoiceParam("Input Column", "sequence"),
        "prefix": StringParam("Output Column Prefix", "mm"),
        "ref": FileParam("Ref FA / Ref MMI", file_types = FILE_TYPES),
        "preset": ChoiceParam("Preset", "sr", choices=MM2_PRESET_CHOICES),
        "min_length": IntegerParam("Minimum Match Length", 0),
        "drop": BooleanParam("Drop Unmatched", False),
    }

    def run_df(self, df: pd.DataFrame, logger: Logger) -> pd.DataFrame:

        assert isinstance(self.parameters['column'], ColumnOrIndexChoiceParam)

        prefix = self.parameters["prefix"].value

        aligner = mappy.Aligner(
            self.parameters["ref"].value, preset=self.parameters["preset"].value
        )

        if not aligner:
            logger.error("ERROR: failed to load/build index file")
            return pd.DataFrame()

        df = df.copy()

        if self.parameters["column"].is_index():
            df["__index"] = df.index
            column_name = "__index"
        else:
            column_name = self.parameters["column"].value

        prefix = self.parameters["prefix"].value
        min_length = self.parameters["min_length"].value

        def process(row: Mapping[str, Any]) -> tuple[Optional[str], int, int, int, Optional[str], Optional[str], Optional[str]]:
            # XXX only returns first match
            x = aligner.map(row[column_name], cs=True)
            for z in x:
                if z.r_en - z.r_st >= min_length:
                    return (z.ctg, z.r_st, z.r_en, z.strand, z.cigar_str, z.cs, \
                            cs_to_hgvs(z.cs, z.r_st+1))
            return (None, 0, 0, 0, None, None, None)

        column_names = [
            prefix + "_ctg",
            prefix + "_r_st",
            prefix + "_r_en",
            prefix + "_r_strand",
            prefix + "_cigar",
            prefix + "_cs",
            prefix + "_hgvs",
        ]
        df[column_names] = df.apply(process, axis=1, result_type="expand")

        if self.parameters["column"].is_index():
            df.drop(columns=["__index"])

        if self.parameters["drop"].value:
            df = df.query("mm_ctg.notnull()")

        return df
