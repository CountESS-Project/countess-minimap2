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

VERSION = "0.0.12"

CS_STRING_RE = r"(=[ACTGTN]+|:[0-9]+|(?:\*[ACGTN][ACGTN])+|\+[ACGTN]+|-[ACGTN]+)"
MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


def cs_to_hgvs(cs_string: str, ctg: str = "", offset: int = 1) -> str:
    """Turn the Minimap2 "difference string" into a HGVS string

    The only reference I have for the "difference string" is from the `minimap2`
    manual page https://lh3.github.io/minimap2/minimap2.html#10 which contains
    the following description:

       The cs tag encodes difference sequences in the short form or the entire
       query AND reference sequences in the  long  form.  It consists of a
       series of operations:

       ┌───┬────────────────────────────┬─────────────────────────────────┐
       │Op │           Regex            │           Description           │
       ├───┼────────────────────────────┼─────────────────────────────────┤
       │ = │ [ACGTN]+                   │ Identical sequence (long form)  │
       │ : │ [0-9]+                     │ Identical sequence length       │
       │ * │ [acgtn][acgtn]             │ Substitution: ref to query      │
       │ + │ [acgtn]+                   │ Insertion to the reference      │
       │ - │ [acgtn]+                   │ Deletion from the reference     │
       │ ~ │ [acgtn]{2}[0-9]+[acgtn]{2} │ Intron length and splice signal │
       └───┴────────────────────────────┴─────────────────────────────────┘

    >>> cs_to_hgvs(":10*at:10")
    'g.11A>T'

    >>> cs_to_hgvs(":10-c:10")
    'g.11del'

    >>> cs_to_hgvs(":10-ttt:10")
    'g.11_13del'

    >>> cs_to_hgvs(":10+gac:10")
    'g.10_11insGAC'

    >>> cs_to_hgvs(":10*at*at*gc*cg")
    'g.11_14delinsTTCG'
    """

    # XXX doesn't support '~'.

    hgvs_ops = []
    prefix = "g."
    if ctg and ctg != 'N/A':
        prefix = ctg + ":" + prefix

    for op in re.findall(CS_STRING_RE, cs_string.upper()):
        if op[0] == ":":
            offset += int(op[1:])
        elif op[0] == "=":
            offset += len(op) - 1
        elif op[0] == "*":
            # regex can match multiple operations like "*AT*AT*GC"
            if len(op) > 3:
                hgvs_ops.append(f"{offset}_{offset+len(op)//3-1}delins{op[2::3]}")
            else:
                hgvs_ops.append(f"{offset}{op[1]}>{op[2]}")
            offset += len(op) // 3
        elif op[0] == "+":
            hgvs_ops.append(f"{offset-1}_{offset}ins{op[1:]}")
            offset += 1
        elif op[0] == "-":
            if len(op) > 2:
                hgvs_ops.append(f"{offset}_{offset+len(op)-2}del")
            else:
                hgvs_ops.append(f"{offset}del")
            offset += len(op)
    if len(hgvs_ops) == 0:
        return prefix + "="
    elif len(hgvs_ops) == 1:
        return prefix + hgvs_ops[0]
    else:
        return prefix + "[" + ";".join(hgvs_ops) + "]"


class MiniMap2Plugin(PandasTransformSingleToDictPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    # XXX what is up with the CIGAR string not showing all variants?

    name = "MiniMap2 Plugin"
    description = "Finds variants using Minimap2"
    additional = "Note that the CIGAR string doesn't always show all variants."
    version = VERSION
    link = "https://github.com/CountESS-Project/countess-minimap2#readme"
    tags = ['bioinformatics']

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
        "location": BooleanParam("Output Location Columns", True),
        "cigar": BooleanParam("Output Cigar String", True),
        "cs": BooleanParam("Output CS String", False),
        "hgvs": BooleanParam("Output HGVS String", False),
    }

    # XXX a shared-memory cache would make a lot of sense
    # here ...
    aligner = None

    def prepare(self, sources: list[str], row_limit: Optional[int]):
        if self.parameters["seq"].value:
            self.aligner = mappy.Aligner(seq=self.parameters["seq"].value, preset=self.parameters["preset"].value)
        elif self.parameters["ref"].value:
            self.aligner = mappy.Aligner(self.parameters["ref"].value, preset=self.parameters["preset"].value)
        else:
            self.aligner = None

    def output_dict(self, alignment):
        d = {}
        prefix = self.parameters["prefix"].value
        if self.parameters["location"].value:
            d.update(
                {
                    prefix + "_ctg": alignment.ctg if alignment else None,
                    prefix + "_r_st": alignment.r_st if alignment else None,
                    prefix + "_r_en": alignment.r_en if alignment else None,
                    prefix + "_strand": alignment.strand if alignment else None,
                }
            )
        if self.parameters["cigar"].value:
            d[prefix + "_cigar"] = alignment.cigar_str if alignment else None
        if self.parameters["cs"].value:
            d[prefix + "_cs"] = alignment.cs if alignment else None
        if self.parameters["hgvs"].value:
            d[prefix + "_hgvs"] = cs_to_hgvs(alignment.cs, alignment.ctg, alignment.r_st + 1) if alignment else None
        return d

    def process_value(self, value: str, logger: Logger):
        assert isinstance(self.parameters["column"], ColumnChoiceParam)
        if not self.aligner:
            return None

        prefix = self.parameters["prefix"].value
        min_length = self.parameters["min_length"].value

        # XXX only returns first match
        calculate_cs = self.parameters["cs"].value or self.parameters["hgvs"].value
        x = self.aligner.map(value, cs=calculate_cs)
        for z in x:
            if z.r_en - z.r_st >= min_length:
                return self.output_dict(z)

        if self.parameters["drop"].value:
            return None
        else:
            return self.output_dict(None)
