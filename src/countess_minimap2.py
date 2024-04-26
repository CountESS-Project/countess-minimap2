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
from countess.utils.variant import find_variant_string

VERSION = "0.0.13"

MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


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
        "cigar": BooleanParam("Output Cigar String", False),
        "cs": BooleanParam("Output CS String", False),
        "hgvs": BooleanParam("Output HGVS", False),
    }

    # XXX a shared-memory cache would make a lot of sense
    # here ...
    aligner = None

    def prepare(self, sources: list[str], row_limit: Optional[int]):
        if self.parameters["seq"].value:
            self.aligner = mappy.Aligner(seq=self.parameters["seq"].value, preset=self.parameters["preset"].value)
        elif self.parameters["ref"].value:
            self.aligner = mappy.Aligner(self.parameters["ref"].value, preset=self.parameters["preset"].value)
            # TODO check file load successful: self.aligner.seq_names is not None?
        else:
            self.aligner = None

    def output_dict(self, value, alignment):
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
            d[prefix + "_hgvs_g"] = find_variant_string("g.", alignment.tseq, value, alignment.r_st) if alignment else None
            d[prefix + "_hgvs_p"] = find_variant_string("p.", alignment.tseq, value, alignment.r_st) if alignment else None
        return d

    def process_value(self, value: str, logger: Logger):
        assert isinstance(self.parameters["column"], ColumnChoiceParam)
        if not self.aligner:
            return None

        prefix = self.parameters["prefix"].value
        min_length = abs(self.parameters["min_length"].value)

        # XXX only returns first match
        calculate_cs = self.parameters["cs"].value
        calculate_tseq = self.parameters["hgvs"].value
        x = self.aligner.map(value, cs=calculate_cs, tseq=calculate_tseq)
        for z in x:
            if abs(z.r_en - z.r_st) >= min_length:
                return self.output_dict(value, z)

        if self.parameters["drop"].value:
            return None
        else:
            return self.output_dict(value, None)
