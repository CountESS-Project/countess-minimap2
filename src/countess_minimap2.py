""" CountESS Minimap2 Plugin"""

import logging
from typing import Optional

import mappy  # type: ignore
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

logger = logging.getLogger(__name__)

VERSION = "0.0.15"

MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


class MiniMap2Plugin(PandasTransformSingleToDictPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    # XXX what is up with the CIGAR string not showing all variants?

    name = "MiniMap2 Plugin"
    description = "Finds variants using Minimap2"
    additional = "Note that the CIGAR string doesn't always show all variants."
    version = VERSION
    link = "https://github.com/CountESS-Project/countess-minimap2#readme"
    tags = ["bioinformatics"]

    FILE_TYPES = [("MMI", "*.mmi"), ("FASTA", "*.fa *.fasta *.fa.gz *.fasta.gz")]
    CHARACTER_SET = set(["A", "C", "G", "T"])

    column = ColumnChoiceParam("Input Column", "sequence")
    prefix = StringParam("Output Column Prefix", "mm")
    ref = FileParam("Ref FA / Ref MMI", file_types=FILE_TYPES)
    seq = StringCharacterSetParam("*OR* Ref Sequence", character_set=CHARACTER_SET)
    preset = ChoiceParam("Preset", "sr", choices=MM2_PRESET_CHOICES)
    min_length = IntegerParam("Minimum Match Length", 0)
    drop = BooleanParam("Drop Unmatched", False)
    location = BooleanParam("Output Location Columns", True)
    cigar = BooleanParam("Output Cigar String", False)
    cs = BooleanParam("Output CS String", False)
    hgvs = BooleanParam("Output HGVS", False)

    # XXX a shared-memory cache would make a lot of sense
    # here ...
    aligner = None

    def prepare(self, sources: list[str], row_limit: Optional[int] = None):
        if self.seq:
            self.aligner = mappy.Aligner(seq=self.seq.value, preset=self.preset.value)
        elif self.ref:
            self.aligner = mappy.Aligner(self.ref.value, preset=self.preset.value)
            # TODO check file load successful: self.aligner.seq_names is not None?
        else:
            self.aligner = None

    def output_dict(self, value, alignment):
        d = {}
        if self.location:
            d.update(
                {
                    self.prefix + "_ctg": alignment.ctg if alignment else None,
                    self.prefix + "_r_st": alignment.r_st if alignment else None,
                    self.prefix + "_r_en": alignment.r_en if alignment else None,
                    self.prefix + "_strand": alignment.strand if alignment else None,
                }
            )
        if self.cigar:
            d[self.prefix + "_cigar"] = alignment.cigar_str if alignment else None
        if self.cs:
            d[self.prefix + "_cs"] = alignment.cs if alignment else None
        if self.hgvs:
            if alignment:
                reference = self.aligner.seq(alignment.ctg)[alignment.r_st:alignment.r_en]
                d[self.prefix + "_hgvs_g"] = (
                    find_variant_string("g.", reference, value, offset=alignment.r_st)
                )
                d[self.prefix + "_hgvs_p"] = (
                    find_variant_string("p.", reference, value, offset=alignment.r_st)
                )
            else:
                d[self.prefix + "_hgvs_g"] = None
                d[self.prefix + "_hgvs_p"] = None

        return d

    def process_value(self, value: str):
        if not self.aligner:
            return None

        min_length = abs(self.min_length.value)

        x = self.aligner.map(value, cs=self.cs.value)
        # XXX only returns first match
        for z in x:
            if abs(z.r_en - z.r_st) >= min_length:
                return self.output_dict(value, z)

        if self.drop:
            return None
        else:
            return self.output_dict(value, None)
