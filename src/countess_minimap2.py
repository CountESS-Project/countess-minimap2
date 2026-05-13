""" CountESS Minimap2 Plugin"""

import logging
from typing import Optional, Dict, Any

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
from countess.core.plugins import DuckdbTransformPlugin
from countess.utils.variant import find_variant_string

logger = logging.getLogger(__name__)

VERSION = "0.1.1"

MM2_PRESET_CHOICES = ["sr", "map-pb", "map-ont", "asm5", "asm10", "splice"]


class MiniMap2Plugin(DuckdbTransformPlugin):
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

    def prepare(self, *_):
        if self.seq:
            self.aligner = mappy.Aligner(seq=self.seq.value, preset=self.preset.value)
        elif self.ref:
            self.aligner = mappy.Aligner(self.ref.value, preset=self.preset.value)
            # TODO check file load successful: self.aligner.seq_names is not None?
        else:
            self.aligner = None

    def add_fields(self):
        cols = {}
        if self.location:
            cols.update({
                self.prefix + "_ctg": str,
                self.prefix + "_r_st": int,
                self.prefix + "_r_en": int,
                self.prefix + "_strand": int,
            })
        if self.cigar:
            cols[self.prefix + "_cigar"] = str
        if self.cs:
            cols[self.prefix + "_cs"] = str
        if self.hgvs:
            cols[self.prefix + "_hgvs_g"] = str
            cols[self.prefix + "_hgvs_p"] = str
        return cols

    def transform(self, data: dict[str, Any]) -> Optional[Dict[str, Any]]:
        value = data[self.column.value]
        min_length = int(self.min_length.value)
        alignments = list(self.aligner.map(value, cs=self.cs.value))
        if not alignments:
            return None
        for alignment in alignments:
            if abs(alignment.r_en - alignment.r_st) >= min_length:
                if self.location:
                    data.update({
                        self.prefix + "_ctg": alignment.ctg,
                        self.prefix + "_r_st": alignment.r_st,
                        self.prefix + "_r_en": alignment.r_en,
                        self.prefix + "_strand": alignment.strand,
                    })
                if self.cigar:
                    data[self.prefix + "_cigar"] = alignment.cigar_str
                if self.cs:
                    data[self.prefix + "_cs"] = alignment.cs
                if self.hgvs:
                    reference = self.seq.value or self.aligner.seq(alignment.ctg)[alignment.r_st:alignment.r_en]
                    data[self.prefix + "_hgvs_g"] = find_variant_string("g.", reference, value, offset=alignment.r_st)
                    data[self.prefix + "_hgvs_p"] = find_variant_string("p.", reference, value, offset=alignment.r_st)
                return data
        return None
