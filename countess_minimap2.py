import re

import pandas as pd

import mappy

from countess.core.logger import Logger
from countess.core.parameters import BooleanParam, ChoiceParam, ColumnOrIndexChoiceParam, IntegerParam, StringParam, FileParam
from countess.core.plugins import PandasTransformPlugin

MM2_PRESET_CHOICES = ['sr', 'map-pb', 'map-ont', 'asm5', 'asm  10', 'splice']

class MiniMap2Plugin(PandasTransformPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    name = "MiniMap2 Plugin"
    description = "Finds variants using Minimap2"
    version = '0.0.1'


    parameters = {
        "column": ColumnOrIndexChoiceParam("Input Column"),
        "prefix": StringParam("Output Column Prefix"),
        "ref": FileParam("Ref FA / Ref MMI"),
        "preset": ChoiceParam("Preset", choices=MM2_PRESET_CHOICES),
        "min_length": IntegerParam("Minimum Match Length", 0),
        "drop": BooleanParam("Drop Unmatched", False),
    }

    def run_df(self, df: pd.DataFrame, logger: Logger) -> pd.DataFrame:

        prefix = self.parameters['prefix'].value

        aligner = mappy.Aligner(
            self.parameters['ref'].value,
            preset=self.parameters['preset'].value
        )

        if not aligner:
            logger.error("ERROR: failed to load/build index file")
            return pd.DataFrame()

        df = df.copy()

        if self.parameters["column"].is_index():
            df['__index'] = df.index
            column_name = '__index'
        else:
            column_name = self.parameters['column'].value

        prefix = self.parameters['prefix'].value
        min_length = self.parameters["min_length"].value

        def process(row):
            # XXX only returns first match
            x = aligner.map(row[column_name], cs=True)
            for z in x:
                if z.r_en - z.r_st >= min_length:
                    return [ z.ctg, z.r_st, z.r_en, z.cigar_str, z.cs ]
            return [ None, 0, 0, None, None ]

        column_names = [
            prefix + "_ctg",
            prefix + "_r_st",
            prefix + "_r_en",
            prefix + "_cigar",
            prefix + "_cs"
        ]
        df[column_names] = df.apply(process, axis=1, result_type='expand')

        if self.parameters["column"].is_index():
            df.drop(columns=['__index'])

        if self.parameters["drop"].value:
            df = df.query("mm_ctg.notnull()")

        return df
