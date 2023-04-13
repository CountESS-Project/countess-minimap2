import dask.dataframe as dd
import pandas as pd

import mappy

from countess.core.logger import Logger
from countess.core.parameters import ChoiceParam, ColumnChoiceParam, IntegerParam, StringParam, FileParam
from countess.core.plugins import DaskTransformPlugin
from countess.utils.dask import empty_dask_dataframe


class MiniMap2Plugin(DaskTransformPlugin):
    """Turns a DNA sequence into a HGVS variant code"""

    name = "MiniMap2 Plugin"
    description = "Finds variants using Minimap2"
    version = '0.0.1'

    parameters = {
        "column": ColumnChoiceParam("Input Column"),
        "prefix": StringParam("Output Column Prefix"),
        "ref": FileParam("Ref FA / Ref MMI"),
        "preset": ChoiceParam("Preset", choices=['sr', 'map-pb', 'map-ont', 'asm5', 'asm10', 'splice']),
        # XXX etc
    }

    def run_dask(self, df: pd.DataFrame|dd.DataFrame, logger: Logger) -> pd.DataFrame|dd.DataFrame:

        prefix = self.parameters['prefix'].value

        aligner = mappy.Aligner(self.parameters['ref'].value, preset=self.parameters['preset'].value)
        if not aligner:
            logger.error("ERROR: failed to load/build index file")
            return empty_dask_dataframe()

        column = self.parameters['column'].value
        prefix = self.parameters['prefix'].value

        def process(row):
            # XXX only returns first match
            x = aligner.map(row[column])
            for z in x:
                return [ z.ctg, z.r_st, z.r_en, z.cigar_str ]
            return [ None, 0, 0, None ]

        df = df.copy()
        df[[prefix + "_ctg", prefix + "_r_st", prefix + "_r_en", prefix + "_cigar"]] = df.apply(process, axis=1, result_type='expand', meta={0: str, 1: int, 2: int, 3: str})
        return df
