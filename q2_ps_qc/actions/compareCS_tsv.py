#!/usr/bin/env python

from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat
)

from qiime2.plugin import (
    Metadata
)

def compareCS_tsv(
        ctx,
        metadata_file,
        c_codename_column,
        s_codename_column,
        zscores,
        min_zscore
):
    compareCS = ctx.get_action("ps-q2-ps-qc", "compareCS")

    metadata_file = Metadata.load(metadata_file)

    # import data into an artifact
    zscores = ctx.make_artifact(
        type="FeatureTable[Zscore]",
        view=zscores,
        view_type=PepsirfContingencyTSVFormat
    )

    compareCS_vis, = compareCS(
        metadata_file = metadata_file,
        c_codename_column = c_codename_column,
        s_codename_column = s_codename_column,
        zscores = zscores,
        min_zscore = min_zscore
    )

    return compareCS_vis