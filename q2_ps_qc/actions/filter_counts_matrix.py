#!/usr/bin/env python

import pandas as pd
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat

def filter_counts_matrix(
    input_matrix: PepsirfContingencyTSVFormat,
    count_thresh: float=None,
    max_zeros: int=None
    ) -> PepsirfContingencyTSVFormat:

    filtered_matrix_filepath = PepsirfContingencyTSVFormat()
    
    table = input_matrix.view(pd.DataFrame).transpose()

    row_len = len(table.index.values)

    # set default values if not specified
    if count_thresh is None:
        count_thresh = 2 * row_len
    if max_zeros is None:
        max_zeros = round(0.25 * row_len)

    filtered_cols = list()
    # filter cols
    for col in table.columns:
        col_series = table[col]
        if col_series.value_counts()[0] <= max_zeros and col_series.sum() >= count_thresh:
            filtered_cols.append(col)

    table[filtered_cols].to_csv(str(filtered_matrix_filepath), sep="\t")
    
    return filtered_matrix_filepath