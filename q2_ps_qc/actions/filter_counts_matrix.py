#!/usr/bin/env python

import pandas as pd
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat

def filter_counts_matrix(
    input_matrix: PepsirfContingencyTSVFormat,
    count_thresh: float=None,
    max_zeros: int=None,
    drop_samp_out: str="filtered_out.tsv"
    ) -> PepsirfContingencyTSVFormat:

    filtered_matrix_filepath = PepsirfContingencyTSVFormat()
    
    table = input_matrix.view(pd.DataFrame).transpose()

    row_len = len(table.index.values)

    # set default values if not specified
    if count_thresh is None:
        count_thresh = 2 * row_len
    if max_zeros is None:
        max_zeros = round(0.25 * row_len)

    filtered_samp = list()
    dropped_samp = list()
    # filter cols
    for col in table.columns:
        col_series = table[col]
        num_zeros = col_series.value_counts()[0]
        count = col_series.sum()
        if num_zeros > max_zeros or count < count_thresh:
            dropped_samp.append((col, count, num_zeros))
        else:
            filtered_samp.append(col)

    assert len(filtered_samp) > 0, \
        "None of the samples meet with the given threshold(s)"

    table[filtered_samp].to_csv(str(filtered_matrix_filepath), sep="\t")
    pd.DataFrame(dropped_samp, columns=['Sample name', 'Count', 'NumZeros']).to_csv(drop_samp_out, sep='\t', index=False)
    
    return filtered_matrix_filepath