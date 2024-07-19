from q2_pepsirf.format_types import PepsirfContingencyTSVFormat


def filter_counts_matrix_tsv(
    ctx,
    input_matrix_filepath,
    count_thresh=None,
    max_zeros=None,
    drop_samp_out="filtered_out.tsv"
    ):
    
    filter_counts_matrix = ctx.get_action("q2-ps-qc", "filter_counts_matrix")

    input_matrix = ctx.make_artifact(
        type="FeatureTable[RawCounts]",
        view=input_matrix_filepath,
        view_type=PepsirfContingencyTSVFormat
    )

    filtered_matrix, = filter_counts_matrix(
        input_matrix=input_matrix,
        count_thresh=count_thresh,
        max_zeros=max_zeros,
        drop_samp_out=drop_samp_out
        )

    return filtered_matrix