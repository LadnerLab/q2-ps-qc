from q2_pepsirf.format_types import PepsirfContingencyTSVFormat, MutantReferenceFileFmt

def repScatters_tsv(
    ctx,
    source = None,
    pn_filepath= None,
    plot_log= False,
    zscore_filepath = None,
    col_sum_filepath = None,
    facet_charts = False,
    xy_threshold = None
):
    repScatters = ctx.get_action('ps-plot', 'repScatters')

    # import data into an artifact
    if zscore_filepath:
        zscore = ctx.make_artifact(type='FeatureTable[Zscore]',
                                    view=zscore_filepath,
                                view_type=PepsirfContingencyTSVFormat)
    else:
        zscore = None

    # import data into an artifact
    if col_sum_filepath:
        col_sum = ctx.make_artifact(type='FeatureTable[Normed]',
                                    view=col_sum_filepath,
                                view_type=PepsirfContingencyTSVFormat)
    else:
        col_sum = None

    repScatters_vis, = repScatters(
        source = source,
        pn_filepath = pn_filepath,
        plot_log = plot_log,
        zscore = zscore,
        col_sum = col_sum,
        facet_charts = facet_charts,
        xy_threshold = xy_threshold
    )

    return repScatters_vis

def mutantScatters_tsv(
    ctx,
    source,
    metadata,
    zscore_filepath,
    reference_file_filepath = None,
    peptide_header = "FeatureID",
    reference_header = "Reference",
    x_axis_header = 'Position',
    category_header = 'Category',
    label_header = 'Label',
    x_axis_label = 'Position',
    y_axis_label = 'Mean Zscore',
    min_wobble = -0.4,
    max_wobble = 0,
    wobble = False,
    scatter_only = True,
    scatter_boxplot = False,
    boxplot_only = False
):

    mutantScatters = ctx.get_action('ps-plot', 'mutantScatters')

    # import data into an artifact
    zscore = ctx.make_artifact(type='FeatureTable[Zscore]',
                            view=zscore_filepath,
                            view_type=PepsirfContingencyTSVFormat)
    
    # import data into an artifact
    if reference_file_filepath:
        reference_file = ctx.make_artifact(
            type = 'MutantReference',
            view = reference_file_filepath,
            view_type = MutantReferenceFileFmt
        )

    mutantScatters_vis, = mutantScatters(
        source = source,
        metadata = metadata,
        zscore = zscore,
        reference_file = reference_file,
        peptide_header = peptide_header,
        reference_header = reference_header,
        x_axis_header = x_axis_header,
        category_header = category_header,
        label_header = label_header,
        x_axis_label = x_axis_label,
        y_axis_label = y_axis_label,
        min_wobble = min_wobble,
        max_wobble = max_wobble,
        wobble = wobble,
        scatter_only = scatter_only,
        scatter_boxplot = scatter_boxplot,
        boxplot_only = boxplot_only
    )

    return mutantScatters_vis
