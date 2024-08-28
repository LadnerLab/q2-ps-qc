#!/usr/bin/env python
import importlib
import q2_ps_qc
import q2_ps_qc.actions as actions

from qiime2.plugin import (
    Plugin, SemanticType, model,
    Int, Range, MetadataColumn,
    Categorical, Str, List,
    Visualization, Metadata, Bool, 
    Float
)
from q2_pepsirf.format_types import (
    Normed, Zscore, InfoSumOfProbes,
    PairwiseEnrichment, InfoSNPN, ProteinAlignment,
    MutantReference, RawCounts 
)
from q2_types.feature_table import FeatureTable, BIOMV210DirFmt


# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin(
    "q2-ps-qc", version = q2_ps_qc.__version__,
    website = "https://github.com/LadnerLab/q2-ps-qc",
    description = "Qiime2 Plug-in for the creation of correlation visualizations from PepSIRF outputs."
)

plugin.pipelines.register_function(
    function = actions.generate_corr_matrix,
    inputs = {},
    input_descriptions = None,
    parameters = {
        "data": Str,
        "samples": Str,
        "log_normalization": Bool,
        "correlation_threshold": Float
    },
    parameter_descriptions = {
		"data": "Name of input file.",
        "samples": "The name of the tab-delimited file containing sample"
            " information, denoting which samples, in the input matrices, are"
            " replicates. This file must be tab-delimited with each line"
            " containing a set of replicates.",
        "log_normalization": "Run a log normalization on each of the sets of"
            " scores before running a correlation test on them.",
        "correlation_threshold": "Set a threshold value; anything below the"
            " value will be considered a bad correlation score, and anything"
            " above will be considered a good correlation score.",
    },
    outputs = [("bad_output", Visualization), ("good_output", Visualization)],
	output_descriptions = {
        "bad_output": "File name for bad correlation visualization",
        "good_output": "File name for good correlation visualization"
    },
    name = "Generate Correlation Matrix",
    description = "Finds all replicate pairs that have poor correlation and"
        " creates a .qsv file that allows the user to visualize them in a"
        " scatter plot."
)


shared_parameters = {   "count_thresh": Float,
                        "max_zeros": Int,
                        "drop_samp_out": Str
}
shared_descriptions = { "count_thresh": "Minimum sequence count to not be filtered out. If None is provided,"
                            " default is 2x the total number of unique peptides.",
                        "max_zeros": "Maximum number of zero counts a sequence needs to not be filtered out. If"
                            " None is provided, default is 25% of the total number of unique peptides",
                        "drop_samp_out": "Filepath to output sample names that are filtered out of matrix."
}

plugin.methods.register_function(
    function = actions.filter_counts_matrix,
    inputs = {
        "input_matrix": FeatureTable[RawCounts]
    },
    input_descriptions = {
        "input_matrix": "FeatureTable containing raw PepSIRF counts matrix for filtering."
    },
    parameters = {
        **shared_parameters
    },
    parameter_descriptions = {
        **shared_descriptions
    },
    outputs = [("filtered_matrix", FeatureTable[RawCounts])],
    output_descriptions = {
        "filtered_matrix": "File name for filtered counts output"
    },
    name = "Filter Counts Matrix",
    description = "Takes in a raw counts matrix and removes sequences below a specific raw read count threshold "
        "and filters for a number of 0 counts."
)

plugin.pipelines.register_function(
    function = actions.filter_counts_matrix_tsv,
    inputs = {},
    input_descriptions = None,
    parameters = {
        "input_matrix_filepath": Str,
        **shared_parameters
    },
    parameter_descriptions = {
        "input_matrix_filepath": "Filepath to .tsv containing raw PepSIRF counts matrix for filtering.",
        **shared_descriptions
    },
    outputs = [("filtered_matrix", FeatureTable[RawCounts])],
    output_descriptions = {
        "filtered_matrix": "File name for filtered counts output"
    },
    name = "Filter Counts Matrix TSV Pipeline",
    description = "Pipeline that converts .tsv files to .qza files and then runs"
        " Filter Counts Matrix.."
)

