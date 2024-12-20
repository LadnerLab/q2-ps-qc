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
    MutantReference
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

plugin.pipelines.register_function(
    function = actions.compareCS,
    inputs = {},
    input_descriptions = None,
    parameters = {
        "metadata_file": Str,
        "codename_column": Str,
        "parent_codename_column": Str,
        "zscores_file": Str,
        "fasta_file": Str,
        "fullname_column": Str,
        "generate_epitope_data": Bool,
        "min_zscore": Float,
        "min_zscore_diff": Float,
        "max_zscore_diff": Float,
        "pep_seq_len": Int,
        "histogram_bins": Int,
        "min_epitope_size": Int,
        "data_output_dir": Str,
    },
    parameter_descriptions = {
        "metadata_file": "",
        "fullname_column": "Column denoting the full sequence name and alignment position for that peptide. Required for epitope data.",
        "codename_column": "Column denoting the codename of the sequence",
        "parent_codename_column": "Column denoting the codename of the cysteine version of the sequence"
            " for sequences in which cysteines were converted to serines",
        "zscores_file": "",
        "fasta_file": "Fasta file that includes cysteine version peptide sequences."
            " It should also include serine version peptides for epitope data",
        "generate_epitope_data": "Include this parameter to generate data for epitopes of similarly reactive peptides",
        "min_zscore": "",
        "min_zscore_diff": "",
        "max_zscore_diff": "",
        "pep_seq_len": "",
        "histogram_bins": "",
        "min_epitope_size": "",
        "data_output_dir": ""
    },
    outputs = [("reactivity_plot", Visualization), ("c_pos_histogram", Visualization)],
	output_descriptions = {
        "reactivity_plot": "",
        "c_pos_histogram": ""
    },
    name = "Compare C to S substitution reactivity",
    description = ""
)