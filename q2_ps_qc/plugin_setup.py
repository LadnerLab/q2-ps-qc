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
        "align_start_column": Str,
        "align_stop_column": Str,
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
        "metadata_file": "Tab delimited file with columns for the C-version peptide (--p-parent-codename-column),"
            " S-version peptide (--p-codename-column), and optionally (for eptitope-level analysis) the the full"
            " sequence name (--p-fullname-column), the alignment start position (--p-align-start-pos),"
            " and the alignment stop position (--p-align-stop-column).",
        "codename_column": "Column denoting the codename of the sequence",
        "parent_codename_column": "Column denoting the codename of the cysteine version of the sequence"
            " for sequences in which cysteines were converted to serines",
        "zscores_file": "TSV file containing z scores of the normalized"
            " read counts. Fist column header must be 'Sequence Name' as"
            " produced by pepsirf.",
        "fasta_file": "Fasta file that includes cysteine version peptide sequences."
            " It should also include serine version peptides for epitope data",
        "fullname_column": "Column denoting the full sequence name for that peptide. Required for epitope data generation.",
        "align_start_column": "Column denoting the start position in the full sequence for that peptide. Required for epitope data generation.",
        "align_stop_column": "Column denoting the stop position in the full sequence for that peptide. Required for epitope data generation.",
        "generate_epitope_data": "Include this parameter to generate data for epitopes of similarly reactive peptides",
        "min_zscore": "Minimun threshold that either version's z-score must meet to be included",
        "min_zscore_diff": "Minimum difference ratio between the two version's z-scores to be"
            " considered more reactive for a single version.",
        "max_zscore_diff": "Maximum difference ratio between the two version's z-scores to be"
            " considered similary reactive.",
        "pep_seq_len": "Length of each peptide sequence.",
        "histogram_bins": "Number of position bins to include in the histogram visualization.",
        "min_epitope_size": "Minimum length over overlapping similarly reactive peptides to be considered an epitope.",
        "data_output_dir": "Name of the directory to be created to place the output files."
    },
    outputs = [("reactivity_plot", Visualization), ("c_pos_histogram", Visualization)],
	output_descriptions = {
        "reactivity_plot": "A scatter plot with the c-version zscore as the x-axis and the s-version zscore as the y-axis,"
            " and an additional for filtering by the number of cysteines in the parent peptide.",
        "c_pos_histogram": "A histogram which shows the number of cysteins at given positions for peptides that are"
            " similary or more reactive with different versions."
    },
    name = "Compare C to S substitution reactivity",
    description = "Compares the reactivity of C-version (unmutated) and S-version (mutated peptides) "
        " and output a .qzv files for the user to visualize information related to versions with higher reactivity"
        " for each peptide."
)