#!/usr/bin/env python
import importlib
import q2_ps_qc

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

import q2_ps_qc.actions as actions


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
        "reps_source": Str,
        "log_normalization": Bool,
        "correlation_threshold": Float
    },
    parameter_descriptions = {
		"data": "Name of input file.",
        "reps_source": "A file containing sample names which will be"
            " considered as replicates during the process.",
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

