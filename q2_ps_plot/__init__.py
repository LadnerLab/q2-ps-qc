#!/usr/bin/env python
from q2_ps_plot.actions.zenrich import zenrich
from q2_ps_plot.actions.zenrich_tsv import zenrich_tsv
from q2_ps_plot.actions.boxplot import readCountsBoxplot
from q2_ps_plot.actions.scatter import repScatters
from q2_ps_plot.actions.boxplot import enrichmentRCBoxplot
from q2_ps_plot.actions.heatmap import proteinHeatmap

__all__ = ['zenrich', 'zenrich_tsv', 'readCountsBoxplot', 'repScatters', 'enrichmentRCBoxplot', 'proteinHeatmap']

from . import _version
__version__ = _version.get_versions()['version']
