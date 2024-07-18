#!/usr/bin/env python
from q2_ps_qc.actions.generate_corr_matrix import generate_corr_matrix
from q2_ps_qc.actions.filter_counts_matrix import filter_counts_matrix
from q2_ps_qc.actions.filter_counts_matrix_tsv import filter_counts_matrix_tsv

__all__ = ['generate_corr_matrix', 'filter_counts_matrix', 'filter_counts_matrix_tsv']

from . import _version
__version__ = _version.get_versions()['version']
