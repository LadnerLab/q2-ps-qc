#!/usr/bin/env python
from q2_ps_qc.actions.generate_corr_matrix import generate_corr_matrix
from q2_ps_qc.actions.compareCS import compareCS

__all__ = ['generate_corr_matrix', 'compareCS']

from . import _version
__version__ = _version.get_versions()['version']
