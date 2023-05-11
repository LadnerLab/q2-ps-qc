""""#!/usr/bin/env python
import qiime2.plugin.model as model
from qiime2.plugin import SemanticType

class CorrelationMatrixTSVFormat(model.TextFileFormat):
    def _validate_(self, level='min'):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith('Sequence name\t'):
                    raise model.ValidationError(
                        'TSV does not start with "Sequence name"')


CorrelationMatrixDirectoryFormat = model.SingleFileDirectoryFormat(
    'CorrelationMatrixDirectoryFormat',
    'correlation-matrix.tsv',
    CorrelationMatrixTSVFormat)

from q2_ps_qc.plugin_setup import plugin

plugin.register_formats(CorrelationMatrixDirectoryFormat)

CorrelationMatrix = SemanticType('CorrelationMatrix')
plugin.register_semantic_types(CorrelationMatrix)"""

