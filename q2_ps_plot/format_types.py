## this code has been moved to q2-pepsirf ##

#!/usr/bin/env python
# import qiime2.plugin.model as model
# from qiime2.plugin import SemanticType

# from q2_types.feature_table import FeatureTable


# Normed = SemanticType('Normed', variant_of=FeatureTable.field['content'])
# Zscore = SemanticType('Zscore', variant_of=FeatureTable.field['content'])


# class PepsirfContingencyTSVFormat(model.TextFileFormat):
#     def _validate_(self, level='min'):
#         with self.open() as fh:
#             for _, line in zip(range(1), fh):
#                 if not line.startswith('Sequence name\t'):
#                     raise model.ValidationError(
#                         'TSV does not start with "Sequence name"')


# PepsirfContingencyTSVDirFmt = model.SingleFileDirectoryFormat(
#     'PepsirfContingencyTSVDirFmt',
#     'pepsirf-table.tsv',
#     PepsirfContingencyTSVFormat)
