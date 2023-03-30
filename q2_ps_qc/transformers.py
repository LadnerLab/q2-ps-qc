## this code has been moved to q2-pepsirf ##

# #!/usr/bin/env python
# from q2_ps_qc.plugin_setup import plugin
# from q2_ps_qc.format_types import PepsirfContingencyTSVFormat

# from q2_types.feature_table import BIOMV210Format

# import pandas as pd
# import biom


# # Then BIOMV210Format -> BIOMV210DirFmt
# @plugin.register_transformer
# def _0(ff: PepsirfContingencyTSVFormat) -> BIOMV210Format:
#     result = BIOMV210Format()

#     dataframe = pd.read_csv(str(ff), sep='\t', index_col=0)
#     table = biom.Table(dataframe.values, observation_ids=dataframe.index,
#                        sample_ids=dataframe.columns)

#     with result.open() as fh:
#         table.to_hdf5(fh, generated_by="q2-ps-qc for pepsirf")

#     return result



# @plugin.register_transformer
# def _1(ff: BIOMV210Format) -> PepsirfContingencyTSVFormat:
#     result = PepsirfContingencyTSVFormat()

#     with ff.open() as fh:
#         table = biom.Table.from_hdf5(fh)
#     df = table.to_dataframe(dense=True)
#     df.index.name = "Sequence name"
#     df.to_csv(str(result), sep='\t')

#     return result
