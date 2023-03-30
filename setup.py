#!/usr/bin/env python
from setuptools import setup, find_packages
import versioneer


setup(
    name="q2-ps-qc",
    version=versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={},
    author="Johnathan Ray",
    author_email="jdr479@nau.edu",
    description="Visualizations for PepSIRF outputs.",
    license='Apache-2.0',
    url="https://github.com/LadnerLab/q2-ps-qc",
    entry_points={
        'qiime2.plugins': ['q2-ps-qc=q2_ps_qc.plugin_setup:plugin']
    },
    zip_safe=False,
)
