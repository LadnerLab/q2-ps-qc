#!/usr/bin/env python
from setuptools import setup, find_packages
import versioneer


setup(
    name="q2-ps-plot",
    version=versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={},
    author="Annabelle Brown",
    author_email="annabelle811@live.com",
    description="Visualizations for PepSIRF outputs.",
    license='Apache-2.0',
    url="https://github.com/LadnerLab/q2-ps-plot",
    entry_points={
        'qiime2.plugins': ['q2-pl-plot=q2_ps_plot.plugin_setup:plugin']
    },
    zip_safe=False,
)
