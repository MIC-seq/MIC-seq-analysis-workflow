from __future__ import absolute_import

import sys
import setuptools
from MICtools import __version__

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print ("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "MICtools requires setuptools 1.1 higher")

###############################################################
###############################################################
# Define dependencies 
# Perform a MICtools Installation

major, minor1, minor2, s, tmp = sys.version_info

if (major == 3 and minor1 < 6) or major < 3:
    raise SystemExit("""umi_tools requires Python 3.6 or later.""")

MICtools_packages = ["MICtools"]
MICtools_package_dirs = {'MICtools': 'MICtools'}   


install_requires = [
    "future",
    "matplotlib",
    "numpy>=1.7",
    "pandas>=1.4",
    "matplotlib>=3.5",
    "palettable>=3.3.3"]

setuptools.setup(
    name="MICtools",
    version="1.0.0",
    author="Gibbs_Qian",
    author_email="12016020@zju.edu.cn",
    description="MIC-seq analysis pipeline",
    long_description="Single-micobre RNA-seq for microbiome sample",
    long_description_content_type="text",
    # package contents
    #packages=setuptools.find_packages(),
    packages=MICtools_packages,
    package_dir=MICtools_package_dirs,
    include_package_data=True,
    # dependencies
    install_requires = install_requires,
    entry_points={
        'console_scripts': ['MICtools = MICtools.MICtools:main']
    }
)
