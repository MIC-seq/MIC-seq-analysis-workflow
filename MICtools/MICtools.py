'''
MICtools - a single-microbe RNA-seq analysis pipeline
====================================================================

:Release: MIC-seq Group
:Verion:  1.0.0
:Date:    2023.05.30
:Tags:    MIC-seq analysis pipeline for single-microbe RNA-seq datasets

There are mainly 3 tools:

  - anno	taxonomic annotation for each microbe in genus or species level 
  - bac 	pipeline for microbiome transcriptional matrix
  - phage	pipeline for transcrptional host-phage reltionship matrix

To get help on a specific tool, type:

	MICtools <tools> --help

To use a specific tools, type:

	MICtools <tool> [tools options] [tool argument]
'''

from __future__ import absolute_import
import os
import sys
import importlib
from MICtools import __version__

def main(argv = None):

  argv = sys.argv

  path = os.path.abspath(os.path.dirname(__file__))

  if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
    print(globals()["__doc__"])

    return

  elif len(argv) == 1 or argv[1] == "--version" or argv[1] == "-v":
    print("MICtools version: %s" % __version__)

    return

  elif argv[2] in ["--help", "-h", "--help-extended"]:
    print("MICtools: Version %s" % __version__)

  command = argv[1]

  module = importlib.import_module("MICtools." + command, "MICtools")
  ##remove 'MICtools' from sys.argv
  del sys.argv[0]
  module.main(sys.argv)

if __name__ == '__main__':
  sys.exit(main())
