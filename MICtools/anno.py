'''
================================================================
MICtools anno - single microbe bacterial taxonomic identification
================================================================

Main function: 
    This module takes the fastq file with extracted barcode and UMI information 
    adding to the read name as input and is able to process following analysis:
    1) Counting taxonomic index of reads in each microbe and finally determining 
    the genus and species of the single microbe. It also generates the confidence 
    of each prediction level
    2) Examination of taxonomic composition based on each reads for some unindentifiable 
    microbes
    3) Collecting the genus or species level information of each microbe in the sample, 
    and displaying the main taxonomic composition of the sample in a pieplot
    4) Simplified command to process single-microbe taxonomic annotation 
'''

description = '''
Method:

    Counting reads index to determine the taxon of each microbe:

        MICtools anno --module barcode_count --kraken-output [kraken_output] --kraken-report [kraken_report] --prefix [prefix]

    Collecting taxonomic information in a specific microbe and present in pieplot format:

        MICtools anno --module bc_compo --barcode [barcode_info] --kraken-output [kraken_output] --kraken-report [kraken_report]

    Displaying sample taxon composition of cells:

        MICtools anno --module sample_compo --barcode-count [barcode_count] --level [genus/species] --prefix [prefix]

    MIC-seq anno analysis directly:

        MICtools anno --module pipeline --sample [Pair2_filename] --kraken-ref [reference_path] --prefix [prefix]

    Note: The "barcode_count" model is recommended to process kraken output and report file after treated and filtered. It is helpful to 
    filter short reads in kraken output file (<60~80bp).  For kraken report file, bracken analysis pipeline is recommended to process 
    and "kreport2mpa.py" can help know the relationship of each taxonomic level. To konw the detailed step, please read the analysis 
    pipeline
 
'''

import re
import os
import sys
import argparse
from MICtools import barcode_count as bc
from MICtools import sample_composition as sc
from MICtools import bc_composition as bcp
from MICtools import anno_pipeline as ap

def main(argv):
    """script main

    parser command line options in sys.argv, unless *argv* is given
    """

    if argv is None:
        argv = sys.argv

    if argv[0] != "anno" and argv[0] != "anno.py":
        print("Please input appropriate tool name!\nspindent closed")
        
        return 

    parser = argparse.ArgumentParser(
                                    usage = globals()["__doc__"],
                                    description = description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    group = parser.add_argument_group("anno-specific options")

    group.add_argument('-m', '--module', choices = ['barcode_count','sample_compo','bc_compo','pipeline'], 
                       required = True, help=
                        "Choose the model to analysis single microbe RNA-seq data and\n"
                        "detailed command format is in method part")

    group.add_argument('-b', '--barcode', type = str, default = None, help=
                        "Text file composed of barcode information in each line")
    
    group.add_argument('-bc', '--barcode-count', type = str, default = None, help=
                        "Please input the \"\{\}.barcode_count.txt\" filepath")
                                                                                
    group.add_argument('-ko', '--kraken-output', type = str, default = None, help=
                        "Please input the output of 'kraken2 --output' filepath")

    group.add_argument('-kr', '--kraken-report', type = str, default = None, help=
                        "Please input the output of 'kraken2 --report' filepath")   
    
    group.add_argument('-l', '--level' ,choices = ['genus', 'species'],  help =
                        "Choose the specific taxonomic level to display the sample composition")
    
    group.add_argument('-p', '--prefix', type = str, help=
                    "output filename prefix of anno module")

    group.add_argument('-r', '--kraken-reference' ,type = str, help =
                        "Please input the filepath of kraken2 reference datasets folder")

    group.add_argument('-s', '--sample' ,type = str, default = None, help =
                        "Please input the filepath of R2 end smRNA-seq fastq file.")


    args = parser.parse_args()

    if args.module == "barcode_count":
        with open(args.kraken_output,'r') as f1:
            kraken_output = f1.readlines()
        with open(args.kraken_report,'r') as f2:
            kraken_report = f2.readlines()
        file_header = args.prefix
        bc.barcodecount(kraken_output, kraken_report, file_header)        
    elif args.module == "sample_compo":
        with open(args.barcode_count,'r') as f1:
            barcode_count = f1.readlines()
        level = args.level
        prefix = args.prefix
        sc.level(level, barcode_count, prefix)
    elif args.module == "bc_compo":
        with open(args.barcode,'r') as f1:
            barcode_input = f1.readlines()
        with open(args.kraken_output,'r') as f2:
            kraken_output =f2.readlines()
        with open(args.kraken_report,'r') as f3:
            kraken_report = f3.readlines()
        bcp.listprocess(barcode_input,kraken_output,kraken_report)
    else:
        sample = args.sample
        kraken_ref = args.kraken_reference
        prefix = args.prefix
        ap.process(sample,kraken_ref,prefix)


if __name__ == '__main__':
    sys.exit(main(sys.argv))