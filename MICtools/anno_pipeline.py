from shlex import quote
import subprocess as sp
from MICtools import output_filter as of

def process(sample,kraken_ref, prefix):

    sample = quote(sample)
    ref = quote(kraken_ref)
    prefix = quote(prefix)

    command1 = 'kraken2 --db {} --threads 24 --report {}.report --output {}.output {}'.format(ref,prefix,prefix,sample)
    command2 = 'bracken -d {} -t 32 -i {}.report -o {}.bracken -w {}.bracken.report -r 100 -l S'.format(ref,prefix,prefix,prefix)
    command3 = 'MICtools anno -m barcode_count -ko {}.filter.output -kr {}.bracken.report -p {}'.format(prefix,prefix,prefix)
    command4 = 'MICtools anno -m sample_compo -bc {}.barcode_count.txt -l genus -p {}'.format(prefix,prefix)
    command5 = 'MICtools anno -m sample_compo -bc {}.barcode_count.txt -l species -p {}'.format(prefix,prefix)

    c1 = sp.run(command1, shell=True)
    of.filter('{}.output'.format(prefix),prefix)
    c2 = sp.run(command2, shell=True)
    c3 = sp.run(command3, shell=True)
    c4 = sp.run(command4, shell=True)
    c5 = sp.run(command5, shell=True)

    return 0
