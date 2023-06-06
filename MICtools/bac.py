'''
==========================================================
MICtools bac - single microbe bacterial matrix construction
==========================================================

Main function: 
    This module takes fastq file and the reference fasta and gtf file as input. The fastq 
    file is recommended to contain extracted barcode and UMI information adding to the 
    read name and reference datasets recommends the corresponding microbiome genome 
    datasets (like uhgg for human gut microbiome analysis). And then this pipeline can 
    generates transcriptional activity matrix of individual bacterium in a sample. 

'''
description = '''
Method:

    Construction of microbial transcriptional activity matrix:

        MICtools bac --sample [sample_file] --reference [ref_folder] --feature_type [gene/exon/transcript...] --prefix [prefix]

    Note: Sample file is the R2 end file of fastq which contain extracted barcode and UMI information adding to the 
    read name. And it is recommended to contain only reference genome fasta and gtf file in ref_folder. The feature type 
    is the major gene information while presentation format differs in different gtf files in column #3.Some files are genes 
    and some are exon or transcript and so on.
 
'''

import re
import os, sys
import argparse
from shlex import quote
import subprocess as sp

def main(argv):

    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
                                    usage = globals()["__doc__"],
                                    description = description,
                                    formatter_class=argparse.RawTextHelpFormatter)
    
    group = parser.add_argument_group("bac-specific options")

    group.add_argument('-s', '--sample', type = str, default = None, help=
                        "R2 end fastq file path")
    
    group.add_argument('-r', '--reference', type = str, default = None, help=
                        "Folder contains microbial reference genome fasta and its gtf file")
    
    group.add_argument('-f', '--feature_type', type = str, default = None, help=
                        "feture type in gtf file column#3")

    group.add_argument('-p', '--prefix', type = str, default = None, help=
                        "output filename prefix if bac module")
    
    args = parser.parse_args()

    sample = quote(args.sample)
    ref = quote(args.reference)
    feature = quote(args.feature_type)
    prefix = quote(args.prefix)


    ##dependent software path
    star_path = sp.run('which STAR',shell= True,capture_output=True)
    featurecounts_path = sp.run('which featureCounts',shell= True,capture_output=True)
    samtools_path = sp.run('which samtools',shell= True,capture_output=True)
    umitools_path = sp.run('which umi_tools',shell= True,capture_output=True)

    if star_path.returncode == 0:
        path = str(star_path.stdout,'utf-8').strip()
        print('STAR path: {}'.format(path))
    else:
        print('STAR is not found! Please add the STAR software to \$PATH')
        sys.exit()
    
    if featurecounts_path.returncode == 0:
        path = str(featurecounts_path.stdout,'utf-8').strip()
        print('featureCounts path: {}'.format(path))
    else:
        print('featureCounts is not found! Please add the featurecounts software to \$PATH')
        sys.exit()

    if samtools_path.returncode == 0:
        path = str(samtools_path.stdout,'utf-8').strip()
        print('samtools path: {}'.format(path))
    else:
        print('samtools is not found! Please add the samtools software to \$PATH')
        sys.exit()
        
    if umitools_path.returncode == 0:
        path = str(umitools_path.stdout,'utf-8').strip()
        print('umi_tools path: {}'.format(path))
    else:
        print('umi_tools is not found! Please add the umi_tools software to \$PATH')
        sys.exit()


    gtf_file = sp.run('ls -l {}|grep .gtf'.format(ref),shell= True,capture_output=True)
    if gtf_file.returncode != 0:
        print('gtf file is not found in reference folder!')
        sys.exit()
    else:
        gtf_file = str(gtf_file.stdout.strip(),'utf-8')
        gtf_file = re.findall(r'(\S*\.gtf)',gtf_file)[0]
        gtf_file = ref + '/' + gtf_file

    fna_file = sp.run('ls -l {}|grep .fna'.format(ref),shell= True,capture_output=True)
    fa_file = sp.run('ls -l {}|grep .fa'.format(ref),shell= True,capture_output=True)
    if fna_file.returncode != 0 and fa_file.returncode != 0:
        print('fna file is not found in reference folder!')
        sys.exit()
    elif fna_file.returncode == 0:
        fna_file = str(fna_file.stdout.strip(),'utf-8')
        fna_file = re.findall(r'(\S*\.fna\S*)',fna_file)[0]
        fna_file = ref + '/' + fna_file
    else:
        fa_file = str(fa_file.stdout.strip(),'utf-8')
        fa_file = re.findall(r'(\S*\.fa\S*)',fa_file)[0]
        fna_file = re.sub('.fa','.fna',fa_file)
        fna_file = ref + '/' + fna_file

    ###whether retain STAR reference datasets 
    file = ref + '/Genome'
    star_ref = sp.run('ls {}'.format(file),shell= True,capture_output=True)
    if star_ref.returncode != 0 :
        ##build STAR ref
        command0 = 'STAR --runMode genomeGenerate --runThreadN 10 --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} '\
                    '--sjdbGTFfeatureExon {} --limitGenomeGenerateRAM 489348748896 --sjdbOverhang 149'.format(ref,fna_file,gtf_file,feature)
        c0 = sp.run(command0, shell=True)

    if re.findall(r'.gz', sample):
        command1 = 'STAR --runThreadN 10 alignReads --genomeDir {} --readFilesIn {} --readFilesCommand zcat '\
                    '--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {}_ --outReadsUnmapped Fastx '\
                    '--winAnchorMultimapNmax 50 --outFilterMultimapNmax 10 --alignIntronMin 0 --alignIntronMax 0 '\
                    '--outBAMsortingThreadN 39'.format(ref, sample, prefix)
    else:
        command1 = 'STAR --runThreadN 10 alignReads --genomeDir {} --readFilesIn {} '\
                    '--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {}_ --outReadsUnmapped Fastx '\
                    '--winAnchorMultimapNmax 50 --outFilterMultimapNmax 10 --alignIntronMin 0 --alignIntronMax 0 '\
                    '--outBAMsortingThreadN 39'.format(ref, sample, prefix)
    
    command2 = 'featureCounts -T 10 -t {} --extraAttributes gene_id,gene_biotype -a {} '\
                '-o {}_assigned.txt -R BAM {}_Aligned.sortedByCoord.out.bam'.format(feature,gtf_file,prefix,prefix)
    
    command3 = 'samtools sort -@ 20 {}_Aligned.sortedByCoord.out.bam.featureCounts.bam '\
                '-o {}_assigned_sorted.bam'.format(prefix, prefix)
    
    command4 = 'samtools index {}_assigned_sorted.bam'.format(prefix)

    command5 = 'umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell '\
    '--wide-format-cell-counts -I {}_assigned_sorted.bam -S {}_microbe_matrix.txt'.format(prefix,prefix)

    c1 = sp.run(command1, shell=True)
    c2 = sp.run(command2, shell=True)
    c3 = sp.run(command3, shell=True)
    c4 = sp.run(command4, shell=True)
    c5 = sp.run(command5, shell=True)
    print("mission complete!")    
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
