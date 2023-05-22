#!/usr/bin/bash

###phage part ###
#################

sample=A197
##filter rRNA reads
  cp -r {sortmerna_idxpath}/idx ./
  sortmerna --threads 16 \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-bac-16s-id90.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-bac-23s-id98.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-euk-18s-id95.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-euk-28s-id98.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-arc-16s-id95.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/silva-arc-23s-id98.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/rfam-5s-database-id98.fasta \
  --ref /public2/labmember/qianqh/sc_microbio/SortMeRNA_ref_db/silva_rRNA_database/rfam-5.8s-database-id98.fasta \
  --workdir ./ --reads ${sample}_R2_extracted.fq.gz --aligned "split_rrna/${sample}_rRNA" --other "split_rrna/${sample}_non_rRNA" --fastx
  rm -r kvdb/ readb/ idx/

STAR --runThreadN 32 --runMode alignReads --genomeDir ${ref_dir} --readFilesIn split_rrna/${sample}_non_rRNA.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}_ --outReadsUnmapped Fastx --winAnchorMultimapNmax 50 --outFilterMultimapNmax 10 --alignIntronMin 0 --alignIntronMax 0  --outBAMsortingThreadN 39

featureCounts -T 32 -t exon --extraAttributes gene_id,gene_biotype -a ${ref_dir}/GPD_sequences.gtf -o ${sample}_assigned.txt -R BAM ${sample}_Aligned.sortedByCoord.out.bam
samtools sort -@ 20 ${sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam -o ${sample}_assigned_sorted.bam
samtools index ${sample}_assigned_sorted.bam
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I ${sample}_assigned_sorted.bam -S ${sample}_cell_matrix.txt