#!/usr/bin/bash


###microbiome part###
#####################

sample=A197

###uhgg mapping
STAR --runThreadN 32 --runMode alignReads --genomeDir ${ref_dir} --readFilesIn ${fq} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}_ --outReadsUnmapped Fastx --winAnchorMultimapNmax 50 --outFilterMultimapNmax 10 --alignIntronMin 0 --alignIntronMax 0  --outBAMsortingThreadN 39
featureCounts -T 32 -t transcript --extraAttributes gene_id,gene_name -a ${ref_dir}/species_genome.gtf -o ${sample}_assigned.txt -R BAM ${sample}_Aligned.sortedByCoord.out.bam
samtools sort -@ 20 ${sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam -o ${sample}_assigned_sorted.bam
samtools index ${sample}_assigned_sorted.bam
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I ${sample}_assigned_sorted.bam -S ${sample}_cell_matrix.txt