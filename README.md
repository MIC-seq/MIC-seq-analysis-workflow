# MIC-seq analysis workflow

Analysis methods for scRNA-seq data in mammalian systems have been well developed in recent years. However, these methods are not quite suitable to be applied on single microbe RNA-seq in microbial community. To solve this difficulty, We developed a computational pipeline named MIC-seq analysis workflow for single-microbe RNA sequencing of human gut microbiome. Our pipeline has three major modules: MIC-Anno, MIC-Bac and MIC-Phage. The main function of MIC-Anno is processing taxonomic annotation for each microbe.Based on the results of MIC-Anno, we can use MIC-Bac and MIC-Phage to construct bacteria and phage transcriptional matrix.

This repository contains:

1. MIC-Anno: A pipeline for single microbe taxonomic annotaion
2. MIC-Bac : A pipeline for construction of microbial transcriptional expression matrix and example of further personalized analysis 
3. MIC-Phage: A pipeline for construction of host-phage transcriptional relationship matrix and scripts for further analysis 
4. MICtools: A MIC-seq analysis software contains above three pipeline modules and can get both matrixes with three commands input

And Here is the instructions for usage of MICtools


---
## Table of Contents

MICtools: a single-microbe RNA-seq analysis pipeline

- [Preparation](#preparation)
- [Install](#install)
- [Data preprocess](#data-preprocess)
- [MIC-Anno](#MIC-Anno)
- [MIC-Bac](#MIC-Bac)
- [MIC-Phage](#MIC-Phage)
- [Help](#help)


---
## Preparation

MICtools relies on these following softwares, please ensure that they were installed(recommended version).

- Kraken2 2.1.2 (https://github.com/DerrickWood/kraken2)
- Bracken 2.6.0 (https://github.com/jenniferlu717/Bracken)
- Sortmerna 4.3.4 (https://github.com/sortmerna/sortmerna)
- Sortmerna rRNA reference datasets (https://github.com/sortmerna/sortmerna/blob/master/data/rRNA_databases/silva_ids_acc_tax.tar.gz)
- STAR 2.7.10a (https://github.com/alexdobin/STAR)
- FeatureCounts 2.0.3 (https://sourceforge.net/projects/subread)
- Samtools 1.15 (https://github.com/samtools/samtools)
- Umi_tools 1.1.2 (https://github.com/CGATOxford/UMI-tools)


---
## Install

It is recommended to work directly from the git repository:

```sh
$ git clone  git@github.com:MIC-seq/MIC-seq-analysis-workflow.git

$ cd MIC-seq-analysis-workflow
$ python setup.py sdist
$ pip install dist/MICtools-1.0.0.tar.gz
```


---
## Data preprocoess

Given that smRNA-seq data is quite different from scRNA-seq data for animals or plants, it is not entirely reliable to use current methods to automatically identify microbial counts. So setting microbe numbers based on data quality contral plot is recommanded. Here we assume 5000 microbes as an example and extract barcode and UMI information of the sample dataset:

```sh
$ umi_tools whitelist --stdin {sample_R1_end.fastq} --set-cell-number=5000 --method=umis --plot-prefix={plot_prefix} \
--extract-method=regex --bc-pattern="(?P<cell_1>.{20})(?P<umi_1>.{8}).*" --stdout whitelist.txt --log whitelist.log

$ umi_tools extract --extract-method=regex --bc-pattern="(?P<cell_1>.{20})(?P<umi_1>.{8}).*" --stdin ${sample} \
--stdout ${sample_R1_end.fastq} --read2-in ${sample_R2_end.fastq} --read2-out=${sample}_R2_extracted.fq.gz \
--whitelist whitelist.txt --filter-cell-barcode --error-correct-cell
```


---
## MIC-Anno

This module is mainly used to process taxonomic annotation for each microbe 

### Quick use

```sh
MICtools anno --module pipeline -s /path/${sample}_R2_extracted.fq.gz -r [kraken_ref_filepath] -p [output_prefix]
```

### Major MIC-Anno result

- ${prefix}.barcode_count.txt

	Main single microbe taxonomic annotation result and the meanings of each column:

```txt
#barcode	#reads_count	#classified_reads	#genus_id	#genus_purity	#genus_name	#species_id	#species_purity	#species_name
GCTCTCTCCACCTAGGCAAC	65236	64703	1408	0.9876	g_Megamonas	3246	0.9625	s_Megamonas funiformis
CGCCTGAGCTGACTGGCTGC	61925	61164	1408	1.0000	g_Megamonas	3246	0.9745	s_Megamonas funiformis
```

- {prefix}.genus_info.txt

	Genus info of each microbe

- {prefix}.species_info.txt

	Species info of each microbe

- {prefix}_genus.pdf

	Pieplot of sample genus composition

![](https://github.com/MIC-seq/MIC-seq-analysis-workflow/blob/main/fig/test_genus.png)

- {prefix}_species.pdf

	Pieplot of sample specie composition

![](https://github.com/MIC-seq/MIC-seq-analysis-workflow/blob/main/fig/test_species.png)


---
## MIC-Bac

The main function of this module is to construct single microbe transcriptional matrix 

### Matrix construction
```sh
$ MICtools bac -s /path/${sample}_R2_extracted.fq.gz -r [microbiome_ref_folder_path] -f [gene/exon/transcript...]  -p [output_prefix]
```

### Further analysis

The R script MIC-Bac/microbiome_transcriptional_analysis.R is an example for microbial further analysis

### Note

The feature type of gtf file (the input of premeter '-f') differs in different gtf files column#3. Here in the following example gtf file, we can identify 'transcript' and 'CDS' are main feature types, so we can set '-f transcript' or '-f CDS' to process the analysis.

```txt
MGYG000000001_1 Prodigal:2.6    transcript      188     1672    .       -       .       transcript_id "MGYG000000001_00001"; gene_id "MGYG000000001_00001"; gene_name "clsA_1"
MGYG000000001_1 Prodigal:2.6    CDS     188     1672    .       -       0       transcript_id "MGYG000000001_00001"; gene_name "clsA_1";
MGYG000000001_1 Prodigal:2.6    transcript      1820    2689    .       -       .       transcript_id "MGYG000000001_00002"; gene_id "MGYG000000001_00002"; gene_name "focA_1"
MGYG000000001_1 Prodigal:2.6    CDS     1820    2689    .       -       0       transcript_id "MGYG000000001_00002"; gene_name "focA_1";
```


---
## MIC-Phage

The main function of this module is to construct host-phage transcriptional relationship matrix

### Matrix construction

```sh
$ MICtools phage -s /path/${sample}_R2_extracted.fq.gz -rr [sortmerna_rRNA_reference_datasets_path] -pr [phage_ref_folder_path] -f [gene/exon/transcript...]  -p [output_prefix]
```

### Further analysis

The R script MIC-Phage/analysis_pipeline.R is an example for host-phage trnascriptonal relationship further analysis


---
## Help

To get help on MICtools run

```sh
$ MICtools -h
```

To get help on the options for a specific [module], run

```sh
$ MICtools [COMMAND] --help
```