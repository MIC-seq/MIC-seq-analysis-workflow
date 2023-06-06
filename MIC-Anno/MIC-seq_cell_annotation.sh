#!/usr/bin/bash


sample={sample_name}
ref_dir={uhgg_datasets/gpd_datasets file folder path}
sortmerna_idxpath={sortmerna index file folder path}

###Data filters
#10000 microbes for example
zcat ${sample}_1.fq.gz |grep -e "^[ATCG]" |cut -c-20|sort|uniq -c|sort -k1 -nr|awk -v OFS="\t" '{print $2,$1}'|head -10000 > whitelist.txt

umi_tools extract --extract-method=regex --bc-pattern="(?P<cell_1>.{20})(?P<umi_1>.{8}).*" --stdin ${sample}_1.fq --stdout ${sample}_R1_extracted.fq.gz --read2-in ${sample}_2.fq --read2-out=${sample}_R2_extracted.fq.gz --whitelist whitelist.txt --filter-cell-barcode --error-correct-cell


###taxonomic annotation###
##########################
kraken2 --db /public2/labmember/qianqh/sc_microbio/kraken2_db_uhgg_v2.0.1 --threads 24 --report ${sample}.report --output ${sample}.output ${sample}_R2_extracted.fq.gz
awk -F "\t" '($4>=60){print $0}' ${sample}.output > ${sample}.filter.output
bracken -d /public2/labmember/qianqh/sc_microbio/kraken2_db_uhgg_v2.0.1 -t 32 -i ${sample}.report -o ${sample}.bracken -w ${sample}.bracken.report -r 100 -l S
python /public2/labmember/qianqh/sc_microbio/scripts/kreport2mpa.py -r ${sample}.bracken.report -o ${sample}.bracken.report.mpa
perl /public2/labmember/qianqh/sc_microbio/scripts/go_perl_newx.pl ${sample}.bracken.report ${sample}.bracken.report.mpa ${sample}.bracken.report.mpa.newx
python barcode_count.py -ko ${sample}.filter.output -kr ${sample}.bracken.report.mpa.newx -p ${sample}
##genus level annotation
awk -F "\t" '{print $1}' ${sample}.barcode_count.txt > bc_tmp.txt
awk -F "|" '{print $5}' ${sample}.barcode_count.txt| sed -e 's/g_\(.*\)/\1/' > ge_tmp.txt
sed -e 's/^\s*$/Unknown_genus/' ge_tmp.txt > ge_tmp2.txt
paste bc_tmp.txt ge_tmp2.txt > ${sample}_bc_genus_info.txt
rm bc_tmp.txt ge_tmp.txt ge_tmp2.txt
##species level annotation
awk -F "\t" '{print $1}' ${sample}.barcode_count.txt > bc_tmp.txt
awk -F "|" '{print $6}' ${sample}.barcode_count.txt| sed -e 's/s_\(.*\)/\1/' > sp_tmp.txt
cut -b1 sp_tmp.txt > sp_tmp1.txt
awk '{print $2}' sp_tmp.txt > sp_tmp2.txt
paste sp_tmp1.txt sp_tmp2.txt >sp_tmp3.txt
sed -e 's/\t/./' -e 's/_[A-Za-z]*//' -e 's/^.$/Unknown_species/' sp_tmp3.txt > sp_tmp4.txt
paste bc_tmp.txt sp_tmp4.txt > ${sample}_bc_species_info.txt
rm bc_tmp.txt sp_tmp.txt sp_tmp1.txt sp_tmp2.txt sp_tmp3.txt sp_tmp4.txt


