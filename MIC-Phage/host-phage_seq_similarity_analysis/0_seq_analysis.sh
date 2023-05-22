#!/usr/bin/bash

###top20 variable phages genome test###
#######################################
for i in `cat genus_list.txt`
do

for j in `cat ${i}_markers.txt.top20`
do
grep -A1 ${j} /public3/labmember/qianqh/gut_phage_database/GPD_sequences.fa >> ${i}_markers.txt.top20.fna
done

blastn \
-query ${i}_markers.txt.top20.fna \
-db /public3/labmember/qianqh/uhgg_gut_genome_db_v2.0.1/A197_main_genus_genome/g__${i}/genome \
-evalue 1e-3 \
-outfmt 6 \
-num_threads 6 \
-out ${i}.seq_sim.blastn.out


for j in `cat ${i}_markers.txt`
do
grep ${j} ${i}.seq_sim.blastn.out|head -1 >> ${i}.seq_similarity.txt
done
done



###all variable phages genome(logfc>=0.1) test###
#################################################
for i in `cat genus_list.txt`
do 

for j in `cat ${i}_markers.txt`
do
grep -A1 ${j} /public3/labmember/qianqh/gut_phage_database/GPD_sequences.fa >> ${i}_markers_all.fna
done

blastn \
-query ${i}_markers_all.fna \
-db /public3/labmember/qianqh/uhgg_gut_genome_db_v2.0.1/A197_main_genus_genome/g__${i}/genome \
-evalue 1e-3 \
-outfmt 6 \
-num_threads 6 \
-out ${i}.all.seq_sim.blastn.out


for j in `cat ${i}_markers.txt`
do
grep ${j} ${i}.all.seq_sim.blastn.out|head -1 >> ${i}.all.seq_similarity.txt
done
done
