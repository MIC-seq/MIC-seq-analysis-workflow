#!/usr/bin/bash

samtools view nonrna_phage_test_assigned_sorted.bam|cut -f1 > mapping_reads_bam_header.txt

bam_file=mapping_reads_bam_header.txt
whitelist=/public2/labmember/wyc_rawdata/A197/supplement_data/data_analysis/whitelist.txt
bc_info_folder=/public2/labmember/wyc_rawdata/A197/supplement_data/phage_test/bc_reads_mapping_rate

tmp_fifo="$$.fifo"
mkfifo $tmp_fifo
exec 6<>$tmp_fifo
rm $tmp_fifo

thread_num=30

for ((i=0;i<${thread_num};i++))
do
        echo
done >&6

for line in `cut -f1 $whitelist`;
do
	read -u6
	{
		mapping_count=`grep $line ${bam_file}|sort|uniq|wc -l`
		reads_count=`grep $line $whitelist|cut -f2`
		rate=$(printf "%.3f" `echo "scale=3;${mapping_count}/${reads_count}"|bc`)
		echo -e "${line}\t${mapping_count}\t${reads_count}\t${rate}" >>result/A197_bc_phage_rate.txt
		echo >&6
	}&
done

wait

exec 6>&-

echo "mission completed"
