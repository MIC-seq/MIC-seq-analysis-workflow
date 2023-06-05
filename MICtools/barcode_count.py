"""
===================================================
barcodecount - count kraken species identification 
===================================================

This script use the result of kraken calculation and the output file contains the species ident of each bacteria 
"""

import copy
import os, sys
from collections import defaultdict
from MICtools import taxontree
#import taxontree



def output_file(file_header):
	output_filename = file_header + ".barcode_count.txt"
	tmp_filename = file_header + ".compose.txt"
	return [output_filename,tmp_filename]


def IFNOintersect(list1,list2):
	count = 0
	for i in list1:
		if i in list2:
			count = count + 1

	if count == 0 :
		return True
	else:
		return False


def output_process(kraken_output):
	
	dict_sum = {}	#sum of reads in one barcode
	dict_classified = {}	#sum of classified reads in one barcode
	dict_index = defaultdict(list)	#list of taxonomic index of each mapped reads in one barcode

	for line in kraken_output:
		split_str = line.strip().split('\t')
		#Count reads number of every barcode
		seq_header = split_str[1]
		barcode = seq_header.split("_")[1]
		if barcode in dict_sum:
			dict_sum[barcode] = dict_sum[barcode] + 1 
		else:
			dict_sum[barcode] = 1
		#Count classified reads number  and collect taxonomic index of each read
		if line.startswith('C'):
			dict_index[barcode].append(split_str[2])
			if barcode in dict_classified:
				dict_classified[barcode] = dict_classified[barcode] + 1
			else:
				dict_classified[barcode] = 1 

	return dict_sum, dict_classified, dict_index

def barcode_component(filename,dict_index):

	str_index = ""
	with open(filename, 'w') as f1:
		for bc in dict_index.keys():
			f1.write('{}\t'.format(bc))
			for index in dict_index[bc]:
				str_index = str_index + index + ","
			f1.write('{}\n'.format(str_index[:-1]))

	return 0

def barcodecount(kraken_output, kraken_report, file_header):
	
	filename = output_file(file_header)
	f = open(filename[0],'w')

	dict_sum, dict_classified, dict_index = output_process(kraken_output)
	#barcode_component(filename[1],dict_index)

	lvl_dict, gdict, sdict = taxontree.rebuild(kraken_report)

	#The threshold of index reads number which would be count
	if min(dict_classified.values()) >=500:
		threshold = float(5/min(dict_classified.values()))
	else:
		threshold = 0.01

	#sort barcode based on reads sum
	dict_sum = sorted(dict_sum.items(),key = lambda x:x[1], reverse = True)

	for barcode, reads_sum in dict_sum:

		dict_index_count = {} #mapped reads number of each index in one barcode

		#no classified reads in the barcode
		if barcode not in dict_classified:
			f.write('{}\t{}\t0\t0\t0\tUnindentifiable\t0\t0\tUnindentifiable\n'.format(barcode, reads_sum))
			continue
		else:
			classified_sum = dict_classified[barcode]
			#the reads number of index in one barcode all lower than threshold 
			for index in set(dict_index[barcode]):
				index_count = dict_index[barcode].count(index)
				if index_count >= len(dict_index[barcode]) * threshold:
					dict_index_count[index] = index_count
			if dict_index_count == {}:
				f.write('{}\t{}\t{}\t0\t0\tUnindentifiable\t0\t0\tUnindentifiable\n'.format(barcode, reads_sum, classified_sum))
				continue
			#index not in the result of kraken report process
			del_list = []
			for index in dict_index_count.keys():
				if index not in lvl_dict.keys():
					del_list.append(index)
			if del_list:
				for index in del_list:
					del dict_index_count[index]
					
			if dict_index_count == {}:
				f.write('{}\t{}\t{}\t0\t0\tUnindentifiable\t0\t0\tUnindentifiable\n'.format(barcode, reads_sum, classified_sum))
				continue
			#genus and species level annotation
			if IFNOintersect(dict_index_count.keys(), gdict.keys()) and IFNOintersect(dict_index_count.keys(), sdict.keys()):
				f.write('{}\t{}\t{}\t0\t0\tUnindentifiable\t0\t0\tUnindentifiable\n'.format(barcode, reads_sum, classified_sum))
				continue
			elif IFNOintersect(dict_index_count.keys(), sdict.keys()):
				#reads index count based on relationship 
				dict_genus_value = {}
				for id in dict_index_count:
					for genus_id in gdict.keys():
						if id in gdict[genus_id]:
							if genus_id not in dict_genus_value.keys():
								dict_genus_value[genus_id] = dict_index_count[id]
							else:
								dict_genus_value[genus_id] = dict_genus_value[genus_id] + dict_index_count[id]
				genus_value = max(dict_genus_value.values())
				genus_index = max(dict_genus_value, key = lambda x : dict_genus_value[x])
				genus_rate = genus_value/sum(dict_index_count.values())
				genus_name = lvl_dict[genus_index]
				f.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t0\t0\tUnindentifiable\n'.format(barcode, reads_sum, classified_sum, 
									 genus_index, genus_rate, genus_name))
			elif IFNOintersect(dict_index_count.keys(), gdict.keys()):
				dict_species_value = {}
				for id in dict_index_count:
					for species_id in sdict.keys():
						if id in sdict[species_id]:
							if species_id not in dict_species_value.keys():
								dict_species_value[species_id] = dict_index_count[id]
							else:
								dict_species_value[species_id] = dict_species_value[species_id] + dict_index_count[id]	
				species_value = max(dict_species_value.values())
				species_index = max(dict_species_value, key = lambda x : dict_species_value[x])
				species_rate = species_value/sum(dict_index_count.values())
				species_name = lvl_dict[species_index]
				genus_index = sdict[species_index][-2]
				genus_rate = species_rate
				genus_name = lvl_dict[genus_index]
				f.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t{}\t{:.4f}\t{}\n'.format(barcode, reads_sum, classified_sum, 
								  genus_index, genus_rate, genus_name, species_index, species_rate, species_name))																				
			else:
				dict_genus_value = {}
				dict_species_value = {}
				for id in dict_index_count:
					for genus_id in gdict.keys():
						if id in gdict[genus_id]:
							if genus_id not in dict_genus_value.keys():
								dict_genus_value[genus_id] = dict_index_count[id]
							else:
								dict_genus_value[genus_id] = dict_genus_value[genus_id] + dict_index_count[id]				
					for species_id in sdict.keys():
						if id in sdict[species_id]:
							if species_id not in dict_species_value.keys():
								dict_species_value[species_id] = dict_index_count[id]
							else:
								dict_species_value[species_id] = dict_species_value[species_id] + dict_index_count[id]
				genus_value = max(dict_genus_value.values())
				genus_index = max(dict_genus_value, key = lambda x : dict_genus_value[x])
				genus_rate = genus_value/sum(dict_index_count.values())
				genus_name = lvl_dict[genus_index]
				species_value = max(dict_species_value.values())
				species_index = max(dict_species_value, key = lambda x : dict_species_value[x])
				species_rate = species_value/sum(dict_index_count.values())
				species_name = lvl_dict[species_index]
				f.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t{}\t{:.4f}\t{}\n'.format(barcode, reads_sum, classified_sum, 
								  genus_index, genus_rate, genus_name, species_index, species_rate, species_name))

	f.close()
	return 0


if __name__ == '__main__':
	f1 = open("/public2/labmember/wyc_rawdata/F1/analysis_result/taxon_annotation/F1.output", 'r')
	f2 = open("/public2/labmember/wyc_rawdata/F1/analysis_result/taxon_annotation/F1.bracken.report",'r')
	kraken_output = f1.readlines()
	kraken_report = f2.readlines()
	file_header = "F1_test2"
	barcodecount(kraken_output, kraken_report, file_header)

