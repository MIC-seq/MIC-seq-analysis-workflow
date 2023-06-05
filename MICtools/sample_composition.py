

###draw pieplot for sample taxon composition###
#############################################

import argparse
import matplotlib.pyplot as plt
from MICtools import taxon_anno_count as tac
import numpy as np
import pandas as pd
import palettable 
import sys

def pieplot(taxon, num, prefix):

	my_dpi=120
	plt.figure(figsize=(1200/my_dpi,1200/my_dpi),dpi=my_dpi)
	if len(num) % 9 == 1:
		plt.pie(x = num,
	 		labels = taxon,
			colors = palettable.cartocolors.qualitative.Bold_9.mpl_colors[1:],
			autopct='%.1f%%',
			pctdistance=0.8,
			wedgeprops={'alpha':0.6}
		)
	else:
		plt.pie(x = num,
	 		labels = taxon,
			colors = palettable.cartocolors.qualitative.Bold_9.mpl_colors,
			autopct='%.1f%%',
			pctdistance=0.8,
			wedgeprops={'alpha':0.6}
		)

	plt.savefig("{}.pdf".format(prefix))

	return 0


def sample_component(taxon_list):

	sample_info = pd.value_counts(taxon_list)

	thread = int(np.sum(sample_info) * 0.01)
	rare_taxon = sample_info.index[sample_info <= thread]
	rare_taxon_num = np.sum(sample_info[sample_info<=thread])
	sample_info.drop(rare_taxon.tolist(), inplace = True)
	sample_info['others'] = rare_taxon_num
	
	taxon = sample_info.index.tolist()
	num = sample_info.tolist()

	return taxon, num
	
def level(level, barcode_count, prefix):
	
	if level == 'genus':
		tac.genus_anno(barcode_count, prefix)
		compofile = prefix + '.genus_info.txt'
		with open(compofile) as f:
			sample_bc_info = f.readlines()
		taxon_list = [line.split()[1].strip() for line in sample_bc_info]
		prefix = prefix + "_genus"
	else:
		tac.species_anno(barcode_count, prefix)
		compofile = prefix + '.species_info.txt'
		with open(compofile) as f:
			sample_bc_info = f.readlines()
		taxon_list = [line.split("\t")[1].strip() for line in sample_bc_info]
		prefix = prefix + "_species"

	taxon, num = sample_component(taxon_list)
	pieplot(taxon, num, prefix)	

	return 0 

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('-i','--input',type = str, default = None, help = "text file of taxon information")
	parser.add_argument('-p','--prefix',type = str, default = None,help = "prefix of pieplot")

	args = parser.parse_args()

	if len(sys.argv) ==1:
		print("\nplease type -h or --help to check the useage of the scripts")
	else:
		with open(args.input, 'r') as f:
			lines = f.readlines()

		prefix = args.prefix

		taxon,num = sample_component(lines)
		pieplot(taxon, num, prefix)


if __name__ == '__main__':
	sys.exit(main())
