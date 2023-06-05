#!/public/home/qianqh/softwares/anaconda3/bin/python

"""
=======================================================
cellconpoment - showing reads conpoment of one microbe 
=======================================================
"""
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt 
import palettable
import pandas as pd 
import sys


def pieplot(taxon,num,prefix):
	my_dpi = 120
	plt.figure(figsize=(1500/my_dpi,1200/my_dpi),dpi=my_dpi)
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

	plt.savefig("microbe_{}.pdf".format(prefix))

	return 0

def microbe_component(bc,bc_id,id_name):

	bc = bc.strip()
	id_num = bc_id.id[bc_id.bc == bc].value_counts()
	id_list = id_num.index.tolist()
	for id in id_list:
		if id not in id_name.keys():
			id_num.drop(id, inplace = True)
	id_list_modified = id_num.index.tolist()
	taxon = [id_name[id] for id in id_list_modified]
	num = id_num.tolist()

	pieplot(taxon, num, bc)


def listprocess(barcode_input,kraken_output,kraken_report):

	bc_list = [read.strip().split('_')[1] for read in kraken_output if read.startswith('C')]
	id_list = [read.strip().split('\t')[2] for read in kraken_output if read.startswith('C')]
	bc_id = {"bc":bc_list, "id":id_list}
	bc_id = pd.DataFrame(bc_id)

	id_name = {}
	for line in kraken_report:
		id= line.split('\t')[4]
		tax = line.split('\t')[3]
		name = line.split('\t')[5].strip()
		name = tax.lower() + "_" + name
		id_name[id] = name

	for bc in barcode_input:
		microbe_component(bc,bc_id,id_name)

	return 0 

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('-b','--bc',type = str, default = None,help = "text file of barcode list")
	parser.add_argument('-i','--input',type = str, default = None, help = "kraken output")
	parser.add_argument('-r','--report',type = str, default = None,help = "kraken report")

	args = parser.parse_args()

	if len(sys.argv) ==1:
		print("\nplease type -h or --help to check the useage of the scripts")
	else:
		with open(args.bc, 'r') as f:
			barcode_input = f.readlines()		
		with open(args.input, 'r') as f:
			kraken_output = f.readlines()
		with open(args.report, 'r') as f:
			kraken_report = f.readlines()

	listprocess(barcode_input,kraken_output,kraken_report)


if __name__ == '__main__':
	sys.exit(main())
