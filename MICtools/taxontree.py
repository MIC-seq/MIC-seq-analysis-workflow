description = '''
==================================
taxonnomic relationship tree build
==================================
'''

import os, sys
from collections import defaultdict

def process_kraken_report(line):
    split_str = line.strip().split('\t')
    try:
        int(split_str[4])
    except ValueError:
        return []
    #get rate
    rate = float(split_str[0])
    #Get taxon level
    level_type = split_str[3]
    level_id = split_str[4]
    #Get name 
    name = split_str[-1].strip()
    #Determine level based on number of spaces
    lvl_dict = {'R1' : 1,'D': 2, 'K': 2, 'P': 3, 'C': 4, 'O': 5, 'F': 6, 'G': 7, 'S': 8}
    if level_type in lvl_dict.keys():
        level_num = lvl_dict[level_type]
    else:
        level_num = 0
    return [rate, name, level_num, level_type, level_id]


def rebuild(report):

    prev_lvl_num = 0
    lvl_path = []   ###from kingdom to species relationship
    lvl_dict = {}   ###dict for taxon in report, eg.: g_Prevotella:472
    gdict = defaultdict(list)   ###all taxon relationship with a genus 
    sdict = defaultdict(list)   ###all taxon relationship with a species 
    #line processing
    for line in report:
        report_vals = process_kraken_report(line)
        #header line
        if len(report_vals) < 5: 
            continue
        #filter taxon level
        [rate, name, level_num, level_type, level_id] = report_vals
        if rate <= 0.01:
            continue
        if level_num == 0 :
            continue
        if level_type == "K":
            level_type = "D"
        level_str = level_type.lower() + "_" + name
        #build dict and relationship tree
        lvl_dict[level_id] = level_str
        if level_num > prev_lvl_num:
            prev_lvl_num = level_num
            lvl_path.append(level_id)
            if level_num == 7:
                gdict[level_id] = lvl_path.copy()
            if level_num == 8:
                gdict[lvl_path[-2]].append(level_id)
                sdict[level_id] = lvl_path.copy()
        elif level_num == prev_lvl_num:
            lvl_path.pop()
            lvl_path.append(level_id)
            if level_num == 7:
                gdict[level_id] = lvl_path.copy()
            if level_num == 8:
                gdict[lvl_path[-2]].append(level_id)
                sdict[level_id] = lvl_path.copy()
        #level_num < prev_lvl_num
        else:  
            lvl_gap = prev_lvl_num - level_num
            lvl_path = lvl_path[:len(lvl_path)-(lvl_gap + 1)]
            lvl_path.append(level_id)
            if level_num == 7:
                gdict[level_id] = lvl_path.copy()
            prev_lvl_num = level_num
            
    return lvl_dict, gdict, sdict


    """
def rebuild_v1(kraken_report)
    ###kraken map transform foramt
    dict = {}   ###dict for taxon in report, eg.: g_Prevotella:472
    taxon_tree = [] ###taxonnomic name information from mpa file
    id_tree = []    ###taxonnomic id based on taxon_tree info and kraken report

    for line in report:
        taxon_level = line.strip().split("\t")[3]
        if taxon_level not in ["P","C","O","F","G","S"]:
            continue

        taxon_id = line.strip().split("\t")[4]
        taxon_name = line.strip().split("\t")[5].strip()
        taxon_level = taxon_level.lower()
        taxon_info = taxon_level + "_" + taxon_name

        dict[taxon_info] = taxon_id


    for detail in mpa:
        taxon_detail = detail.strip().split("\t")[0]
        id_detail = ""

        tmp_list = taxon_detail.split("|")
        for i in range(0,len(tmp_list)):
            if i != len(tmp_list) - 1:
                id_detail = id_detail + dict[tmp_list[i]] + "|"
            else:
                id_detail = id_detail + dict[tmp_list[i]]

        taxon_tree.append(taxon_detail)
        id_tree.append(id_detail)

    return taxon_tree, id_tree
    """


