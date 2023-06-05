
import re

def genus_anno(barcode_count,prefix):
    bc_list = [line.split()[0].strip() for line in barcode_count]
    genus_info = [line.split()[5].strip() for line in barcode_count]
    
    new_genus_info = []
    for genus in genus_info:
        if genus == "Unindentifiable":
            new_genus_info.append(genus)
        else:
            new_genus_info.append(re.findall(r"g_(.*)",genus)[0])

    outputfile = prefix + '.genus_info.txt'
    with open(outputfile,'w') as f:
        for i in range(0,len(barcode_count)):
            f.write('{}\t{}\n'.format(bc_list[i], new_genus_info[i]))

    return 0


def species_anno(barcode_count,prefix):
    bc_list = [line.split()[0].strip() for line in barcode_count]
    species_info = [line.split("\t")[8].strip() for line in barcode_count]

    new_species_info = []
    for species in species_info:
        if len(species.split(" ")) > 1:
            seg_list = species.split(" ")
            first_part = re.findall(r"s_(.*)",seg_list[0])[0]
            if first_part[0].isalpha():
                first_info = first_part[0] + "."
            else:
                first_info = first_part
            second_info = ""
            for seg in seg_list[1:]:
                second_info = second_info + seg
            species = first_info.strip()+" "+second_info
            new_species_info.append(species)
        elif species == "Unindentifiable":
            new_species_info.append(species)
        else:
            species_tmp = re.findall(r"s_(.*)",species)
            new_species_info.append(species_tmp)

        outputfile = prefix + '.species_info.txt'
    with open(outputfile,'w') as f:
        for i in range(0,len(barcode_count)):
            f.write('{}\t{}\n'.format(bc_list[i], new_species_info[i]))

    return 0


