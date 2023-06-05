
def filter(kraken_output,prefix):
    with open(kraken_output,'r') as f1:
        lines = f1.readlines()
        with open('{}.filter.output'.format(prefix),'w') as f2:
            for line in lines:
                len = line.split("\t")[2]
                if int(len) >= 60:
                    f2.write(line)

    return 0
