import sys
from tkinter.tix import Tree
from recomb_map import TreeIO
def main():
    ts_addr = sys.argv[1]
    border_list_addr = sys.argv[2]
    output_addr = sys.argv[3]
    maf = 0
    if len(sys.argv) > 4:
        maf = float(sys.argv[4])
    print('Single tree sequence file to single plink binary file.')
    tio = TreeIO()
    print(f'Loading the tree sequence in memory: {ts_addr}')
    tio.read_from_file(ts_addr,border_list_addr)
    print('Dividing the TS objects into chromosomes...')
    tio.chr_divider()
    print(f'Writing variants in this file: {output_addr} - Minimum MAF threshold: {maf}')
    tio.write_bed(output_addr,maf=maf)
    print('Done!')

if __name__ == '__main__':
    main()