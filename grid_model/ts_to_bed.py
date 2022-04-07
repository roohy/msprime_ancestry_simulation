import sys
from tkinter.tix import Tree
from recomb_map import TreeIO
def main():
    ts_addr = sys.argv[1]
    border_list_addr = sys.argv[2]
    output_addr = sys.argv[3]
    if len(sys.argv) > 4:
        maf = float(sys.argv[4])
    tio = TreeIO()

    tio.read_from_file(ts_addr,border_list_addr)
    tio.chr_divider()
    tio.write_bed(output_addr,maf=maf)

if __name__ == '__main__':
    main()